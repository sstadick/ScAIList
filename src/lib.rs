use std::cmp::Ordering::{self};
#[derive(Debug)]
pub struct AIList<T> {
    /// The list of intervals
    intervals: Vec<Interval>,
    /// The labels that map to each interval
    labels: Vec<T>,
    /// List of starts of sublists
    header: Vec<usize>,
}

/// Hold the start and stop of each sublist
#[derive(Eq, Debug, Clone, Copy)]
pub struct Interval {
    pub start: u32,
    pub stop: u32,
    pub max_end: u32,
}

impl<T: Clone> AIList<T = Clone> {
    /// Create a new AIList out of the passed in intervals. the `min_cov_len`
    pub fn new(mut intervals_labels: Vec<(Interval, T)>, min_cov_len: Option<usize>) -> Self {
        let min_cov_len = min_cov_len.unwrap_or(20);
        intervals_labels.sort_by(|a, b| a.0.cmp(&b.0)); // sort by start site

        let mut header = vec![0];
        let mut intervals = Vec::with_capacity(intervals_labels.len());
        // Decompose long intervals
        while intervals_labels.len() > 0 {
            let (list_1, list_2) = Self::decompose(intervals_labels, min_cov_len);
            intervals.extend(list_1);
            header.push(intervals.len());
            intervals_labels = list_2;
        }

        // Set the max stop on all intervals
        let (intervals, labels): (Vec<Interval>, Vec<T>) = intervals.into_iter().unzip();
        AIList {
            intervals,
            labels,
            header,
        }
    }

    #[inline]
    fn decompose(
        list: Vec<(Interval, T)>,
        min_cov_len: usize,
    ) -> (Vec<(Interval, T)>, Vec<(Interval, T)>) {
        let mut list_1 = vec![];
        let mut list_2 = vec![];

        for i in 0..list.len() {
            let interval = list[i].0;
            let mut covered = 0;
            for j in i..std::cmp::max(i + (2 * min_cov_len), list.len()) {
                if interval.overlap(list[j].0.start, list[j].0.stop) {
                    covered += 1;
                }
            }
            // check if list[i] covers more than the min coverage
            // TODO: make the clones go away
            if covered > min_cov_len {
                list_2.push(list[i].clone());
            } else {
                list_1.push(list[i].clone());
            }
        }
        // Add the running max to list_1 since that one will be sticking around
        Self::add_running_max(&mut list_1);

        (list_1, list_2)
    }

    #[inline]
    fn add_running_max(list: &mut Vec<(Interval, T)>) {
        let mut max_stop = 0;
        for (iv, _label) in list.iter_mut() {
            if iv.stop >= max_stop {
                max_stop = iv.stop;
            }
            iv.stop = max_stop;
        }
    }

    /// Binary search to find the right most index where interval.start < query.stop
    /// TODO: test this
    #[inline]
    fn upper_bound(stop: u32, intervals: &[Interval]) -> usize {
        let mut size = intervals.len();
        let mut high = size;
        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = high - half;
            let other_high = high - other_half;
            let v = &intervals[probe];
            size = half;
            high = if v.start >= stop { other_high } else { high }
        }
        high
    }

    #[inline]
    fn find(&self, start: u32, stop: u32) -> IterFind<T> {
        IterFind {
            inner: self,
            off: Self::lower_bound(
                stop,
                &self.intervals[..self.header.get(1).unwrap_or(self.intervals.len())],
            ),
            header_offset: 0,
            end: self.intervals.len(),
            start,
            stop,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, T> {
    inner: &'a AIList<T>,
    off: usize,
    header_offset: usize,
    end: usize,
    start: u32,
    stop: u32,
}

impl<'a, T> Iterator for IterFind<'a, T> {
    type Item = &'a IntervalNode<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(h) = self.inner.header.get(self.header_offset) {
            self.off = AIList::upper_bound(
                self.stop,
                &self.inner.intervals[h..self
                    .header
                    .get(self.header_offset + 1)
                    .unwrap_or(self.inner.intervals.len())],
            );
            // TODO: Stopped here, not sure this is the best way to iterate
        }
        while self.off < self.inner.intervals.len() {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.stop) {
                return self.inner.interval_nodes.get(self.off - 1);
            } else if interval.start >= self.stop {
                break;
            }
        }
        None
    }
}
impl Interval {
    /// Create a new interval with default max_end
    pub fn new(start: u32, stop: u32) -> Self {
        Interval {
            start,
            stop,
            max_end: stop,
        }
    }
    /// Compute the intersect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Interval) -> u32 {
        std::cmp::min(self.stop, other.stop)
            .checked_sub(std::cmp::max(self.start, other.start))
            .unwrap_or(0)
    }

    /// Internal version of intersect for working with IntervalNode / Interval
    pub fn intersect_raw(&self, start: u32, stop: u32) -> u32 {
        std::cmp::min(self.stop, stop)
            .checked_sub(std::cmp::max(self.start, start))
            .unwrap_or(0)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: u32, stop: u32) -> bool {
        self.start < stop && self.stop > start
    }
}

impl Ord for Interval {
    #[inline]
    fn cmp(&self, other: &Interval) -> Ordering {
        if self.start < other.start {
            Ordering::Less
        } else if other.start < self.start {
            Ordering::Greater
        } else {
            self.stop.cmp(&other.stop)
        }
    }
}

impl PartialOrd for Interval {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl PartialEq for Interval {
    #[inline]
    fn eq(&self, other: &Interval) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
