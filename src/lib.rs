use std::cmp::Ordering::{self, Equal, Greater, Less};
#[derive(Debug)]
pub struct AIList<T> {
    /// The list of intervals
    intervals: Vec<Interval>,
    /// The labels that map to each interval
    labels: Vec<T>,
    /// List of starts of sublists
    headers: Vec<Header>,
    /// The max number of items to just do a linear search on
    max_lin_search: usize,
}

/// Hold the start and stop of each sublist
#[derive(Eq, Debug, Clone, Copy)]
pub struct Interval {
    pub start: u32,
    pub stop: u32,
    pub max_end: u32,
}

#[derive(Debug)]
struct Header {
    start: usize,
    stop: usize,
}

impl<T: Clone + std::fmt::Debug> AIList<T> {
    /// Create a new AIList out of the passed in intervals. the `min_cov_len`
    pub fn new(
        mut intervals_labels: Vec<(Interval, T)>,
        min_cov_len: Option<usize>,
        max_lin_search: Option<usize>,
    ) -> Self {
        let min_cov_len = min_cov_len.unwrap_or(20);
        let max_lin_search = max_lin_search.unwrap_or(15);
        intervals_labels.sort_by(|a, b| a.0.cmp(&b.0)); // sort by start site

        let mut headers = vec![];
        let mut intervals = Vec::with_capacity(intervals_labels.len());
        // Decompose long intervals
        let mut header_start = 0;
        while intervals_labels.len() > 0 {
            let (list_1, list_2) = Self::decompose(intervals_labels, min_cov_len);
            intervals.extend(list_1);
            headers.push(Header {
                start: header_start,
                stop: intervals.len(),
            });
            header_start = intervals.len();
            intervals_labels = list_2;
        }

        let (intervals, labels): (Vec<Interval>, Vec<T>) = intervals.into_iter().unzip();
        AIList {
            intervals,
            labels,
            headers,
            max_lin_search,
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
            let mut j = i;
            let j_end = std::cmp::min(i + (2 * min_cov_len), list.len());
            while j < j_end && interval.stop < list[j].0.start {
                if interval.start < list[j].0.stop {
                    covered += 1;
                }
                j += 1;
            }
            //for j in i..std::cmp::min(i + (2 * min_cov_len), list.len()) {
            //if interval.overlap(list[j].0.start, list[j].0.stop) {
            //covered += 1;
            //}
            //}
            // check if list[i] covers more than the min coverage
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
            iv.max_end = max_stop;
        }
    }

    /// Binary search to find the right most index where interval.start < query.stop
    #[inline]
    pub fn upper_bound(stop: u32, intervals: &[Interval]) -> Option<usize> {
        let mut right = intervals.len() - 1;
        let mut left = 0;

        if intervals[right - 1].start < stop {
            // last start pos is less than the stop, then return the last pos
            return Some(right);
        } else if intervals[left].start >= stop {
            // first start pos > stop, not in this cluster at all
            return None;
        }

        while left < right - 1 {
            let mid = (left + right) / 2;
            if intervals[mid].start >= stop {
                right = mid - 1;
            } else {
                left = mid;
            }
        }
        if intervals[right].start < stop {
            Some(right)
        } else if intervals[left].start < stop {
            Some(left)
        } else {
            None
        }
    }

    #[inline]
    pub fn iter(&self) -> IterAIList<T> {
        IterAIList {
            inner: self,
            pos: 0,
        }
    }

    #[inline]
    pub fn find(&self, start: u32, stop: u32) -> IterFind<T> {
        let mut result = vec![];
        for h_idx in 0..self.headers.len() {
            let header = &self.headers[h_idx];
            if header.stop - header.start > 15 {
                let mut offset =
                    match Self::upper_bound(stop, &self.intervals[header.start..header.stop]) {
                        Some(n) => n,
                        None => continue,
                    };
                while offset >= header.start && self.intervals[offset].max_end > start {
                    if self.intervals[offset].stop > start {
                        result.push((&self.intervals[offset], &self.labels[offset]));
                    }
                    offset = match offset.checked_sub(1) {
                        Some(n) => n,
                        None => break,
                    }
                }
            } else {
                for offset in header.start..header.stop {
                    if self.intervals[offset].start < stop && self.intervals[offset].stop > start {
                        result.push((&self.intervals[offset], &self.labels[offset]));
                    }
                }
            }
        }
        // Gather a list of starting points for each header
        IterFind { curr: 0, result }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, T> {
    curr: usize,
    result: Vec<(&'a Interval, &'a T)>,
}

impl<'a, T: Clone + std::fmt::Debug> Iterator for IterFind<'a, T> {
    type Item = (&'a Interval, &'a T);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.result.len() {
            self.curr += 1;
            return Some(self.result[self.curr - 1]);
        } else {
            return None;
        }
    }
}

/// AIList Iterator
pub struct IterAIList<'a, T>
where
    T: Clone + std::fmt::Debug + 'a,
{
    inner: &'a AIList<T>,
    pos: usize,
}

impl<'a, T: Clone + std::fmt::Debug> Iterator for IterAIList<'a, T> {
    type Item = (&'a Interval, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            Some((
                &self.inner.intervals[self.pos - 1],
                &self.inner.labels[self.pos - 1],
            ))
        }
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
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval;
    fn setup_nonoverlapping() -> AIList<u32> {
        let data: Vec<(Iv, u32)> = (0..100)
            .step_by(20)
            .map(|x| (Iv {
                start: x,
                stop: x + 10,
                max_end: 0,
            }, 0))
            .collect();
        let ailist= AIList::new(data, Some(3), Some(5));
       ailist 
    }
    fn setup_overlapping() -> AIList<u32> {
        let data: Vec<(Iv, u32)> = (0..100)
            .step_by(10)
            .map(|x| (Iv {
                start: x,
                stop: x + 15,
                max_end: 0,
            }, 0))
            .collect();
        let ailist= AIList::new(data, Some(1), Some(5));
       ailist 
    }
    fn setup_badlapper() -> AIList<u32> {
        let data: Vec<(Iv, u32)> = vec![
            (Iv{start: 70, stop: 120, max_end: 0}, 0), // max_len = 50
            (Iv{start: 10, stop: 15, max_end: 0}, 0),
            (Iv{start: 10, stop: 15, max_end: 0}, 0), // exact overlap
            (Iv{start: 12, stop: 15, max_end: 0}, 0), // inner overlap
            (Iv{start: 14, stop: 16, max_end: 0}, 0), // overlap end
            (Iv{start: 40, stop: 45, max_end: 0}, 0),
            (Iv{start: 50, stop: 55, max_end: 0}, 0),
            (Iv{start: 60, stop: 65, max_end: 0}, 0),
            (Iv{start: 68, stop: 71, max_end: 0}, 0), // overlap start
            (Iv{start: 70, stop: 75, max_end: 0}, 0),
        ];
        let ailist= AIList::new(data, Some(3), Some(5));
       ailist 
    }
    fn setup_single() -> AIList<u32> {
        let data: Vec<(Iv, u32)> = vec![(Iv {
            start: 10,
            stop: 35,
            max_end: 0,
        }, 0)];
        let ailist= AIList::new(data, Some(3), Some(5));
       ailist 
    }

    // Test that a query stop that hits an interval start returns no interval
    #[test]
    fn test_query_stop_interval_start() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        //assert_eq!(None, lapper.seek(15, 20, &mut cursor).next())
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_stop() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        //assert_eq!(None, lapper.seek(30, 35, &mut cursor).next())
    }

    // Test that a query that overlaps the start of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_start() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = (&Iv {
            start: 20,
            stop: 30,
            max_end: 0,
        }, &0);
        assert_eq!(Some(expected), lapper.find(15, 25).next());
        //assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next())
    }

    // Test that a query that overlaps the stop of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_stop() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = (&Iv {
            start: 20,
            stop: 30,
            max_end: 0,
        }, &0);
        assert_eq!(Some(expected), lapper.find(25, 35).next());
        //assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next())
    }

    // Test that a query that is enveloped by interval returns interval
    #[test]
    fn test_interval_envelops_query() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = (&Iv {
            start: 20,
            stop: 30,
            max_end: 0,
        }, &0);
        assert_eq!(Some(expected), lapper.find(22, 27).next());
        //assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next())
    }

    // Test that a query that envolops an interval returns that interval
    #[test]
    fn test_query_envolops_interval() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = (&Iv {
            start: 20,
            stop: 30,
            max_end: 0,
        }, &0);
        assert_eq!(Some(expected), lapper.find(15, 35).next());
        //assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next())
    }

    #[test]
    fn test_overlapping_intervals() {
        let lapper = setup_overlapping();
        //let mut cursor = 0;
        let e1 = (&Iv {
            start: 0,
            stop: 15,
            max_end: 0,
        }, &0);
        let e2 = (&Iv {
            start: 10,
            stop: 25,
            max_end: 0,
        }, &0);
        let mut found = lapper.find(8, 20).collect::<Vec<(&Iv,&u32)>>();
        found.sort_by(|a, b| a.0.start.cmp(&b.0.start));
        assert_eq!(vec![e1, e2], found);
        //assert_eq!(
            //vec![&e1, &e2],
            //lapper.seek(8, 20, &mut cursor).collect::<Vec<&Iv>>()
        //);
    }
}
