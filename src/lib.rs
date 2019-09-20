//! This module provides an implementation of an AIList, but with a dynamic scaling for the number
//! of sublists.
//! ## Features
//! - Extremely consistant. The way the input intervals are decomposed diminishes the effects of
//! super containment.
//! - Parallel friendly. Queries are on an immutable structure, even for seek
//! - Consumer / Adapter paradigm, an iterator is returned.
//!
//! ## Details:
//!
//! Please see the [paper](https://www.biorxiv.org/content/10.1101/593657v1).
//!
//! Most interaction with this crate will be through the [`ScAIList``](struct.ScAIList.html)` struct
//! The main methods is [`find`](struct.Lapper.html#method.find).
//!
//! The overlap function for this assumes a zero based genomic coordinate system. So [start, stop)
//! is not inclusive of the stop position for neither the queries, nor the Intervals.
//!
//! ScAIList is composed of four primary parts. A main interval list, which holds all the
//! intervals after they have been decomposed. A component index's list, which holds the start
//! index of each sublist post-decomposition, A component lengths list, which holds the length of
//! each component, and finally a max_ends list, which holds the max end releative to a sublist up
//! to a given point for each interval.
//!
//! The decomposition step is achieved by walking the list of intervals and recursively (with a
//! cap) extracting intervals that overlap a given number of other intervals within a certain
//! distance from it. The unique development in this implementation is to make the cap dynamic.
//!
//! # Examples
//!
//! ```rust
//!    use scailist::{Interval, ScAIList};
//!    use std::cmp;
//!    type Iv = Interval<u32>;
//!
//!    // create some fake data
//!    let data: Vec<Iv> = (0..20).step_by(5).map(|x| Iv{start: x, end: x + 2, val: 0}).collect();
//!    println!("{:#?}", data);
//!
//!    // make lapper structure
//!    let laps = ScAIList::new(data, None);
//!    assert_eq!(laps.find(6, 11).next(), Some(&Iv{start: 5, end: 7, val: 0}));
//!    
//!    let mut sim: u32= 0;
//!    // Calculate the overlap between the query and the found intervals, sum total overlap
//!    for i in (0..10).step_by(3) {
//!        sim += laps
//!            .find(i, i + 2)
//!            .map(|iv| cmp::min(i + 2, iv.end) - cmp::max(i, iv.start))
//!            .sum::<u32>();
//!    }
//!    assert_eq!(sim, 4);
//! ```
use std::cmp::Ordering;

/// This is the main object of this repo, see associated methods
#[derive(Debug)]
pub struct ScAIList<T: Clone + Eq + std::fmt::Debug> {
    /// The list of intervals
    intervals: Vec<Interval<T>>,
    // number of comps total
    num_comps: usize,
    // list of lengths of comps
    comp_lens: Vec<usize>,
    // start index's of comps
    comp_idxs: Vec<usize>,
    // vec of max ends, pairs with intervals
    max_ends: Vec<u32>,
}

/// Hold the start and stop of each sublist
#[derive(Eq, Debug, Clone, Copy)]
pub struct Interval<T: Clone + Eq + std::fmt::Debug> {
    pub start: u32,
    pub end: u32,
    pub val: T,
}

/// The ScAIList itself
impl<T: Clone + Eq + std::fmt::Debug> ScAIList<T> {
    /// Create a new ScAIList out of the passed in intervals. The min_cov_len should probably be
    /// left as default, which is 20. It dictates how far ahead to look from a given point to
    /// determine if that interval covers enough other intervals to be moved to a sublist. The
    /// number of intervals it has to cover is equal to min_cov_len / 2. The number of sublists
    /// that might be fored is capped at intervals.len().log2(), but if there aren't many overlaps,
    /// fewer will be created.
    pub fn new(mut input_intervals: Vec<Interval<T>>, min_cov_len: Option<usize>) -> Self {
        let max_comps = (input_intervals.len() as f64).log2().floor() as usize + 1;
        let min_cov_len = min_cov_len.unwrap_or(20); // number of elements ahead to check for cov
        let min_cov = min_cov_len / 2; // the number of elemnts covered to trigger an extraction

        let num_comps; // number of components
        let mut comp_lens = vec![]; // lengths of each component
        let mut comp_idxs = vec![]; // start positions of each component
        let mut max_ends = vec![]; // list of the max end positions

        let min_comp_len = std::cmp::max(64, min_cov_len); // min length of a component

        let input_len = input_intervals.len();
        input_intervals.sort();

        let mut decomposed = vec![];
        if input_len <= min_comp_len {
            num_comps = 1;
            comp_lens.push(input_len);
            comp_idxs.push(0);
            decomposed.append(&mut input_intervals);
        } else {
            // TODO: I suspect that in here is where the decomp goes badly and we basically end up
            // with a glorified linear search
            let mut curr_comp = 0;
            while curr_comp < max_comps && input_len - decomposed.len() > min_comp_len {
                let mut list1 = vec![];
                let mut list2 = vec![];
                for i in 0..input_intervals.len() {
                    let interval = &input_intervals[i];
                    let mut j = 1;
                    let mut cov = 1;
                    while j < min_comp_len && cov < min_cov && j + i < input_intervals.len() {
                        if input_intervals[j + i].end >= interval.end {
                            cov += 1;
                        }
                        j += 1;
                    }
                    if cov < min_cov {
                        list1.push(input_intervals[i].clone());
                    } else {
                        list2.push(input_intervals[i].clone())
                    }
                }
                // Add the component info to ScAIList
                comp_idxs.push(decomposed.len());
                comp_lens.push(list1.len());
                curr_comp += 1;

                if list2.len() <= min_comp_len || curr_comp == max_comps - 2 {
                    // exit: add L2 to the end
                    if list2.is_empty() {
                        decomposed.append(&mut list1);
                        comp_idxs.push(decomposed.len());
                        comp_lens.push(list2.len());
                        decomposed.append(&mut list2);
                        curr_comp += 1;
                    }
                } else {
                    decomposed.append(&mut list1);
                    input_intervals = list2;
                }
            }
            num_comps = curr_comp;
        }

        // Augment with maxend
        for j in 0..num_comps {
            let comp_start = comp_idxs[j];
            let comp_end = comp_start + comp_lens[j];
            let mut max_end = decomposed[comp_start].end;
            max_ends.push(max_end);
            for iv in decomposed[comp_start + 1..comp_end].iter() {
                if iv.end > max_end {
                    max_end = iv.end;
                }
                max_ends.push(max_end);
            }
        }

        ScAIList {
            num_comps,
            comp_idxs,
            comp_lens,
            max_ends,
            intervals: decomposed,
        }
    }

    /// Binary search to find the right most index where interval.start < query.stop
    #[inline]
    pub fn upper_bound(stop: u32, intervals: &[Interval<T>]) -> Option<usize> {
        let mut right = intervals.len();
        let mut left = 0;

        if intervals[right - 1].start < stop {
            // last start pos is less than the stop, then return the last pos
            return Some(right - 1);
        } else if intervals[left].start >= stop {
            // first start pos > stop, not in this cluster at all
            return None;
        }

        while right > 0 {
            let half = right / 2;
            let other_half = right - half;
            let probe = left + half;
            let other_left = left + other_half;
            let v = &intervals[probe];
            right = half;
            left = if v.start < stop { other_left } else { left }
        }
        // Guarded at the top from ending on either extreme
        if intervals[left].start >= stop {
            Some(left - 1)
        } else {
            Some(left)
        }
    }

    #[inline]
    pub fn iter(&self) -> IterScAIList<T> {
        IterScAIList {
            inner: self,
            pos: 0,
        }
    }

    #[inline]
    // TODO: Convert to a real iterator
    // TODO: Add seek... multiple cursors? how would that work?
    pub fn find(&self, start: u32, stop: u32) -> IterFind<T> {
        let mut result = vec![];
        for comp_num in 0..self.num_comps {
            let comp_start = self.comp_idxs[comp_num];
            let comp_end = comp_start + self.comp_lens[comp_num];
            if self.comp_lens[comp_num] > 15 {
                let mut offset =
                    match Self::upper_bound(stop, &self.intervals[comp_start..comp_end]) {
                        Some(n) => n,
                        None => continue,
                    };
                offset += comp_start;
                while offset >= comp_start && self.max_ends[offset] > start {
                    if self.intervals[offset].end > start {
                        result.push(&self.intervals[offset]);
                    }
                    offset = match offset.checked_sub(1) {
                        Some(n) => n,
                        None => break,
                    };
                }
            } else {
                for offset in comp_start..comp_end {
                    if self.intervals[offset].start < stop && self.intervals[offset].end > start {
                        result.push(&self.intervals[offset]);
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
pub struct IterFind<'a, T: Clone + Eq + std::fmt::Debug> {
    curr: usize,
    result: Vec<&'a Interval<T>>,
}

impl<'a, T: Clone + Eq + std::fmt::Debug> Iterator for IterFind<'a, T> {
    type Item = &'a Interval<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.result.len() {
            self.curr += 1;
            Some(self.result[self.curr - 1])
        } else {
            None
        }
    }
}

/// ScAIList Iterator
pub struct IterScAIList<'a, T>
where
    T: Clone + Eq + std::fmt::Debug + 'a,
{
    inner: &'a ScAIList<T>,
    pos: usize,
}

impl<'a, T: Clone + Eq + std::fmt::Debug> Iterator for IterScAIList<'a, T> {
    type Item = &'a Interval<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            Some(&self.inner.intervals[self.pos - 1])
        }
    }
}

impl<T: Clone + Eq + std::fmt::Debug> Interval<T> {
    /// Compute the intersect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Self) -> u32 {
        std::cmp::min(self.end, other.end)
            .checked_sub(std::cmp::max(self.start, other.start))
            .unwrap_or(0)
    }

    /// Internal version of intersect for working with IntervalNode / Interval
    pub fn intersect_raw(&self, start: u32, end: u32) -> u32 {
        std::cmp::min(self.end, end)
            .checked_sub(std::cmp::max(self.start, start))
            .unwrap_or(0)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: u32, end: u32) -> bool {
        self.start < end && self.end > start
    }
}

impl<T: Clone + Eq + std::fmt::Debug> Ord for Interval<T> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        if self.start < other.start {
            Ordering::Less
        } else if other.start < self.start {
            Ordering::Greater
        } else {
            self.end.cmp(&other.end)
        }
    }
}

impl<T: Clone + Eq + std::fmt::Debug> PartialOrd for Interval<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl<T: Clone + Eq + std::fmt::Debug> PartialEq for Interval<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval<u32>;
    fn setup_nonoverlapping() -> ScAIList<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                end: x + 10,
                val: 0,
            })
            .collect();
        let lapper = ScAIList::new(data, None);
        lapper
    }
    fn setup_overlapping() -> ScAIList<u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                end: x + 15,
                val: 0,
            })
            .collect();
        let lapper = ScAIList::new(data, None);
        lapper
    }
    fn setup_badlapper() -> ScAIList<u32> {
        let data: Vec<Iv> = vec![
            Iv{start: 70, end: 120, val: 0}, // max_len = 50
            Iv{start: 10, end: 15, val: 0},
            Iv{start: 10, end: 15, val: 0}, // exact overlap
            Iv{start: 12, end: 15, val: 0}, // inner overlap
            Iv{start: 14, end: 16, val: 0}, // overlap end
            Iv{start: 40, end: 45, val: 0},
            Iv{start: 50, end: 55, val: 0},
            Iv{start: 60, end: 65, val: 0},
            Iv{start: 68, end: 71, val: 0}, // overlap start
            Iv{start: 70, end: 75, val: 0},
        ];
        let lapper = ScAIList::new(data, None);
        lapper
    }
    fn setup_single() -> ScAIList<u32> {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            end: 35,
            val: 0,
        }];
        let lapper = ScAIList::new(data, None);
        lapper
    }

    // Test that a query end that hits an interval start returns no interval
    #[test]
    fn test_query_end_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        //assert_eq!(None, lapper.seek(15, 20, &mut cursor).next())
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_end() {
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
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 25).next());
        //assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next())
    }

    // Test that a query that overlaps the end of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_end() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(25, 35).next());
        //assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next())
    }

    // Test that a query that is enveloped by interval returns interval
    #[test]
    fn test_interval_envelops_query() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(22, 27).next());
        //assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next())
    }

    // Test that a query that envolops an interval returns that interval
    #[test]
    fn test_query_envolops_interval() {
        let lapper = setup_nonoverlapping();
        //let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 35).next());
        //assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next())
    }

    #[test]
    fn test_overlapping_intervals() {
        let lapper = setup_overlapping();
        //let mut cursor = 0;
        let e1 = Iv {
            start: 0,
            end: 15,
            val: 0,
        };
        let e2 = Iv {
            start: 10,
            end: 25,
            val: 0,
        };
        assert_eq!(vec![&e1, &e2], lapper.find(8, 20).collect::<Vec<&Iv>>());
        //assert_eq!(
            //vec![&e1, &e2],
            //lapper.seek(8, 20, &mut cursor).collect::<Vec<&Iv>>()
        //);
    }

    //#[test]
    //fn test_merge_overlaps() {
        //let mut lapper = setup_badlapper();
        //let expected: Vec<&Iv> = vec![
            //&Iv{start: 10, end: 16, val: 0},
            //&Iv{start: 40, end: 45, val: 0},
            //&Iv{start: 50, end: 55, val: 0},
            //&Iv{start: 60, end: 65, val: 0},
            //&Iv{start: 68, end: 120, val: 0}, // max_len = 50
        //];
        //lapper.merge_overlaps();
        //assert_eq!(expected, lapper.iter().collect::<Vec<&Iv>>())
        
    //}

    //#[test]
    //fn test_lapper_cov() {
        //let mut lapper = setup_badlapper();
        //let before = lapper.cov();
        //lapper.merge_overlaps();
        //let after = lapper.cov();
        //assert_eq!(before, after);

        //let mut lapper = setup_nonoverlapping();
        //lapper.set_cov();
        //assert_eq!(lapper.cov(), 50);
    //}

    #[test]
    fn test_interval_intersects() {
        let i1 = Iv{start: 70, end: 120, val: 0}; // max_len = 50
        let i2 = Iv{start: 10, end: 15, val: 0};
        let i3 = Iv{start: 10, end: 15, val: 0}; // exact overlap
        let i4 = Iv{start: 12, end: 15, val: 0}; // inner overlap
        let i5 = Iv{start: 14, end: 16, val: 0}; // overlap end
        let i6 = Iv{start: 40, end: 50, val: 0};
        let i7 = Iv{start: 50, end: 55, val: 0};
        let i_8 = Iv{start: 60, end: 65, val: 0};
        let i9 = Iv{start: 68, end: 71, val: 0}; // overlap start
        let i10 = Iv{start: 70, end: 75, val: 0};

        assert_eq!(i2.intersect(&i3), 5); // exact match
        assert_eq!(i2.intersect(&i4), 3); // inner intersect
        assert_eq!(i2.intersect(&i5), 1); // end intersect
        assert_eq!(i9.intersect(&i10), 1); // start intersect
        assert_eq!(i7.intersect(&i_8), 0); // no intersect
        assert_eq!(i6.intersect(&i7), 0); // no intersect end = start
        assert_eq!(i1.intersect(&i10), 5); // inner intersect at start
    }

    //#[test]
    //fn test_union_and_intersect() {
        //let data1: Vec<Iv> = vec![
            //Iv{start: 70, end: 120, val: 0}, // max_len = 50
            //Iv{start: 10, end: 15, val: 0}, // exact overlap
            //Iv{start: 12, end: 15, val: 0}, // inner overlap
            //Iv{start: 14, end: 16, val: 0}, // overlap end
            //Iv{start: 68, end: 71, val: 0}, // overlap start
        //];
        //let data2: Vec<Iv> = vec![

            //Iv{start: 10, end: 15, val: 0},
            //Iv{start: 40, end: 45, val: 0},
            //Iv{start: 50, end: 55, val: 0},
            //Iv{start: 60, end: 65, val: 0},
            //Iv{start: 70, end: 75, val: 0},
        //];
        
        //let (mut lapper1, mut lapper2) = (ScAIList::new(data1), ScAIList::new(data2)) ;
        //// Should be the same either way it's calculated
        //let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        //assert_eq!(intersect, 10);
        //assert_eq!(union, 73);
        //let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        //assert_eq!(intersect, 10);
        //assert_eq!(union, 73);
        //lapper1.merge_overlaps();
        //lapper1.set_cov();
        //lapper2.merge_overlaps();
        //lapper2.set_cov();

        //// Should be the same either way it's calculated
        //let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        //assert_eq!(intersect, 10);
        //assert_eq!(union, 73);
        //let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        //assert_eq!(intersect, 10);
        //assert_eq!(union, 73);
    //}

    #[test]
    fn test_find_overlaps_in_large_intervals() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, end: 8, val: 0},
            Iv{start: 1, end: 10, val: 0}, 
            Iv{start: 2, end: 5, val: 0}, 
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 4, end: 7, val: 0},
            Iv{start: 5, end: 8, val: 0},
            Iv{start: 8, end: 8, val: 0},
            Iv{start: 9, end: 11, val: 0},
            Iv{start: 10, end: 13, val: 0},
            Iv{start: 100, end: 200, val: 0},
            Iv{start: 110, end: 120, val: 0},
            Iv{start: 110, end: 124, val: 0},
            Iv{start: 111, end: 160, val: 0},
            Iv{start: 150, end: 200, val: 0},
        ];
        let lapper = ScAIList::new(data1, None);
        let found = lapper.find(8, 11).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 1, end: 10, val: 0}, 
            &Iv{start: 9, end: 11, val: 0},
            &Iv{start: 10, end: 13, val: 0},
        ]);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 100, end: 200, val: 0},
            &Iv{start: 111, end: 160, val: 0},
            &Iv{start: 150, end: 200, val: 0},
        ]);

    }

    //#[test]
    //fn test_depth_sanity() {
        //let data1: Vec<Iv> = vec![
            //Iv{start: 0, end: 10, val: 0},
            //Iv{start: 5, end: 10, val: 0}
        //];
        //let lapper = ScAIList::new(data1);
        //let found = lapper.depth().collect::<Vec<Interval<u32>>>();
        //assert_eq!(found, vec![
                   //Interval{start: 0, end: 5, val: 1},
                   //Interval{start: 5, end: 10, val: 2}
        //]);
    //}

    //#[test]
    //fn test_depth_hard() {
        //let data1: Vec<Iv> = vec![
            //Iv{start: 1, end: 10, val: 0},
            //Iv{start: 2, end: 5, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 5, end: 8, val: 0},
            //Iv{start: 9, end: 11, val: 0},
        //];
        //let lapper = ScAIList::new(data1);
        //let found = lapper.depth().collect::<Vec<Interval<u32>>>();
        //assert_eq!(found, vec![
                   //Interval{start: 1, end: 2, val: 1},
                   //Interval{start: 2, end: 3, val: 2},
                   //Interval{start: 3, end: 8, val: 5},
                   //Interval{start: 8, end: 9, val: 1},
                   //Interval{start: 9, end: 10, val: 2},
                   //Interval{start: 10, end: 11, val: 1},
        //]);
    //}
    //#[test]
    //fn test_depth_harder() {
        //let data1: Vec<Iv> = vec![
            //Iv{start: 1, end: 10, val: 0},
            //Iv{start: 2, end: 5, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 3, end: 8, val: 0},
            //Iv{start: 5, end: 8, val: 0},
            //Iv{start: 9, end: 11, val: 0},
            //Iv{start: 15, end: 20, val: 0},
        //];
        //let lapper = ScAIList::new(data1);
        //let found = lapper.depth().collect::<Vec<Interval<u32>>>();
        //assert_eq!(found, vec![
                   //Interval{start: 1, end: 2, val: 1},
                   //Interval{start: 2, end: 3, val: 2},
                   //Interval{start: 3, end: 8, val: 5},
                   //Interval{start: 8, end: 9, val: 1},
                   //Interval{start: 9, end: 10, val: 2},
                   //Interval{start: 10, end: 11, val: 1},
                   //Interval{start: 15, end: 20, val: 1},
        //]);
    //}
    // BUG TESTS - these are tests that came from real life

    // Test that it's not possible to induce index out of bounds by pushing the cursor past the end
    // of the lapper.
    //#[test]
    //fn test_seek_over_len() {
        //let lapper = setup_nonoverlapping();
        //let single = setup_single();
        //let mut cursor: usize = 0;

        //for interval in lapper.iter() {
            //for o_interval in single.seek(interval.start, interval.end, &mut cursor) {
                //println!("{:#?}", o_interval);
            //}
        //}
    //}

    // Test that if lower_bound puts us before the first match, we still return a match
    #[test]
    fn test_find_over_behind_first_match() {
        let lapper = setup_badlapper();
        let e1 = Iv {start: 50, end: 55, val: 0};
        let found = lapper.find(50, 55).next();
        assert_eq!(found, Some(&e1));
    }

    // When there is a very long interval that spans many little intervals, test that the little
    // intevals still get returne properly
    #[test]
    fn test_bad_skips() {
        let data = vec![
            Iv{start:25264912, end: 25264986, val: 0},	
            Iv{start:27273024, end: 27273065	, val: 0},
            Iv{start:27440273, end: 27440318	, val: 0},
            Iv{start:27488033, end: 27488125	, val: 0},
            Iv{start:27938410, end: 27938470	, val: 0},
            Iv{start:27959118, end: 27959171	, val: 0},
            Iv{start:28866309, end: 33141404	, val: 0},
        ];
        let lapper = ScAIList::new(data, None);

        let found = lapper.find(28974798, 33141355).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start:28866309, end: 33141404	, val: 0},
        ])

    }
}
