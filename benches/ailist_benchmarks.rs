#[macro_use]
extern crate criterion;
//extern crate AIList;
extern crate rand;

//use criterion::black_box;
use criterion::Criterion;
//use rand::prelude::*;
use rand::Rng;
//use std::time::Duration;
use ailist::{AIList, Interval};

type Iv = (Interval, u32);

fn randomi(imin: u32, imax: u32) -> u32 {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: u32, range_max: u32, size_min: u32, size_max: u32) -> Vec<Iv> {
    let mut result = Vec::with_capacity(n as usize);
    for i in 0..n {
        let s = randomi(0, range_max);
        let e = s + randomi(size_min, size_max);
        result.push((
            Interval {
                start: s,
                stop: e,
                max_end: 0,
            },
            0,
        ));
    }
    result
}

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>) {
    //let n = 3_000_000;
    let n = 10_000;
    let chrom_size = 100_000_000;
    let min_interval_size = 500;
    let max_interval_size = 80000;
    let intervals = make_random(n, chrom_size, min_interval_size, max_interval_size);
    let other_intervals = make_random(n, 10 * chrom_size, 1, 2);
    (intervals, other_intervals)
}

pub fn query(c: &mut Criterion) {
    let (intervals, other_intervals) = make_interval_set();
    let mut bad_intervals = intervals.clone();

    // Make Lapper intervals
    let lapper = AIList::new(intervals, Some(10), None);
    let other_lapper = AIList::new(other_intervals, Some(10), None);
    bad_intervals.push((
        Interval {
            start: 0,
            stop: 90_000_000,
            max_end: 0,
        },
        0,
    ));
    let bad_lapper = AIList::new(bad_intervals, None, None);

    let mut group = c.benchmark_group("AIList find overlaps");
    group.bench_function("find with 100% hit rate", |b| {
        b.iter(|| {
            for (x, _) in lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });

    group.bench_function("find with below 100% hit rate", |b| {
        b.iter(|| {
            for (x, _) in other_lapper.iter() {
                lapper.find(x.start, x.stop).count();
            }
        });
    });

    group.bench_function(
        "find with 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for (x, _) in lapper.iter() {
                    bad_lapper.find(x.start, x.stop).count();
                }
            });
        },
    );
    group.finish();
}

criterion_group!(benches, query);
criterion_main!(benches);
