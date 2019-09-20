#[macro_use]
extern crate criterion;
//extern crate ScAIList;
extern crate rand;

//use criterion::black_box;
use criterion::Criterion;
//use rand::prelude::*;
use rand::Rng;
//use std::time::Duration;
use scailist::{Interval, ScAIList};

type Iv = Interval<u32>;

fn randomi(imin: u32, imax: u32) -> u32 {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: u32, range_max: u32, size_min: u32, size_max: u32) -> Vec<Iv> {
    let mut result = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let s = randomi(0, range_max);
        let e = s + randomi(size_min, size_max);
        result.push(Interval {
            start: s,
            end: e,
            val: 0,
        });
    }
    result
}

fn make_interval_set() -> (Vec<Iv>, Vec<Iv>) {
    //let n = 3_000_000;
    let n = 1_000;
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
    let ivs = intervals.clone();
    let ovs = other_intervals.clone();
    let lapper = ScAIList::new(intervals, None);
    bad_intervals.push(Interval {
        start: 0,
        end: 90_000_000,
        val: 0,
    });
    let bad_lapper = ScAIList::new(bad_intervals, None);

    let mut group = c.benchmark_group("ScAIList find overlaps");
    group.bench_function("find with 100% hit rate", |b| {
        b.iter(|| {
            for x in ivs.iter() {
                lapper.find(x.start, x.end).count();
            }
        });
    });

    group.bench_function("find with below 100% hit rate", |b| {
        b.iter(|| {
            for x in ovs.iter() {
                lapper.find(x.start, x.end).count();
            }
        });
    });

    group.bench_function(
        "find with 100% hit rate - chromosome spanning interval",
        |b| {
            b.iter(|| {
                for x in ivs.iter() {
                    bad_lapper.find(x.start, x.end).count();
                }
            });
        },
    );
    group.finish();
}

criterion_group!(benches, query);
criterion_main!(benches);
