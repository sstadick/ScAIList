![docs](https://docs.rs/scailist/badge.svg)
![crates.io](https://img.shields.io/crates/v/scailist.svg)

# ScAIList

This is rust implementation of the AIList algorithm as described
[here](https://www.biorxiv.org/content/10.1101/593657v1). The biggest
difference is that this implementation dynamicaly determines the max
number of components instead of capping at 10. One might call it a
Scaled Augmented Interval List. It takes the log2 of the input
element lengths to be the max number of components and then decomposes
into that. 

It seems to be very fast. As fast as rust-lapper in all easy cases with
spread our intervals, and faster when things get nested. Benchmarks will
be added as the `interval_bakeoff` project moves along.

[Documentation](https://docs.rs/scailist)

[Crates.io](https://crates.io/crates/scailist)
