# AIList

This is rust implementation of the AIList algorithm as described
[here](https://www.biorxiv.org/content/10.1101/593657v1). There are some
differences in construction, and I believe some bugs have been removed,
but it's hard to tell. This one passes tests. Perhaps the biggest
difference is that this implementation dynamicaly determines the max
number of components instead of capping at 10. One might call it a
Balance Augmented Interval List ... or BAIList. It the log2 of the input
element lengths to be the max number of components and then decomposes
into that. 

It seems to be very fast. As fast as rust-lapper in all easy cases with
spread our intervals, and faster when things get nested. Benchmarks will
be added as the `interval_bakeoff` project moves along.

## TODOs

- Update docs
- Add seek method and other helper methods?
- Add benchmarks
- Make find a real iterator

## Install

## Usage
