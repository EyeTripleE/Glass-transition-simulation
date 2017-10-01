#!/bin/sh

grep "Rank" timings.tiles_cutoff.hybrid.1024 > tiles_cutoff.hybrid.1024
grep "Rank" timings.tiles_cutoff.hybrid.2048 > tiles_cutoff.hybrid.2048
grep "Rank" timings.tiles_cutoff.hybrid.4096 > tiles_cutoff.hybrid.4096
grep "Rank" timings.tiles_cutoff.hybrid.8192 > tiles_cutoff.hybrid.8192

grep "Rank" timings.tiles.hybrid.1024 > tiles.hybrid.1024
grep "Rank" timings.tiles.hybrid.2048 > tiles.hybrid.2048
grep "Rank" timings.tiles.hybrid.4096 > tiles.hybrid.4096
grep "Rank" timings.tiles.hybrid.8192 > tiles.hybrid.8192

grep "Rank" timings.strips_cutoff.hybrid.8192 > strips_cutoff.hybrid.8192
grep "Rank" timings.strips_cutoff.hybrid.4096 > strips_cutoff.hybrid.4096
grep "Rank" timings.strips_cutoff.hybrid.2048 > strips_cutoff.hybrid.2048
grep "Rank" timings.strips_cutoff.hybrid.1024 > strips_cutoff.hybrid.1024

grep "Rank" timings.strips.hybrid.8192 > strips.hybrid.8192
grep "Rank" timings.strips.hybrid.4096 > strips.hybrid.4096
grep "Rank" timings.strips.hybrid.2048 > strips.hybrid.2048
grep "Rank" timings.strips.hybrid.1024 > strips.hybrid.1024

grep "Rank" timings.bh.hybrid.8192 > bh.hybrid.8192
grep "Rank" timings.bh.hybrid.4096 > bh.hybrid.4096
grep "Rank" timings.bh.hybrid.2048 > bh.hybrid.2048
grep "Rank" timings.bh.hybrid.1024 > bh.hybrid.1024

grep "Rank" timings.bh1.hybrid.8192 > bh1.hybrid.8192
grep "Rank" timings.bh1.hybrid.4096 > bh1.hybrid.4096
grep "Rank" timings.bh1.hybrid.2048 > bh1.hybrid.2048
grep "Rank" timings.bh1.hybrid.1024 > bh1.hybrid.1024

grep "Rank" timings.bh3.hybrid.8192 > bh3.hybrid.8192
grep "Rank" timings.bh3.hybrid.4096 > bh3.hybrid.4096
grep "Rank" timings.bh3.hybrid.2048 > bh3.hybrid.2048
grep "Rank" timings.bh3.hybrid.1024 > bh3.hybrid.1024

grep "Rank" timings.tiles_cutoff.mpi.1024 > tiles_cutoff.mpi.1024
grep "Rank" timings.tiles_cutoff.mpi.2048 > tiles_cutoff.mpi.2048
grep "Rank" timings.tiles_cutoff.mpi.4096 > tiles_cutoff.mpi.4096
grep "Rank" timings.tiles_cutoff.mpi.8192 > tiles_cutoff.mpi.8192

grep "Rank" timings.tiles.mpi.1024 > tiles.mpi.1024
grep "Rank" timings.tiles.mpi.2048 > tiles.mpi.2048
grep "Rank" timings.tiles.mpi.4096 > tiles.mpi.4096
grep "Rank" timings.tiles.mpi.8192 > tiles.mpi.8192

grep "Rank" timings.strips_cutoff.mpi.8192 > strips_cutoff.mpi.8192
grep "Rank" timings.strips_cutoff.mpi.4096 > strips_cutoff.mpi.4096
grep "Rank" timings.strips_cutoff.mpi.2048 > strips_cutoff.mpi.2048
grep "Rank" timings.strips_cutoff.mpi.1024 > strips_cutoff.mpi.1024

grep "Rank" timings.strips.mpi.8192 > strips.mpi.8192
grep "Rank" timings.strips.mpi.4096 > strips.mpi.4096
grep "Rank" timings.strips.mpi.2048 > strips.mpi.2048
grep "Rank" timings.strips.mpi.1024 > strips.mpi.1024

grep "Rank" timings.bh.mpi.8192 > bh.mpi.8192
grep "Rank" timings.bh.mpi.4096 > bh.mpi.4096
grep "Rank" timings.bh.mpi.2048 > bh.mpi.2048
grep "Rank" timings.bh.mpi.1024 > bh.mpi.1024

grep "Rank" timings.bh1.mpi.8192 > bh1.mpi.8192
grep "Rank" timings.bh1.mpi.4096 > bh1.mpi.4096
grep "Rank" timings.bh1.mpi.2048 > bh1.mpi.2048
grep "Rank" timings.bh1.mpi.1024 > bh1.mpi.1024

grep "Rank" timings.bh3.mpi.8192 > bh3.mpi.8192
grep "Rank" timings.bh3.mpi.4096 > bh3.mpi.4096
grep "Rank" timings.bh3.mpi.2048 > bh3.mpi.2048
grep "Rank" timings.bh3.mpi.1024 > bh3.mpi.1024
