#!/bin/sh

plot  'bh.hybrid.1024' using 3:7 with linespoints,   \
      'bh.hybrid.2048' using 3:7 with linespoints,   \
      'bh.hybrid.4096' using 3:7 with linespoints,   \
      'bh.hybrid.8192' using 3:7 with linespoints,   \
     'bh1.hybrid.1024' using 3:7 with linespoints,   \
     'bh1.hybrid.2048' using 3:7 with linespoints,   \
     'bh1.hybrid.4096' using 3:7 with linespoints,   \
     'bh1.hybrid.8192' using 3:7 with linespoints,   \
     'bh3.hybrid.1024' using 3:7 with linespoints,   \
     'bh3.hybrid.2048' using 3:7 with linespoints,   \
     'bh3.hybrid.4096' using 3:7 with linespoints,   \
     'bh3.hybrid.8192' using 3:7 with linespoints,   \
     'strips_cutoff.hybrid.1024' using 3:7 with linespoints,   \
     'strips_cutoff.hybrid.2048' using 3:7 with linespoints,   \
     'strips_cutoff.hybrid.4096' using 3:7 with linespoints,   \
     'strips_cutoff.hybrid.8192' using 3:7 with linespoints,   \
            'strips.hybrid.1024' using 3:7 with linespoints,   \
            'strips.hybrid.2048' using 3:7 with linespoints,   \
            'strips.hybrid.4096' using 3:7 with linespoints,   \
            'strips.hybrid.8192' using 3:7 with linespoints,   \
      'tiles_cutoff.hybrid.1024' using 3:7 with linespoints,   \
      'tiles_cutoff.hybrid.2048' using 3:7 with linespoints,   \
      'tiles_cutoff.hybrid.4096' using 3:7 with linespoints,   \
      'tiles_cutoff.hybrid.8192' using 3:7 with linespoints,   \
             'tiles.hybrid.1024' using 3:7 with linespoints,   \
             'tiles.hybrid.2048' using 3:7 with linespoints,   \
             'tiles.hybrid.4096' using 3:7 with linespoints,   \
             'tiles.hybrid.8192' using 3:7 with linespoints,   \
      'bh.mpi.1024' using 3:7 with linespoints,   \
      'bh.mpi.2048' using 3:7 with linespoints,   \
      'bh.mpi.4096' using 3:7 with linespoints,   \
      'bh.mpi.8192' using 3:7 with linespoints,   \
     'bh1.mpi.1024' using 3:7 with linespoints,   \
     'bh1.mpi.2048' using 3:7 with linespoints,   \
     'bh1.mpi.4096' using 3:7 with linespoints,   \
     'bh1.mpi.8192' using 3:7 with linespoints,   \
     'bh3.mpi.1024' using 3:7 with linespoints,   \
     'bh3.mpi.2048' using 3:7 with linespoints,   \
     'bh3.mpi.4096' using 3:7 with linespoints,   \
     'bh3.mpi.8192' using 3:7 with linespoints,   \
     'strips_cutoff.mpi.1024' using 3:7 with linespoints,   \
     'strips_cutoff.mpi.2048' using 3:7 with linespoints,   \
     'strips_cutoff.mpi.4096' using 3:7 with linespoints,   \
     'strips_cutoff.mpi.8192' using 3:7 with linespoints,   \
            'strips.mpi.1024' using 3:7 with linespoints,   \
            'strips.mpi.2048' using 3:7 with linespoints,   \
            'strips.mpi.4096' using 3:7 with linespoints,   \
            'strips.mpi.8192' using 3:7 with linespoints,   \
      'tiles_cutoff.mpi.1024' using 3:7 with linespoints,   \
      'tiles_cutoff.mpi.2048' using 3:7 with linespoints,   \
      'tiles_cutoff.mpi.4096' using 3:7 with linespoints,   \
      'tiles_cutoff.mpi.8192' using 3:7 with linespoints,   \
             'tiles.mpi.1024' using 3:7 with linespoints,   \
             'tiles.mpi.2048' using 3:7 with linespoints,   \
             'tiles.mpi.4096' using 3:7 with linespoints,   \
             'tiles.mpi.8192' using 3:7 with linespoints 


set logscale xy
set xlabel 'MPI Ranks'
set ylabel 'Time (s)'
set title 'Time taken for different algorithms'

replot

pause -1
