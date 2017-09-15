#!/bin/sh -x


function run_mpi() {
  ranks=(1 2 4 8 16 32 64)
  sizes=(1024 2048 4096 8192)
  case=$1

  for size in ${sizes[@]}; do
    for rnk in ${ranks[@]}; do
       mpirun -np $rnk n_body_simulation >> timings.${case}.mpi.${size}
    done
  done
}

function run_hybrid() {
  ranks=(1 2 4 8 16 32 64)
  sizes=(1024 2048 4096 8192)
  case=$1

  for size in ${sizes[@]}; do
    for rnk in ${ranks[@]}; do
       export KMP_HOT_TEAMS_MODE=1
       export KMP_HOT_TEAMS_MAX_LEVEL=2
       export KMP_PLACE_THREADS = 1T
       export KMP_AFFINITY=compact,granularity=fine
       mpirun -np $rnk n_body_simulation >> timings.${case}.hybrid.${size}
    done
  done
}

targets=(bh bh1 bh3 strips strips_cutoff tiles tiles_cutoff)

for t in ${targets[@]}; do
  make clean
  echo "Building target $t"
  make $t &
  pid=$!
  wait $pid

  echo 'Running MPI benchmarks'
  run_mpi $t

  make clean
  echo "Building target $t with MPI+OpenMP"
  make USE_OMP=1 ${t} &
  pid=$!
  wait $pid

  echo 'Running OpenMP benchmarks'
  run_hybrid $t
done
