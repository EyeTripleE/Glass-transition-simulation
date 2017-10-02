#!/bin/sh -x

MPI=mpiexec

function run_mpi() {
  ranks=(1 2 4 8 16 32 64)
  sizes=(1024 2048 4096 8192)
  case=$1

  for size in ${sizes[@]}; do
    for rnk in ${ranks[@]}; do
       $MPI -np $rnk ./n_body_simulation $size >> timings.${case}.mpi.${size}
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
       if [ $rnk -gt 8 ]; then
         export KMP_HOT_TEAMS_MAX_LEVEL=3
       else
         export KMP_HOT_TEAMS_MAX_LEVEL=2
       fi
       export KMP_PLACE_THREADS=1T
       export KMP_AFFINITY=compact,granularity=fine
       $MPI -np $rnk ./n_body_simulation $size >> timings.${case}.hybrid.${size}
    done
  done
}

#targets=(bh bh1 bh3 strips strips_cutoff tiles tiles_cutoff)
targets=(strips_cutoff tiles_cutoff)

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
