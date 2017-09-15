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
       export OMP_NUM_THREADS=$(( 64/rnk ))
       export OMP_THREAD_LIMIT=$OMP_NUM_THREADS
       export OMP_NESTED=TRUE
       echo "Thilina: $OMP_NUM_THREADS"
       mpirun -np $rnk n_body_simulation >> timings.${case}.hybrid.${size}
    done
  done
}

targets=(bh strips strips_cutoff tiles tiles_cutoff)

for t in ${targets[@]}; do
##  make clean
##  echo "Building target $t"
##  make $t &
##  pid=$!
##  wait $pid
##
##  echo 'Running MPI benchmarks'
##  run_mpi $t

  make clean
  echo "Building target $t with MPI+OpenMP"
  make USE_OMP=1 ${t} &
  pid=$!
  wait $pid

  echo 'Running OpenMP benchmarks'
  run_hybrid $t
done
