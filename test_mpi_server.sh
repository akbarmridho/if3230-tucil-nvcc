#!/bin/bash

echo "OpenMPI server testing"

for case in 32 64 128 256 512 1024 2048
do
    for n_core in 2 4 8
    do
        echo "Case $case x $case with n $n_core";
        time -f %e sh -c "cat ./test_cases/$case.txt | mpirun --hostfile hostfile -n $n_core ./bin/mpi > /dev/null";
    done;
done