# Tugas Kecil - Paralel Inverse Matrix

Tugas Kecil IF3230 Sistem Paralel dan Terdistribusi

- Rachel Gabriela Chen (13521044)
- Jeffrey Chow (13521046)
- Akbar Maulana Ridho (13521093)

## How to Compile

### Compile all

Run `make` on root folder.

### Compile serial

Run `make serial` or `g++ src/serial/serial.cpp -o bin/serial` on root folder.

### Compile OpenMPI

Make sure OpenMPI is installed. Run `sudo apt install openmpi-bin libopenmpi-dev` to install OpenMPI.

Run `make mpi` or `mpicc src/open-mpi/mpi.c -o bin/mpi` on root folder.

### Compile OpenMP

Run `make omp` or `gcc src/open-mp/open-mp.c -fopenmp -o bin/omp` on root folder.

### Compile CUDA

Run `make cuda` or `nvcc src/cuda/cuda.cu -o bin/cuda` on root folder.

## How to Run

### Serial

Manual

`cat ./test_cases/{test_case}.txt | ./bin/serial`

Script

`sh ./test_serial.sh`

### OpenMPI

Manual

`cat ./test_cases/{test_case}.txt | mpirun -n {n_core} ./bin/mpi`

Script

`sh ./test_mpi.sh`

### OpenMP

Manual

`cat ./test_cases/{test_case}.txt | ./bin/omp {n_core}`

Script

`sh ./test_omp.sh`