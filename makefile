OUTPUT_FOLDER = bin

all: serial mpi omp

mpi:
	mpicc src/open-mpi/mpi.c -o $(OUTPUT_FOLDER)/mpi

omp:
	gcc src/open-mp/open-mp.c --openmp -o $(OUTPUT_FOLDER)/omp

serial:
	g++ src/serial/serial.cpp -o $(OUTPUT_FOLDER)/serial