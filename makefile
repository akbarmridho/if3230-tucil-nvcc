OUTPUT_FOLDER = bin

all: serial mpi omp cuda

mpi:
	mpicc src/open-mpi/mpi.c -o $(OUTPUT_FOLDER)/mpi

omp:
	gcc src/open-mp/open-mp.c -fopenmp -o $(OUTPUT_FOLDER)/omp

serial:
	g++ src/serial/serial.cpp -o $(OUTPUT_FOLDER)/serial

cuda:
	nvcc src/cuda/cuda.cu -o $(OUTPUT_FOLDER)/cuda