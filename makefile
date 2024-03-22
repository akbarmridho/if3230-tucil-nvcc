OUTPUT_FOLDER = bin

all: serial mpi

mpi:
	mpicc src/open-mpi/mpi.c -o $(OUTPUT_FOLDER)/mpi

serial:
	g++ src/serial/serial.cpp -o $(OUTPUT_FOLDER)/serial