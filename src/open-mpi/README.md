# Parallelized Program for Matrix Inversion using OpenMPI

This program is a parallel implementation designed to find the inverse of a matrix. It utilizes the Gaussian elimination method and the Open Message Passing Interface (OpenMPI) to efficiently perform computations in parallel.

## How it Works

1. **Initialization**: The program starts by initializing the MPI environment and getting the total number of processes (`size`) and the rank of the current process (`rank`).

2. **Reading the matrix dimension and broadcasting it**: If the current process is the root process (rank 0), it reads the dimension of the matrix (`dim`) from the standard input. It then broadcasts this value to all other processes in the MPI communicator using the `MPI_Bcast` function.

3. **Calculating the number of rows each process will handle**: Each process calculates the number of rows it will handle (`rows_cnt`). This is done by dividing the total number of rows (`dim`) by the number of processes (`size`). If the total number of rows is not evenly divisible by the number of processes, an extra row is added to some processes to ensure all rows are handled.

4. **Initializing the matrix chunk and pivot row**: Each process allocates memory for its chunk of the matrix and a pivot row. The chunk is a portion of the matrix that the process will work on, and the pivot row is used during the Gaussian elimination process.

5. **Inputting the matrix**: If the current process is the root process, it inputs the coefficients of the matrix and initializes the right-hand side to the identity matrix. This is done because the program is designed to find the inverse of the matrix. If any row in the matrix has a zero diagonal element, that row is swapped with a row where the corresponding element is non-zero. This is done to prevent division by zero during the Gaussian elimination process.

6. **Distributing the matrix chunks to each process**: The `MPI_Scatter` function is used to distribute equal chunks of the matrix to each process. This is a communication operation where data is sent from one process (the root) to all processes in a communicator, including itself.

7. **Performing Gaussian elimination on each chunk**: Each process performs Gaussian elimination on its chunk of the matrix. If the pivot row is in the current process's chunk, it normalizes the pivot row (divides every element by the pivot element) and sends this pivot row to all other processes using `MPI_Isend`. If the pivot row is not in the current process's chunk, it receives the pivot row from the appropriate process using `MPI_Recv` and uses it to eliminate the corresponding column in its chunk.

8. **Collecting the chunks from all processes**: After all processes have finished their computations, the `MPI_Gather` function is used to collect the chunks from all processes back to the root process.

This process is repeated for all rows in the matrix. Because each process is working on a different set of rows at the same time, the Gaussian elimination is performed in parallel.

## Data Distribution Scheme

The program uses a row-wise block-striped partitioning scheme to distribute the rows of the matrix among the processes. This scheme was chosen because it allows each process to independently perform computations on its chunk of the matrix, which improves parallel efficiency. Moreover, this scheme enhances data locality as each process works on a contiguous block of rows from the matrix. This improves cache performance and reduces the time spent transferring data, leading to faster execution times. Additionally, this scheme minimizes communication overhead between processes, further improving the efficiency of the program.

## Performance Test Results

The results, expressed in seconds, were obtained from tests conducted on a server.

| Matrix Size | Serial | OpenMPI 2-core | OpenMPI 4-core |
| ----------- | ------ | -------------- | -------------- |
| 32 x 32     | 0.00   | 1.46           | 0.70           |
| 64 x 64     | 0.02   | 0.69           | 0.70           |
| 128 x 128   | 0.06   | 0.70           | 0.72           |
| 256 x 256   | 0.42   | 0.82           | 0.84           |
| 512 x 512   | 2.82   | 1.44           | 1.25           |
| 1024 x 1024 | 21.38  | 5.48           | 3.96           |
| 2048 x 2048 | 160.74 | 36.71          | 23.04          |
