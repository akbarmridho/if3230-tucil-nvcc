# Parallelized Program for Matrix Inversion using OpenMP

This program is a parallel implementation designed to find the inverse of a matrix. It utilizes the Gaussian elimination method and the OpenMP to efficiently perform computations in parallel.

## How it Works

1. **Initialization**: The program starts by reading the number of process that are going to used. The default value is `n_num = 4`.

2. **Reading the matrix dimension**: It will reads the dimension of the matrix (`dim`) from the standard input.

3. **Inputting the matrix**: The process will inputs the coefficients of the matrix and initializes the right-hand side to the identity matrix. This is done because the program is designed to find the inverse of the matrix. If any row in the matrix has a zero diagonal element, that row is swapped with a row where the corresponding element is non-zero. This is done to prevent division by zero during the Gaussian elimination process.

4. **Performing Gaussian elimination several threads**: Each threads performs Gaussian elimination on its own part of the matrix.

## Data Distribution Scheme

The program uses a row-wise block-striped partitioning scheme to distribute the rows of the matrix among the threads. This scheme was chosen because it allows each threads to independently perform computations on its part of the matrix, which improves parallel efficiency. Moreover, this scheme enhances data locality as each process works on a contiguous block of rows from the matrix. This improves cache performance, leading to faster execution times.

## Performance Test Results

The results, expressed in seconds, were obtained from tests conducted on a server.

| Matrix Size | Serial | OpenMP 2-core  | OpenMP 4-core  |
| ----------- | ------ | -------------- | -------------- |
| 32 x 32     | 0.00   | 0.00           | 0.00           |
| 64 x 64     | 0.02   | 0.00           | 0.00           |
| 128 x 128   | 0.06   | 0.02           | 0.03           |
| 256 x 256   | 0.42   | 0.14           | 0.14           |
| 512 x 512   | 2.82   | 0.92           | 0.94           |
| 1024 x 1024 | 21.38  | 6.67           | 6.71           |
| 2048 x 2048 | 160.74 | 51.35          | 51.30          |

