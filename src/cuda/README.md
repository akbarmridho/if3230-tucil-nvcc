# Parallelized Program for Matrix Inversion using OpenMP

This program is a parallel implementation designed to find the inverse of a matrix. It utilizes the Gaussian elimination method and the OpenMP to efficiently perform computations in parallel.

## How it Works

1. **Reading the matrix dimension**: It will reads the dimension of the matrix (`dim`) from the standard input.

2. **Inputting the matrix**: The process will inputs the coefficients of the matrix and initializes the right-hand side to the identity matrix. This is done because the program is designed to find the inverse of the matrix. If any row in the matrix has a zero diagonal element, that row is swapped with a row where the corresponding element is non-zero. This is done to prevent division by zero during the Gaussian elimination process.

3. **Initializing Cuda**: The process will allocate memory for matrix and two other variable. After allocating memory, the process will copy the matrix from host to device.

4. **Performing Gaussian elimination**: This algorithm wil be divided into two parts. The first part is looping for dimension times to perform diagonal normalization and gaussian elimination. This will result in upper triangle matrix. The second part is back substitution that loop for dimension times. Each diagonal normalization and gaussian elimination is done in parallel on the GPU.

## Data Distribution Scheme

The program use global GPU memory to store matrix and result. The data distribution is configured with thread index, block dimension, and block index so that each thread can perform computation on its own part.

## Performance Test Results

Test environment: AMD Ryzen 7 4800H and Nvidia Geforce RTX 3050 on Ubuntu 22.04 (WSL2) and Cuda 12.3.

The results, expressed in seconds, were obtained from tests conducted on a server.

| Matrix Size | Serial | Cuda GPU       |
| ----------- | ------ | -------------- |
| 32 x 32     | 0.00   | 0.47           |
| 64 x 64     | 0.00   | 0.38           |
| 128 x 128   | 0.04   | 0.40           |
| 256 x 256   | 0.26   | 0.60           |
| 512 x 512   | 1.92   | 0.94           |
| 1024 x 1024 | 14.19  | 3.73           |
| 2048 x 2048 | 109.46 | 19.46          |

