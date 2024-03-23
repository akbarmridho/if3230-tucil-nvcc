#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void printMat(double *mat, int n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = n; j < 2 * n; ++j)
        {
            printf("%lf ", mat[i * 2 * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void eliminate_col_from_pivot(int row_start, int row_end, double *pivot_row, double *chunk, int col_start, int col_size)
{
    for (int i = row_start; i < row_end; i++)
    {
        double scale = chunk[i * col_size + col_start];
        for (int col = col_start; col < col_size; col++)
        {
            chunk[i * col_size + col] -= (pivot_row[col] * scale);
        }
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    double start_time;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i = 0, j = 0, k = 0, dim = 0;
    double d = 0.0;

    if (rank == 0)
    {
        scanf("%d", &dim);
    }

    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Column size
    int col_size = 2 * dim;
    
    // Rows to be processed by each processes
    int rows_cnt = dim / size; 
    if (dim % size != 0)
    {
        rows_cnt += 1;
    }

    // Initialize matrix
    double *mat = (double *)malloc(rows_cnt * size * col_size * sizeof(double));
    int start_row = rank * rows_cnt;
    int end_row = start_row + rows_cnt;

    // Store the chunk processed by each processes
    double *chunk = (double *)malloc(col_size * rows_cnt * sizeof(double));

    // Pivot row for each processes
    double *pivot_row = (double *)malloc(col_size * sizeof(double));

    // Input Matrix
    if (rank == 0)
    {
        // Inputs the coefficients of the matrix
        for (i = 0; i < dim; ++i)
        {
            for (j = 0; j < dim; ++j)
            {
                scanf("%lf", &mat[i * col_size + j]);
            }
        }

        // Initializing Right-hand side to identity matrix
        for (i = 0; i < dim; ++i)
        {
            for (j = dim; j < col_size; ++j)
            {
                if (j == (i + dim))
                {
                    mat[i * col_size + j] = 1;
                }
                else
                {
                    mat[i * col_size + j] = 0;
                }
            }
        }

        start_time = MPI_Wtime();

        // Begin swapping if needed
        for (int i = 0; i < dim; i++)
        {
            if (mat[i * col_size + i] == 0)
            {
                // Swap nearest row where mat[j][i] != 0
                for (int j = i + 1; j < dim; j++)
                {
                    if (mat[j * col_size + i] != 0.0)
                    {
                        for (int l = 0; l < col_size; l++)
                        {
                            double *row_a = &mat[i * col_size];
                            double *row_b = &mat[j * col_size];
                            double temp = row_a[l];
                            row_a[l] = row_b[l];
                            row_b[l] = temp;
                        }

                        break;
                    }
                    if (j == size - 1)
                    {
                        printf("Inverse does not exist for this matrix");
                        exit(0);
                    }
                }
            }
        }
    }

    // Partition
    MPI_Scatter(mat, col_size * rows_cnt, MPI_DOUBLE, chunk,
                col_size * rows_cnt, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Request *requests = (MPI_Request *)malloc(size * sizeof(MPI_Request));

    for (int row = 0; row < dim; row++)
    {
        int row_rank = row / rows_cnt;

        // If the row is supposed to be processed by current process
        if (rank == row_rank)
        {
            // curr_row is the row relative to the process rows
            int curr_row = row % rows_cnt; 
            double pivot = chunk[curr_row * col_size + row];

            for (int col = row; col < col_size; col++)
            {
                chunk[curr_row * col_size + col] /= pivot;
            }

            // Send the pivot to every other rows
            for (int i = 0; i < size; i++)
            {
                MPI_Isend(chunk + dim * curr_row * 2, col_size, MPI_DOUBLE, i, 0,
                          MPI_COMM_WORLD, &requests[i]);
            }

            // Perform Gaussian elimination on rows before the pivot row
            eliminate_col_from_pivot(0, curr_row, chunk + (curr_row * col_size), chunk, row, col_size);

            // Perform Gaussian elimination on rows after the pivot row
            eliminate_col_from_pivot(curr_row + 1, rows_cnt, chunk + (curr_row * col_size), chunk, row, col_size);

            for (int i = row_rank + 1; i < size; i++)
            {
                MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
            }
        }
        else
        {
            MPI_Recv(pivot_row, col_size, MPI_DOUBLE, row_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            eliminate_col_from_pivot(0, rows_cnt, pivot_row, chunk, row, col_size);
        }
    }

    MPI_Gather(chunk, rows_cnt * col_size, MPI_DOUBLE, mat, rows_cnt * col_size,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // Get the end time
        double end_time = MPI_Wtime();

        // Compute the elapsed time
        double elapsed_time = end_time - start_time;

        printf("Elapsed time: %f seconds\n", elapsed_time);

        printf("\n=============RESULT FROM PARALEL=============\n");
        printMat(mat, dim);
    }

    // Destructing
    free(mat);

    MPI_Finalize();

    return 0;
}
