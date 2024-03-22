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

    // case when dim = 3
    // core = 2
    int col_size = 2 * dim;
    // initialize matrix

    int rows_cnt = dim / size; // the numbber of rows to be processed by each process
    if (dim % size != 0)
    {
        rows_cnt += 1;
    }

    double *mat = (double *)malloc(rows_cnt * size * col_size * sizeof(double));
    int start_row = rank * rows_cnt;
    int end_row = start_row + rows_cnt;

    // // the variable to store the chunk processed by each proceess
    double *chunk = (double *)malloc(col_size * rows_cnt * sizeof(double));

    // // the pivot row for each process
    double *pivot_row = (double *)malloc(col_size * sizeof(double));

    // // input the m in main process

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
    }

    // partition
    MPI_Scatter(mat, col_size * rows_cnt, MPI_DOUBLE, chunk,
                col_size * rows_cnt, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Request *requests = (MPI_Request *)malloc(size * sizeof(MPI_Request));

    for (int row = 0; row < dim; row++)
    {
        int row_rank = row / rows_cnt;
        // if the row is supposed to be processed by current rank
        if (rank == row_rank)
        {
            int curr_row = row % rows_cnt; // curr row is the row relative to the process rows
            double pivot = chunk[curr_row * col_size + row];
            for (int col = row; col < col_size; col++)
            {
                chunk[curr_row * col_size + col] /= pivot;
            }

            // send the pivot to every other rows

            for (int i = 0; i < size; i++)
            {
                MPI_Isend(chunk + dim * curr_row * 2, col_size, MPI_DOUBLE, i, 0,
                          MPI_COMM_WORLD, &requests[i]);
            }
            // eliminate the rows before it
            eliminate_col_from_pivot(0, curr_row, chunk + (curr_row * col_size), chunk, row, col_size);

            // eliminate the rows after it
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

    // end_row = start_row;

    // for (int row = dim - 1; row >= end_row; row--)
    // {
    //     // row_rank is the process that is supposed to be dealing with row
    //     int row_rank = row / rows_cnt;
    //     // if the row is supposed to be processed by current rank
    //     if (rank == row_rank)
    //     {
    //         int curr_row = row % rows_cnt; // curr row is the row relative to the process rows
    //         // send the pivot to every other rows

    //         for (int i = row_rank - 1; i >= 0; i--)
    //         {
    //             MPI_Isend(chunk + dim * curr_row * 2, col_size, MPI_DOUBLE, i, 0,
    //                       MPI_COMM_WORLD, &requests[i]);
    //         }

    //         // set the cols of the pivot to 0
    //         eliminate_col_from_pivot(0, curr_row, chunk + (curr_row * col_size), chunk, row, col_size);

    //         for (int i = row_rank - 1; i >= 0; i--)
    //         {
    //             MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
    //         }
    //     }
    //     else
    //     {
    // MPI_Recv(pivot_row, col_size, MPI_DOUBLE, row_rank, 0, MPI_COMM_WORLD,
    //          MPI_STATUS_IGNORE);
    //         eliminate_col_from_pivot(0, rows_cnt, pivot_row, chunk, row, col_size);
    //     }
    // }

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

    // // Destructing
    free(mat);

    MPI_Finalize();

    return 0;
}
