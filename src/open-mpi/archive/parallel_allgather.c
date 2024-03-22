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

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    double start_time;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i = 0, j = 0, k = 0, n = 0;
    double *mat = NULL;
    double d = 0.0;

    if (rank == 0)
    {
        scanf("%d", &n);
    }

    // Broadcast n
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int partition_size = n / size;

    // Allocating memory for matrix array in all processes
    mat = (double *)malloc(n * 2 * n * sizeof(double));

    if (rank == 0)
    {
        // Inputs the coefficients of the matrix
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                scanf("%lf", &mat[i * 2 * n + j]);
            }
        }

        // Initializing Right-hand side to identity matrix
        for (i = 0; i < n; ++i)
        {
            for (j = n; j < 2 * n; ++j)
            {
                if (j == (i + n))
                {
                    mat[i * 2 * n + j] = 1;
                }
            }
        }

        start_time = MPI_Wtime();
    }

    // Broadcast mat
    MPI_Bcast(mat, 2 * n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int start, end;

    // Reducing To Diagonal Matrix
    for (i = 0; i < n; ++i)
    {
        start = rank * partition_size;
        end = start + partition_size;

        // printf("\nRank: %d, Start:%d, End:%d\n", rank, start, end);

        for (j = start; j < end; ++j)
        {
            if (j != i)
            {
                d = mat[j * 2 * n + i] / mat[i * 2 * n + i];

                for (k = 0; k < n * 2; k++)
                {
                    mat[j * 2 * n + k] -= mat[i * 2 * n + k] * d;
                }
            }
        }

        // row partisi

        // Reduce row i
        if (i >= start && i < end)
        {
            d = mat[i * 2 * n + i];
            for (j = 0; j < 2 * n; ++j)
            {
                mat[i * 2 * n + j] = mat[i * 2 * n + j] / d;
            }
        }

        // row i aja

        // // Broadcast mat
        // for (j = 0; j < size; ++j)
        // {
        //     start = j * partition_size * 2 * n;

        //     MPI_Bcast(&mat[start], 2 * n * partition_size, MPI_DOUBLE, j, MPI_COMM_WORLD);
        // }

        MPI_Allgather(
            &mat[2 * n * start],    // send buffer
            2 * n * partition_size, // send count
            MPI_DOUBLE,             // send type
            mat,                    // recv buffer
            2 * n * partition_size, // recv count
            MPI_DOUBLE,             // recv datatype
            MPI_COMM_WORLD);
    }

    if (rank == 0)
    {
        // Get the end time
        double end_time = MPI_Wtime();

        // Compute the elapsed time
        double elapsed_time = end_time - start_time;

        printf("Elapsed time: %f seconds\n", elapsed_time);

        printf("\n=============RESULT=============\n");

        printMat(mat, n);
    }

    // Destructing
    free(mat);

    MPI_Finalize();

    return 0;
}
