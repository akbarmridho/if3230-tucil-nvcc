#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void printMat(double **mat, int n) {
    for(int i = 0; i < n; ++i)
    {
        for(int j = n; j < 2*n; ++j)
        {
            printf("%lf ", mat[i][j]);
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
    double **mat = NULL;
    double d = 0.0;

    if (rank == 0) {
        scanf("%d", &n);
    }

    // Broadcast n
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int partition_size = n / size;
    int remainder = n % size;

    // Allocating memory for matrix array in all processes
    mat = (double**)malloc(n * sizeof(double*));
    for(i = 0; i < n; i++)
    {
        mat[i] = (double*)malloc(2 * n * sizeof(double));
    }

    if (rank == 0) {
        // Inputs the coefficients of the matrix
        for(i = 0; i < n; ++i)
        {
            for(j = 0; j < n; ++j)
            {
                scanf("%lf", &mat[i][j]);
            }
        }

        start_time = MPI_Wtime();
        
        // Initializing Right-hand side to identity matrix
        for(i = 0; i < n; ++i)
        {
            for(j = n; j < 2*n; ++j)
            {
                if(j == (i+n))
                {
                    mat[i][j] = 1;
                }
            }
        }
    }

    // Broadcast mat
    for(i = 0; i < n; i++)
    {
        MPI_Bcast(mat[i], 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    int start, end;

    // Reducing To Diagonal Matrix
    for(i = 0; i < n; ++i)
    {
        start = rank * partition_size + ((rank < remainder) ? rank : remainder);
        end = start + partition_size + ((rank < remainder) ? 1 : 0);

        for(j = start; j < end; ++j)
        {
            if(j != i)
            {
                d = mat[j][i] / mat[i][i];

                for(k = 0; k < n*2; k++)
                {
                    mat[j][k] -= mat[i][k]*d;
                }
            }
        }

        // Reduce row i
        if (i >= start && i < end) {
            d = mat[i][i];
            for(j = 0; j < 2*n; ++j) {
                mat[i][j] = mat[i][j]/d;
            }
        }

        // Broadcast mat[start until end]
        for(j = 0; j < n; j++)
        {
            int source_rank = j / partition_size;
            if (source_rank < remainder) source_rank++;

            MPI_Bcast(mat[j], 2 * n, MPI_DOUBLE, source_rank, MPI_COMM_WORLD);
        }
        
        // Synchronize all processes before the next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        // Get the end time
        double end_time = MPI_Wtime();

        // Compute the elapsed time
        double elapsed_time = end_time - start_time;

        printf("Elapsed time: %f seconds\n", elapsed_time);

        printf("\n=============RESULT=============\n");
        printMat(mat, n);
    }
    
    // Destructing
    for(i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);

    MPI_Finalize();

    return 0;
}
