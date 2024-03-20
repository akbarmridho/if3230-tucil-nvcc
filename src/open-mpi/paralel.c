#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

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

void printArr(double* arr, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%lf ", arr[i]);
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

    // Allocating memory for matrix array in all processes
    mat = (double**)malloc(n * sizeof(double*));
    for(i = 0; i < n; ++i)
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

        // Partial pivoting, TODO: check again, suspicious
        for(i = n - 1; i > 0; --i)
        {
            if(mat[i-1][1] < mat[i][1])
            {
                for(j = 0; j < 2*n; ++j)
                {
                    d = mat[i][j];
                    mat[i][j] = mat[i-1][j];
                    mat[i-1][j] = d;
                }
            }
        }
    }

    // Broadcast mat
    for(i = 0; i < n; ++i)
    {
        MPI_Bcast(mat[i], 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Reducing To Diagonal Matrix
    MPI_Request *requests = (MPI_Request*)malloc(n * sizeof(MPI_Request));
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < n; ++j)
        {
            requests[j] = MPI_REQUEST_NULL;
        }

        for(j = 0; j < n; ++j)
        {
            if(j != i && j % size == rank)
            {
                d = mat[j][i] / mat[i][i];

                for(k = 0; k < n*2; k++)
                {
                    mat[j][k] -= mat[i][k]*d;
                }

                // Non-blocking send of the changed row mat[j] to root
                if (rank != 0) {
                    MPI_Isend(mat[j], 2 * n, MPI_DOUBLE, 0, j, MPI_COMM_WORLD, &requests[j]);
                }
            }
        }

        // Wait for all non-blocking send operations to complete
        MPI_Waitall(n, requests, MPI_STATUSES_IGNORE);

        if (rank == 0) {
            // If root, non-blocking receive of all the changed data
            for(j = 0; j < n; ++j)
            {
                if (j != i && j % size != rank) {
                    MPI_Irecv(mat[j], 2 * n, MPI_DOUBLE, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, &requests[j]);
                }
            }

            // Wait for all non-blocking receive operations to complete
            MPI_Waitall(n, requests, MPI_STATUSES_IGNORE);
        }
        
        // Then rebroadcast it for next iteration
        for(j = 0; j < n; ++j)
        {
            MPI_Bcast(mat[j], 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        // Synchronize all processes before the next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }
    free(requests);

    if (rank == 0) {
        for(i = 0; i < n; ++i) {
            d = mat[i][i];
            for(j = 0; j < 2*n; ++j) {
                mat[i][j] = mat[i][j]/d;
            }
        }

        // Get the end time
        double end_time = MPI_Wtime();

        // Compute the elapsed time
        double elapsed_time = end_time - start_time;

        printf("Elapsed time: %f seconds\n", elapsed_time);

        printf("\n=============RESULT=============\n");
        printMat(mat, n);
    }
    
    // Destructing
    for(i = 0; i < n; ++i)
    {
        free(mat[i]);
    }
    free(mat);

    MPI_Finalize();

    return 0;
}
