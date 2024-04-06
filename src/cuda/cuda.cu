#include <stdio.h>
#include <sys/time.h>

#define blocksize 64

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

__global__ void swap_zero_diagonal_and_get_scale(double *mat, int i, int dim, double *scale, int *status)
{
    int col_size = 2 * dim;

    if (mat[i * col_size + 1] == 0)
    {
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
            if (j == dim - 1)
            {
                *status = 1;
                return;
            }
        }
    }

    *status = 0;
    *scale = mat[i * 2 * dim + i];
}

__global__ void normalize_diagonal(double *mat, int i, int dim, double *scale)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;

    if (x < dim)
    {
        mat[i * 2 * dim + x] /= *scale;
        mat[i * 2 * dim + dim + x] /= *scale;
    }
}

__global__ void perform_elimination(int row_start, int row_end, int pivot_idx, double *mat, int col_start, int dim)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int i = x + row_start;

    if (i < row_end)
    {
        double scale = mat[i * 2 * dim + col_start];

        for (int col = col_start; col < 2 * dim; col++)
        {
            mat[i * 2 * dim + col] -= mat[2 * dim * pivot_idx + col] * scale;
        }
    }
}

void eliminate_col_from_pivot(int row_start, int row_end, int pivot_idx, double *mat, int col_start, int dim)
{
    int row_size = row_end - row_start;

    int gridsize = row_size / blocksize;

    if (row_size % blocksize != 0)
        gridsize++;

    perform_elimination<<<blocksize, gridsize>>>(row_start, row_end, pivot_idx, mat, col_start, dim);
}

int main(int argc, char *argv[])
{
    int dim = 0;

    scanf("%d", &dim);

    // case when dim = 3
    // core = 2
    int col_size = 2 * dim;
    int row_size = dim;

    // initialize matrix
    double *mat = (double *)malloc(row_size * col_size * sizeof(double));

    // scan matrix
    for (int i = 0; i < row_size; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            scanf("%lf", &mat[i * col_size + j]);
        }
    }

    // Initializing Right-hand side to identity matrix
    for (int i = 0; i < dim; ++i)
    {
        for (int j = dim; j < col_size; ++j)
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

    // initialize cuda

    double *mat_d;
    cudaMalloc(&mat_d, row_size * col_size * sizeof(double));
    cudaMemcpy(mat_d, mat, row_size * col_size * sizeof(double), cudaMemcpyHostToDevice);

    double *scale_d;
    cudaMalloc(&scale_d, sizeof(double));

    int *status_d;
    int status;
    cudaMalloc(&status_d, sizeof(int));

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    for (int i = 0; i < dim; i++)
    {
        swap_zero_diagonal_and_get_scale<<<1, 1>>>(mat_d, i, dim, scale_d, status_d);

        cudaMemcpy(&status, status_d, sizeof(int), cudaMemcpyDeviceToHost);

        if (status != 0)
        {
            printf("No inverse exist\n");
            exit(0);
        }

        int gridsize = dim / blocksize;

        if (dim % blocksize != 0)
            gridsize++;

        normalize_diagonal<<<blocksize, gridsize>>>(mat_d, i, dim, scale_d);

        if (i == dim - 1)
            continue;

        eliminate_col_from_pivot(i + 1, dim, i, mat_d, i, dim);
    }

    for (int i = dim - 1; i >= 1; i--)
    {
        eliminate_col_from_pivot(0, i, i, mat_d, i, dim);
    }

    cudaDeviceSynchronize();

    cudaMemcpy(mat, mat_d, row_size * col_size * sizeof(double), cudaMemcpyDeviceToHost);

    gettimeofday(&t2, NULL);
    double elapsed_time = t2.tv_sec - t1.tv_sec;

    printf("Elapsed time: %f seconds\n", elapsed_time);

    printf("\n=============RESULT FROM CUDA=============\n");
    printMat(mat, dim);

    free(mat);
    cudaFree(mat_d);
    cudaFree(scale_d);
    cudaFree(status_d);

    return 0;
}