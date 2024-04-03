#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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

#pragma omp for
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
    int n_num = 4;

    if (argc == 2)
    {
        n_num = atoi(argv[1]);
    }

    int i = 0, j = 0, k = 0, dim = 0;

    scanf("%d", &dim);

    // case when dim = 3
    // core = 2
    int col_size = 2 * dim;
    int row_size = dim;

    // initialize matrix
    double *mat = (double *)malloc(row_size * col_size * sizeof(double));

    // scan matrix
    for (i = 0; i < row_size; ++i)
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

    omp_set_num_threads(n_num);

    double start_time = omp_get_wtime();

#pragma omp parallel
    {
        for (int i = 0; i < dim; i++)
        {
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
                        printf("Inverse does not exist for this matrix");
                        exit(0);
                    }
                }
            }

            double scale = mat[i * col_size + i];
            for (int j = 0; j < col_size; j++)
            {
                mat[i * col_size + j] /= scale;
            }

            if (i == dim - 1)
                continue;
            eliminate_col_from_pivot(i + 1, dim, mat + i * col_size, mat, i, col_size);
        }
        for (int i = dim - 1; i >= 1; i--)
        {

            eliminate_col_from_pivot(0, i, mat + i * col_size, mat, i, col_size);
        }
    }
    double end_time = omp_get_wtime();

    double elapsed_time = end_time - start_time;

    printf("Elapsed time: %f seconds\n", elapsed_time);

    printf("\n=============RESULT FROM PARALEL=============\n");
    printMat(mat, dim);

    // // Destructing
    free(mat);

    return 0;
}
