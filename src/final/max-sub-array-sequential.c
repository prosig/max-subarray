#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>

typedef struct tr
{
    int max_sum;
    int left;
    int right;
    int top;
    int bottom;
} thread_ret;

long
get_usecs(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return t.tv_sec*1000000+t.tv_usec;
}

void
usage(const char* app_name)
{
    printf("Argument error! Usage: %s <input_file 1> ... <input_file n>  \n", app_name);
    exit(0);
}

void
clear(
    int* a,
    int len)
{
    memset((void*)a, 0,len*sizeof(int));
}

/**
 * allocate matrix of size NxM on the heap
 */
int**
alloc_matrix(
    int n,
    int m)
{
    int **matrix;

    matrix = malloc(sizeof(*matrix) * n);

    if(!matrix)
    {
        printf("[ERROR] could not allocate memory for matrix\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            matrix[i] = malloc(sizeof *matrix[i] * m);
        }
    }

    return matrix;
}

/**
 * free matrix
 */
void
free_matrix(
    int** matrix,
    int n)
{
    if(!matrix)
    {
        printf("[ERROR] could not free memory for matrix - already freed\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            free(matrix[i]);
        }
    }

    free(matrix);
}

/**
 * return the matrix (NxN matrix) dimension as stated in the input file
 */
int
get_mat_dim(FILE* input_file)
{
    int dim = 0;

    fscanf(input_file, "%u\n", &dim);

    return dim;
}

/**
 * read the matrix from input file
 */
void
read_matrix(
    int **mat,
    int m,
    int n,
    FILE *input_file)
{
    int element = 0;

    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            if (j != (n-1))
            {
                fscanf(input_file, "%d\t", &element);
            }
            else
            {
                fscanf(input_file, "%d\n",&element);
            }
            mat[i][j] = element;
        }
    }
}

/**
 * print matrix
 */
void
print_matrix(
    int **matrix,
    int m,
    int n)
{
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

/**
 * precompiles matrix' vertical prefix sum
 */
void
precomp_matrix(
    int **mat,
    int **ps,
    int dim)
{
    /* precompute vertical prefix sum */
    for(int j=0; j<dim; j++)
    {
        ps[0][j] = mat[0][j];
        for (int i=1; i<dim; i++)
        {
            ps[i][j] = ps[i-1][j] + mat[i][j];
        }
    }
}

thread_ret
max_sub_arr(
    int **mat,
    int **ps,
    int **outmat,
    int dim)
{
    int max_sum = mat[0][0];
    int top = 0, left = 0, bottom = 0, right = 0;

    /* auxilliary variables */
    int sum[dim];
    int pos[dim];
    int local_max;

    /* precompute vertical prefix sum */
    precomp_matrix(mat, ps, dim);

    for (int i=0; i<dim; i++)
    {
        for (int k=i; k<dim; k++)
        {
            /* Kandane over all columns with the i..k rows */
            clear(sum, dim);
            clear(pos, dim);
            local_max = 0;

            /* we keep track of the position of the max value over each
             * Kandane's execution
             * -> Notice that we do not keep track of the max value,
             *    but only its position */
            sum[0] = ps[k][0] - (i==0 ? 0 : ps[i-1][0]);
            for (int j=1; j<dim; j++)
            {
                if (sum[j-1] > 0)
                {
                    sum[j] = sum[j-1] + ps[k][j] - (i==0 ? 0 : ps[i-1][j]);
                    pos[j] = pos[j-1];
                }
                else
                {
                    sum[j] = ps[k][j] - (i==0 ? 0 : ps[i-1][j]);
                    pos[j] = j;
                }
                if (sum[j] > sum[local_max])
                {
                    local_max = j;
                }
            } /* Kandane ends here */

            if (sum[local_max] > max_sum)
            {
                /* sum[local_max] is the new max value
                 * the corresponding submatrix goes from rows i..k.
                 * and from columns pos[local_max]..local_max */

                max_sum = sum[local_max];
                top = i;
                left = pos[local_max];
                bottom = k;
                right = local_max;
            }
        }
    }

    thread_ret ret;
    ret.max_sum = max_sum;
    ret.top = top;
    ret.left = left;
    ret.bottom = bottom;
    ret.right = right;

    return ret;
}

int
main(int argc, char* argv[])
{
    int dim = 0;
    int **mat=NULL, **ps=NULL, **outmat=NULL;
    FILE* input_file;

    if(argc < 2)
    {
        usage(argv[0]);
    }

    for(int i=1; i<argc; i++)
    {

        /* open files */
        input_file = fopen(argv[i], "r");
        if(input_file == NULL)
        {
            usage(argv[0]);
        }

        /* read matrix dimension */
        dim = get_mat_dim(input_file);

        /* init matrices */
        mat = alloc_matrix(dim, dim);
        ps = alloc_matrix(dim, dim);

        /* read matrix */
        read_matrix(mat, dim, dim, input_file);

        thread_ret ret = max_sub_arr(mat, ps, outmat, dim);
        printf("%d %d %d %d\n",ret.left,ret.top,ret.right,ret.bottom);

        /* release resources */
        fclose(input_file);
        free_matrix(mat, dim);
        free_matrix(ps, dim);
    }

    return EXIT_SUCCESS;
}
