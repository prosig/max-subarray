#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <pthread.h>

#define NUM_THREADS 48


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
    printf("Argument error! Usage: %s <input_file 1> ... <input_file n>\n", app_name);
    exit(0);
}

void
clear(
    int* a,
    int len)
{
    memset((void*)a, 0, len);
}

int**
alloc_matrix(
    int n,
    int m)
{
    int **matrix;

    matrix = malloc(sizeof(*matrix) * n);

    if(!matrix)
    {
        printf("[ERROR] could not allocate memory for %dx%d-matrix\n",n,m);
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

int
get_mat_dim(FILE* input_file)
{
    int dim = 0;

    fscanf(input_file, "%u\n", &dim);

    return dim;
}

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

typedef struct precomp_ta
{
    int threadId;
    int start;
    int end;
    int dim;
    int** mat;
    int** ps;
} precomp_thread_args;

void* precomp_pthread_func(void* args)
{
    precomp_thread_args targs = *(precomp_thread_args*)args;
    int** mat = targs.mat;
    int** ps = targs.ps;
    int dim = targs.dim;
    int start = targs.start;
    int end = targs.end;

    for(int j=start; j<end; j++)
    {
        ps[0][j] = mat[0][j];
        for (int i=1; i<dim; i++)
        {
            ps[i][j] = ps[i-1][j] + mat[i][j];
        }
    }
}

void precomp_matrix(int **mat, int **ps, int dim, int num_threads)
{
    int iterations = dim;
    int n=0;
    int step = iterations / num_threads;
    int over = iterations % num_threads;
    pthread_t threads[num_threads];
    precomp_thread_args thread_args[num_threads];

    for(int i=0; i<iterations; i+=step,n++)
    {
        thread_args[n].threadId = n;
        thread_args[n].dim = dim;
        thread_args[n].start = i;
        if(over-->0) i++;
        thread_args[n].end = i + step > iterations ? iterations : i + step;
        thread_args[n].ps = ps;
        thread_args[n].mat = mat;

        int error = pthread_create (&threads[n], 0, precomp_pthread_func, (void *)&thread_args[n]);
        if(error)
        {
            printf("[ERROR %d] Couldn't create thread with threadId=%d,dimension=%d,start=%d,end=%d\n",error,thread_args[n].threadId,thread_args[n].dim,thread_args[n].start,thread_args[n].end);
            exit(EXIT_FAILURE);
        }
    }

    for(int i=0; i<n; i++)
    {
        int error = pthread_join(threads[i],0);
        if(error)
        {
            printf("[ERROR %d] Couldn't join thread with threadId=%d\n",error,i);
            exit(EXIT_FAILURE);
        }
    }
}

typedef struct tr
{
    int max_sum;
    int left;
    int right;
    int top;
    int bottom;
} thread_ret;

typedef struct ta
{
    int threadId;
    int dimension;
    int start;
    int end;
    int** ps;
    int* iter_i;
    int mat0;
    thread_ret result;
} thread_args;


thread_ret maxArray_pthread(int** ps, int** mat, int dim, int num_threads);
void* maxArray_pthread_func(void* tdata);

thread_ret maxArray_pthread(int** ps, int** mat, int dim, int num_threads)
{
    int iterations = dim * (dim + 1) / 2;

    int iter_i[dim];
    int sum = 0;
    for(int i=0; i<dim; sum+=dim-i,i++) iter_i[i] = sum;

    pthread_t threads[num_threads];
    thread_args thread_args[num_threads];
    int n=0;
    int step = iterations / num_threads + 1;
    for(int i=0; i<iterations; i+=step,n++)
    {
        thread_args[n].threadId = n;
        thread_args[n].dimension = dim;
        thread_args[n].start = i;
        thread_args[n].end = i + step > iterations ? iterations : i + step;
        thread_args[n].iter_i = iter_i;
        thread_args[n].ps = ps;
        thread_args[n].mat0 = mat[0][0];

        int error = pthread_create (&threads[n], 0, maxArray_pthread_func, (void *)&thread_args[n]);
        if(error)
        {
            printf("[ERROR %d] Couldn't create thread with threadId=%d,dimension=%d,start=%d,end=%d\n",error,thread_args[n].threadId,thread_args[n].dimension,thread_args[n].start,thread_args[n].end);
            exit(EXIT_FAILURE);
        }
    }

    thread_ret ret;
    ret.max_sum = mat[0][0];
    ret.top = 0;
    ret.left = 0;
    ret.bottom = 0;
    ret.right = 0;
    for(int i=0; i<n; i++)
    {
        int error = pthread_join(threads[i],0);
        if(error)
        {
            printf("[ERROR %d] Couldn't join thread with threadId=%d\n",error,i);
            exit(EXIT_FAILURE);
        }

        if(thread_args[i].result.max_sum > ret.max_sum) ret = thread_args[i].result;
    }

    return ret;
}

void* maxArray_pthread_func(void* tdata)
{
    thread_args targs = *(thread_args*)tdata;
    int tid = targs.threadId;
    int dim = targs.dimension;
    int start = targs.start;
    int end = targs.end;
    int* iter_i = targs.iter_i;
    int** ps = targs.ps;

    int max_sum = targs.mat0;
    int top = 0, left = 0, bottom = 0, right = 0;

    //Auxilliary variables
    int sum[dim];
    int pos[dim];
    int local_max;

    int i = 0;
    for(; i<dim&&iter_i[i]<=start; i++);
    int next_i = i<dim ? iter_i[i] : end;
    i--;

    for(int n=start; n<end; n++)
    {
        if(n >= next_i)
        {
            i++;
            next_i = i<dim-1 ? iter_i[i+1] : end;
        }

        int k = i + n - iter_i[i];

        // Kandane over all columns with the i..k rows
        clear(sum, dim);
        clear(pos, dim);
        local_max = 0;

        // We keep track of the position of the max value over each Kandane's execution
        // Notice that we do not keep track of the max value, but only its position

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
        } //Kandane ends here

        if (sum[local_max] > max_sum)
        {
            // sum[local_max] is the new max value
            // the corresponding submatrix goes from rows i..k.
            // and from columns pos[local_max]..local_max

            max_sum = sum[local_max];
            top = i;
            left = pos[local_max];
            bottom = k;
            right = local_max;
        }
    }

    thread_ret ret;
    ret.max_sum = max_sum;
    ret.top = top;
    ret.left = left;
    ret.bottom = bottom;
    ret.right = right;

    ((thread_args*)tdata)->result = ret;

    return 0;
}

int
main(int argc, char* argv[])
{
    int dim = 0;
    int **mat, **ps, **outmat;
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

        precomp_matrix(mat, ps, dim, NUM_THREADS);

        thread_ret ret = maxArray_pthread(ps, mat, dim, NUM_THREADS);
        printf("%d %d %d %d\n",ret.left,ret.top,ret.right,ret.bottom);


        /* release resources */
        fclose(input_file);
        free_matrix(mat, dim);
        free_matrix(ps, dim);
    }

    return EXIT_SUCCESS;
}
