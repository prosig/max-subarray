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

#define get_data(type_t,t_arg) (*(type_t*)((*(thread_arg*)t_arg).data))
#define get_start(t_arg) ((*(thread_arg*)t_arg).start)
#define get_end(t_arg) ((*(thread_arg*)t_arg).end)
#define get_tid(t_arg) ((*(thread_arg*)t_arg).tid)
#define write_result(type_t,t_arg,data) (*(type_t*)((*(thread_arg*)t_arg).result)=data)
#define divide(iterations,num_threads,thread_func,args) divide_and_reduce(iterations,num_threads,thread_func,(void*)&args,reduce_nothing,0,0)
#define divide_reduce(iterations,num_threads,thread_func,args,reduce_func,result_out) divide_and_reduce(iterations,num_threads,thread_func,(void*)&args,reduce_func,sizeof(result_out),(void*)&result_out)

void* recursive_reduce(void*(*reduce_func)(void*,void*),void* data, int n1, int n2, int element_size)
{
    if(n2 <= n1) return (void*)((char*)data+n1*element_size);
    else return reduce_func(recursive_reduce(reduce_func,data,n1,n1+(n2-n1)/2,element_size),recursive_reduce(reduce_func,data,1+n1+(n2-n1)/2,n2,element_size));
}

void* reduce_nothing(void* n1,void* n2)
{
    return 0;
}

typedef struct
{
    int tid;
    int start;
    int end;
    void* data;
    void* result;
} thread_arg;

void divide_and_reduce(int iterations, int num_threads, void*(*thread_func)(void*), void* data, void*(*reduce_func)(void*,void*), int result_element_size, void* result_out)
{
    int step = iterations / num_threads;
    int over = iterations % num_threads;
    pthread_t threads[num_threads];
    thread_arg thread_args[num_threads];
    char results[result_element_size*num_threads];
    int n=0;
    for(int i=0; i<iterations; i+=step, n++)
    {
        thread_args[n].tid = n;
        thread_args[n].data = data;
        thread_args[n].result = (void*)(results+n*result_element_size);
        thread_args[n].start = i;
        if(over-->0) i++;
        thread_args[n].end = i + step > iterations ? iterations : i + step;
        int error = pthread_create(&threads[n],0,thread_func,(void*)&thread_args[n]);
        if(error)
        {
            printf("[ERROR %d] Couldn't create thread with tid=%d,start=%d,end=%d\n",error,thread_args[n].tid,thread_args[n].start,thread_args[n].end);
            exit(EXIT_FAILURE);
        }
    }
    for(int i=0; i<n; i++)
    {
        int error = pthread_join(threads[i],0);
        if(error)
        {
            printf("[ERROR %d] Couldn't join thread with tid=%d\n",error,i);
            exit(EXIT_FAILURE);
        }
    }
    if(result_out)
    {
        memcpy(result_out,recursive_reduce(reduce_func,(void*)results,0,n-1,result_element_size),result_element_size);
    }
}

typedef struct
{
    int dim;
    int** mat;
    int** ps;
} precomp_thread_args;

void* precomp_pthread_func(void* args)
{
    precomp_thread_args targs = get_data(precomp_thread_args,args);
    int** mat = targs.mat;
    int** ps = targs.ps;
    int dim = targs.dim;

    for(int j=get_start(args); j<get_end(args); j++)
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
    precomp_thread_args args;
    args.dim = dim;
    args.mat = mat;
    args.ps = ps;
    divide(dim,num_threads,precomp_pthread_func,args);
}

typedef struct
{
    int max_sum;
    int left;
    int right;
    int top;
    int bottom;
} maxarray_thread_ret;

typedef struct
{
    int dimension;
    int** ps;
    int* iter_i;
    int mat0;
} maxarray_thread_args;

void* maxArray_pthread_func(void* tdata)
{
    maxarray_thread_args targs = get_data(maxarray_thread_args,tdata);
    int tid = get_tid(tdata);
    int dim = targs.dimension;
    int start = get_start(tdata);
    int end = get_end(tdata);
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
            next_i = ++i<dim-1 ? iter_i[i+1] : end;

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

    maxarray_thread_ret ret;
    ret.max_sum = max_sum;
    ret.top = top;
    ret.left = left;
    ret.bottom = bottom;
    ret.right = right;

    write_result(maxarray_thread_ret,tdata,ret);
}

void* reduce_maxarray(void* d1, void* d2)
{
    if((*(maxarray_thread_ret*)d1).max_sum > (*(maxarray_thread_ret*)d2).max_sum) return d1;
    return d2;
}

maxarray_thread_ret maxarray_pthread(int** ps, int** mat, int dim, int num_threads)
{
    maxarray_thread_args args;

    int iter_i[dim];
    for(int i=0, sum=0; i<dim; sum+=dim-i,i++) iter_i[i] = sum;

    args.iter_i = iter_i;
    args.dimension = dim;
    args.ps = ps;
    args.mat0 = mat[0][0];

    maxarray_thread_ret ret;
    divide_reduce(dim * (dim + 1) / 2,num_threads,maxArray_pthread_func,args,reduce_maxarray,ret);
    return ret;
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

        maxarray_thread_ret ret = maxarray_pthread(ps, mat, dim, NUM_THREADS);
        printf("%d %d %d %d\n",ret.left,ret.top,ret.right,ret.bottom);


        /* release resources */
        fclose(input_file);
        free_matrix(mat, dim);
        free_matrix(ps, dim);
    }

    return EXIT_SUCCESS;
}
