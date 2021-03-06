#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

//#define	_PRINT_INFO 
//#define PRINT_RESULT		1
//#define STOP_TIME			1

#define	NUM_ROUNDS			10
#define	MAX_THREADS			48

static int thread_nr = 8;
static double runtime = 0;

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
    printf("Argument error! Usage: %s <input_file> \n", app_name);
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

	if(!matrix) {
		printf("[ERROR] could not allocate memory for matrix\n");
		exit(EXIT_FAILURE);
	} else {
		for (int i=0; i<n; i++) {
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
	if(!matrix) {
		printf("[ERROR] could not free memory for matrix - already freed\n");
		exit(EXIT_FAILURE);
	} else {
		for (int i=0; i<n; i++) {
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

	for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            if (j != (n-1)) { 
                fscanf(input_file, "%d\t", &element);
			} else { 
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
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void 
precomp_matrix(
	int **mat, 
	int **ps, 
	int dim) 
{
#ifdef STOP_TIME
	double start, time;
	start = omp_get_wtime();
#endif /* STOP_TIME */

	/* precompute vertical prefix sum */
	#pragma omp parallel for \
		default(none) \
		firstprivate(dim) \
		shared(ps, mat) \
		if(dim > 8) \
		schedule(static) \
		num_threads(thread_nr)
	for(int j=0; j<dim; j++) {
        ps[0][j] = mat[0][j];
        for (int i=1; i<dim; i++) {
            ps[i][j] = ps[i-1][j] + mat[i][j];
        }
    }

#ifdef STOP_TIME
	time = omp_get_wtime() - start;
	printf("[TIME] - precomp_matrix(): time needed %f sec\n", time); 
#endif /* STOP_TIME */
}

void 
max_sub_arr(
	int **mat, 
	int **ps, 
	int **outmat, 
	int dim) 
{
	int max_sum = mat[0][0];
    int top = 0, left = 0, bottom = 0, right = 0; 
	long alg_start, alg_end;

    /* auxilliary variables */
    int sum[dim];
    int pos[dim];
    int local_max;

    alg_start = get_usecs();
    
	/* precompute vertical prefix sum */
	precomp_matrix(mat, ps, dim);

#ifdef _PRINT_INFO
    /* Print the matrix */
    printf("[INFO] Vertical prefix sum matrix [%d]\n", dim);
	print_matrix(ps, dim, dim); 
#endif

	thread_ret tr[thread_nr];
	memset((void*)tr, 0, sizeof(thread_ret) * thread_nr);
	/*
	for(int t=0; t<thread_nr; t++) {
		tr[t].max_sum = 0;
	}
	*/

    for(int i=0; i<dim; i++) {
		#pragma omp parallel for \
			default(none) \
			schedule(static) \
			firstprivate(dim, i) \
			private(sum, pos, local_max) \
			shared(ps, tr) \
			num_threads(thread_nr)
/*			firstprivate(dim, i) \
			private(sum, pos, local_max) \
			shared(max_sum, top, bottom, left, right, ps) \
*/
		for(int k=i; k<dim; k++) {
			int tid = omp_get_thread_num();

            /* Kandane over all columns with the i..k rows */
            clear(sum, dim);
            clear(pos, dim);
            local_max = 0;

			/* we keep track of the position of the max value over each 
			 * Kandane's execution 
             * -> Notice that we do not keep track of the max value, 
			 *    but only its position */
            sum[0] = ps[k][0] - (i==0 ? 0 : ps[i-1][0]);
            for(int j=1; j<dim; j++) {
                if(sum[j-1] > 0){
                    sum[j] = sum[j-1] + ps[k][j] - (i==0 ? 0 : ps[i-1][j]);
                    pos[j] = pos[j-1];
                } 
                else {
                    sum[j] = ps[k][j] - (i==0 ? 0 : ps[i-1][j]);
                    pos[j] = j;
                }
                if(sum[j] > sum[local_max]) {
                    local_max = j;
                }
            } /* Kandane ends here */

			if(sum[local_max] > tr[tid].max_sum) {
                /* sum[local_max] is the new max value
                 * the corresponding submatrix goes from rows i..k.
                 * and from columns pos[local_max]..local_max */
                tr[tid].max_sum = sum[local_max];
                tr[tid].top = i;
                tr[tid].left = pos[local_max];
                tr[tid].bottom = k;
                tr[tid].right = local_max;
			}


/* WARNING: old implementation: DO NOT DELETE */
#if 0
			/* use flush because it is faster than critical - this way
			 * the critical section will be accessed less often */
			#pragma omp flush (max_sum)
			if(sum[local_max] > max_sum) {
				#pragma omp critical
				if(sum[local_max] > max_sum) { 
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
#endif
		}

		/* executed sequentially */
		for(int t=0; t<thread_nr; t++) {
			if(tr[t].max_sum > max_sum) {
               	max_sum = tr[t].max_sum;
                top = tr[t].top;
                left = tr[t].left;
				bottom = tr[t].bottom;
                right = tr[t].right;
			}
		}
    }

	int outmat_row_dim = bottom - top + 1;
    int outmat_col_dim = right - left + 1;
	outmat = alloc_matrix(outmat_row_dim, outmat_col_dim);

#ifdef STOP_TIME
	double start, time;
	start = omp_get_wtime();
#endif /* STOP_TIME */

    /* TODO - Question: Do we need to compute the output matrix? */
	/* Compose the output matrix */
	/* TODO - Question: how to parallelize??? does it make sense to parallelize
	 * copying operations??? */

	/* calculate output matrix */
    for(int i=top, k=0; i<=bottom; i++, k++) {
		int j=left;

		memcpy((void*) outmat[k], (void*)(&mat[i][j]), 
				sizeof(**mat)*(right-left+1));
	}

#ifdef STOP_TIME
	time = omp_get_wtime() - start;
	printf("[TIME] - outmat(): time needed %f sec\n", time); 
#endif /* STOP_TIME */

    alg_end = get_usecs();
	runtime = (double)(alg_end-alg_start)/1000000;

#ifdef PRINT_RESULT
    printf("[RESULT] Threads=%d sub-matrix [%dX%d] in %f sec\n", 
			thread_nr, outmat_row_dim, outmat_col_dim, runtime);
	printf(" -> max sum=%d, left=%d, top=%d right=%d, bottom=%d\n", 
			 max_sum, left, top, right, bottom);
#endif

#ifdef _PRINT_INFO
    /* print output matrix */
	printf("[INFO] output matrix: \n");
	print_matrix(outmat, outmat_row_dim, outmat_col_dim); 
#endif

	free_matrix(outmat, outmat_row_dim); 
}

int 
main(int argc, char* argv[]) 
{
	int dim = 0;
	int **mat=NULL, **ps=NULL, **outmat=NULL;
	FILE* input_file;
    
	if(argc != 2) {
        usage(argv[0]);
    }

    /* open files */
    input_file = fopen(argv[1], "r");
    if(input_file == NULL) {
        usage(argv[0]);
    }

    /* read matrix dimension */
    dim = get_mat_dim(input_file);  
	
	/* init matrices */
	mat = alloc_matrix(dim, dim); 
	ps = alloc_matrix(dim, dim);

	/* read matrix */
	read_matrix(mat, dim, dim, input_file); 

#ifdef _PRINT_INFO
    /* Print the matrix */
    printf("[INFO] Input matrix [%d]\n", dim);
	print_matrix(mat, dim, dim); 
#endif

	printf("Evaluation:\n");
	printf("===========\n");
	
	/* evaluate the program using multiple number of threads/cores */
	for(int i=1; i<=MAX_THREADS; i*=2) {
		thread_nr = i;
		double avg_time = 0;

		/* compute an average value of the time needed to run the program */
		for(int j=0; j<NUM_ROUNDS; j++) {
			max_sub_arr(mat, ps, outmat, dim);
			/* runtime of the algorithm is modified by the algorithm itself */
			avg_time += runtime; 
		}

		avg_time /= NUM_ROUNDS;
    
		/* print stats */
    	printf("[STAT] Threads=%d: runtime=%f sec\n", thread_nr, avg_time);
	}

	/* last iteration */
	thread_nr = MAX_THREADS;
	double avg_time = 0;
	/* compute an average value of the time needed to run the program */
	for(int j=0; j<NUM_ROUNDS; j++) {
		max_sub_arr(mat, ps, outmat, dim);
		/* runtime of the algorithm is modified by the algorithm itself */
		avg_time += runtime; 
	}

	avg_time /= NUM_ROUNDS;
    
	/* print stats */
    printf("[STAT] Threads=%d: runtime=%f sec\n", thread_nr, avg_time);




    /* release resources */
    fclose(input_file);
	free_matrix(mat, dim); 
	free_matrix(ps, dim); 

    return EXIT_SUCCESS;
}
