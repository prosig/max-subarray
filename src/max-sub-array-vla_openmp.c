#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

//#define _PRINT_INFO 

static int thread_nr = 8;

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

/* FIXME - do we need parallelism here? is equally fast */
void 
clear(
	int* a, 
	int len) 
{
#if STOP_TIME
//	double start, time;
//	start = omp_get_wtime();
#endif /* STOP_TIME */

	/* precompute vertical prefix sum */

#if OMP
	/* FIXME - does not make sense to parallelize??? */
//	#pragma omp parallel for \
//		if(len > 512) \
//		schedule(static) \
//		num_threads(4)
#endif
#if 0
    for (int index=0; index<len; index+=4) {
		/* use loop unrolling */
		int i,j,k,l;
		
		i=index;
		j=index+1;
		k=index+2;
		l=index+3;
		
		a[i] = 0;
		a[j] = 0;
		a[k] = 0;
		a[l] = 0;
	}
#endif
memset((void*)a, 0, len); 

#if STOP_TIME
//	time = omp_get_wtime() - start;
//	printf("[TIME] - clear(): time needed %f sec\n", time); 
#endif /* STOP_TIME */
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
#if STOP_TIME
	double start, time;
	start = omp_get_wtime();
#endif /* STOP_TIME */

	/* precompute vertical prefix sum */
#if OMP
	#pragma omp parallel for \
		if(dim > 8) \
		schedule(static)  
//		num_threads(thread_nr)
#endif
	for(int j=0; j<dim; j++) {
        ps[0][j] = mat[0][j];
        for (int i=1; i<dim; i++) {
            ps[i][j] = ps[i-1][j] + mat[i][j];
        }
    }

#if STOP_TIME
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

    for(int i=0; i<dim; i++) {
#if OMP
		#pragma omp parallel for \
			schedule(static) \
			private(sum, pos, local_max) \
			shared(max_sum, top, bottom, left, right)
//			num_threads(8)
#endif
		for(int k=i; k<dim; k++) {
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
	
			/* use flush because it is faster than critical - this way
			 * the critical section will be accessed less often */
			#pragma omp flush (max_sum)
			if(sum[local_max] > max_sum) {
				#pragma omp critical
				if(sum[local_max] > max_sum) { 
                	/* sum[local_max] is the new max value
                	 * the corresponding submatrix goes from rows i..k.
                	 * and from columns pos[local_max]..local_max */
//printf("TEEEST hallo welt, wie geht es dir heute??? mir gehts sseeeeeehr guuuut!!! und was ist mit dir thread nur irgendwas anderes??? möchtest du mich unterbrechen???[THREAD=%d]\n", omp_get_thread_num());
                	max_sum = sum[local_max];
                	top = i;
                	left = pos[local_max];
                	bottom = k;
                	right = local_max;
				}
            }
        }
    }

    /* FIXME - Question: Do we need to compute the output matrix? */
	/* Compose the output matrix */

	int outmat_row_dim = bottom - top + 1;
    int outmat_col_dim = right - left + 1;
//	int xoutmat[outmat_row_dim][outmat_col_dim];
	outmat = alloc_matrix(outmat_row_dim, outmat_col_dim);

#if 0
int k=0;
	#pragma omp parallel for \
		schedule(static) \
		firstprivate(k)
//		num_threads(thread_nr)
    for(int i=top/*, k=0*/; i<=bottom; i++/*, k++*/) {
        for(int j=left, l=0; j<=right; j++, l++) {
            outmat[k][l] = mat[i][j];
        }
		k++;
    }
int l=0;
    for(int i=top, k=0; i<=bottom; i++, k++) {
		#pragma omp parallel for \
			schedule(static) \
			firstprivate(l)
//			num_threads(thread_nr)
        for(int j=left/*, l=0*/; j<=right; j++/*, l++*/) {
            outmat[k][l] = mat[i][j];
			l++;
        }
	}
#endif
	/* TODO - Question: how to parallelize??? maybe only a small nr of threads? */
    for(int i=top, k=0; i<=bottom; i++, k++) {
        for(int j=left, l=0; j<=right; j++, l++) {
            outmat[k][l] = mat[i][j];
        }
	}

    alg_end = get_usecs();

    /* print output matrix */
    printf("Sub-matrix [%dX%d] with max sum = %d, left = %d, top = %d,"
			" right = %d, bottom = %d\n", 
			outmat_row_dim, outmat_col_dim, max_sum, left, top, right, bottom);

#ifdef _PRINT_INFO
	printf("[INFO] output matrix: \n");
	print_matrix(outmat, outmat_row_dim, outmat_col_dim); 
#endif

    /* print stats */
    printf("%s,%f sec\n", "CHECK_NOT_PERFORMED", 
			((double)(alg_end-alg_start))/1000000);
	
	printf("[DEBUG] freeing outmat\n");
//	free_matrix(outmat, outmat_row_dim); 
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

	max_sub_arr(mat, ps, outmat, dim);

    /* release resources */
    fclose(input_file);
	printf("[DEBUG] freeing mat\n");
	free_matrix(mat, dim); 
	printf("[DEBUG] freeing ps\n");
	free_matrix(ps, dim); 

    return EXIT_SUCCESS;
}
