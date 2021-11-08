#include <assert.h>

/* errno */
#include <errno.h>

/* fopen, fscanf, fprintf, fclose */
#include <stdio.h>

/* EXIT_SUCCESS, EXIT_FAILURE, malloc, free */
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

static int create_mat(size_t const nrows, size_t const ncols, double ** const matp)
{
    double * mat=NULL;
    if (!(mat = (double*) malloc(nrows*ncols*sizeof(*mat)))) {
        goto cleanup;
    }

    /** Initialize matrix with random values **/
    for(size_t i = 0; i < nrows; i++){
        for (size_t j = 0; j < ncols; j++){
            mat[(i * ncols) + j] = (double)(rand() % 1000) / 353.0;
        }
    }
    /** End random initialization **/

    *matp = mat;

    return 0;

    cleanup:
    free(mat);
    return -1;
}

static int mult_mat(size_t const n, size_t const m, size_t const p,
                    double const * const A, double const * const B,
                    double ** const Cp)
{
  size_t i, j, k;
  double sum;
  double * C = NULL;
  int avail_threads;
  if (!(C = (double*) malloc(n*p*sizeof(*C)))) {
    goto cleanup;
  }

  //avail_threads = omp_get_num_threads();
  //printf("Available threads : %d", avail_threads);
  #pragma omp parallel 
  {
	  avail_threads = omp_get_num_threads();
	  printf("Available threads : %d", avail_threads);

  #pragma omp for collapse(2) schedule(dynamic) 
  
for (i=0; i<n; ++i) {
    for (j=0; j<p; ++j) {
      
      for (k=0, sum=0.0; k<m; ++k) {
        sum += A[i*m+k] * B[k*p+j];
      }
      C[i*p+j] = sum;
     }
  }
}

  *Cp = C;

  return 0;

  cleanup:
  free(C);

  /*failure:*/
  return -1;
}

void matrix_matrix_mult_tile(
  double ** dst, double * src1, double * src2,
  int nr, int nc, int nq,
  int rstart, int rend, int cstart, int cend,
  int qstart, int qend) {
  int r, c, q;
  double sum;

	double* F = NULL;
	  if (!(F = (double*) malloc(nr*nq*sizeof(*F)))) {
  //  goto cleanup;
  }

//  #pragma omp parallel for collapse(2) schedule(dynamic)
  for (r = rstart; r <= rend; r++) {
    for (c = cstart; c <= cend; c++) {
      if (qstart == 0) F[r*nc+c] = 0.0;
      for (q = qstart = 0; q <= qend; q++) {
        F[r*nc+c] += src1[r* nc + q] * src2[q * nq + c];
      } /* for q */
		} /* for c */
  } /* for r */
	
       *dst = F;
} /* matrix_matrix_mult_tile */

void matrix_matrix_mult_by_tiling(
  double ** dst, double * src1, double * src2,
  int nr, int nc, int nq,
  int rtilesize, int ctilesize, int qtilesize) {
  /* matrix_matrix_mult_by_tiling */
  int rstart, rend, cstart, cend, qstart, qend;

  double * C = NULL;

  if (!(C = (double*) malloc(nr*nq*sizeof(*C)))) {
    goto cleanup;
  }
  #pragma omp parallel for schedule(dynamic) 
  for (rstart = 0; rstart < nr; rstart += rtilesize) {
    rend = rstart + rtilesize + 1;
    if (rend >= nr) rend = nr - 1;
    for (cstart = 0; cstart < nc; cstart += ctilesize) {
      cend = cstart + ctilesize + 1;
      if (cend >= nc) cend = nc - 1;
      for (qstart = 0; qstart < nq; qstart += qtilesize) {
        qend = qstart + qtilesize + 1;
        if (qend >= nq) qend = nq - 1;
        matrix_matrix_mult_tile(&C, src1, src2, nr, nc, nq,
          rstart, rend, cstart, cend, qstart, qend);
	     } /* for qstart */
    } /* for cstart */
  
} /* for rstart */

	
	
  *dst = C;
	cleanup:
     	  free(C);
} /* matrix_matrix_mult_by_tiling */

int main(int argc, char * argv[])
{
  // size_t stored an unsigned integer
  size_t nrows, ncols, ncols2;
  double * A=NULL, * B=NULL, * C=NULL, * D = NULL;
  double startTime, endTime, startTime2, endTime2;
//	printf("Global variable is %f", x);	
	
  if (argc != 4) {
   fprintf(stderr, "usage: matmult nrows ncols ncols2\n");
   goto failure;
 }

 nrows = atoi(argv[1]);
 ncols = atoi(argv[2]);
 ncols2 = atoi(argv[3]);

 if (create_mat(nrows, ncols, &A)) {
   perror("error");
   goto failure;
 }

 if (create_mat(ncols, ncols2, &B)) {
   perror("error");
   goto failure;
 }
 startTime = omp_get_wtime();
 if (mult_mat(nrows, ncols, ncols2, A, B, &C)) {
   perror("error");
   goto failure;
 }
 endTime = omp_get_wtime();
printf("\n---------------------------------------- \n "); 
printf("\n Without tiling: %f \n ", (endTime - startTime) * pow(10,6));
 printf(" \n --------------------------------------- \n ");

 startTime2 = omp_get_wtime();
 matrix_matrix_mult_by_tiling(&D, A, B, nrows, ncols, ncols2, 200, 200, 200);
endTime2 = omp_get_wtime();

 printf(" \n --------------------------------------- \n ");
 printf("\n With tiling: %f \n", (endTime2 - startTime2) * pow(10,6));
 printf(" \n --------------------------------------- \n ");
 free(A);
 free(B);
 free(C);
 //free(D);
 return EXIT_SUCCESS;

 failure:
 if(A){
   free(A);
 }
 if(B){
   free(B);
 }
 if(C){
   free(C);
 }
 if(D){
  free(D);
 }
  return EXIT_FAILURE;
}





