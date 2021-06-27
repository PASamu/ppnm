#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<assert.h>

void GS(gsl_matrix* A, gsl_matrix* R){
	assert(A->size2==R->size1 && R->size2==R->size1);
	double normsq;
	double prod;
	for(int i=0;i<A->size2;i++){
		gsl_vector* ai = gsl_vector_alloc(A->size1);
		gsl_matrix_get_col(ai,A,i);
		gsl_blas_ddot(ai,ai,&normsq);
		gsl_matrix_set(R,i,i,sqrt(normsq));
		gsl_vector_scale(ai,(double)1/sqrt(normsq));
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<A->size2;j++){
			gsl_vector* aj=gsl_vector_alloc(A->size1);
			gsl_matrix_get_col(aj,A,j);
			gsl_blas_ddot(ai,aj,&prod);
			gsl_vector_scale(aj,1./prod);
			gsl_vector_sub(aj,ai);
			gsl_vector_scale(aj,prod);
			gsl_matrix_set_col(A,j,aj);
			gsl_matrix_set(R,i,j,prod);
			gsl_matrix_set(R,j,i,0);
			gsl_vector_free(aj);
		}
		gsl_vector_free(ai);
	}
}

void backsub(gsl_matrix* U, gsl_vector* c){
	for(int i=c->size-1;i>=0;i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k<c->size;k++) s-=gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	backsub(R,x);
}
