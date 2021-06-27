#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

void qr_decomp(gsl_matrix* A, gsl_matrix* R){
	int m=A->size2;
	for(int i=0;i<m;i++){
		gsl_vector_view ai = gsl_matrix_column(A,i);
		double rii = gsl_blas_dnrm2(&ai.vector);
		gsl_matrix_set(R,i,i,rii);
		gsl_vector_scale(&ai.vector,1/rii);
		for(int j=i+1;j<m;j++){
			gsl_vector_view aj = gsl_matrix_column(A,j);
			double rij;
			gsl_blas_ddot(&ai.vector,&aj.vector,&rij);
			gsl_blas_daxpy(-rij,&ai.vector,&aj.vector);
			gsl_matrix_set(R,i,j,rij);
			gsl_matrix_set(R,j,i,0);
		}
	}
}

void qr_back(gsl_matrix *R,gsl_vector *x){
	int m=R->size1;
	for(int i=m-1;i>=0;i--){
		double s=gsl_vector_get(x,i);
		for(int k=i+1;k<m;k++)
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));		
	}
}

void qr_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	qr_back(R,x);
}

void qr_inv(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
	int n=R->size1;
	gsl_vector *b = gsl_vector_calloc(n);
	gsl_vector *x = gsl_vector_calloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(b,i,1.0);
		qr_solve(Q,R,b,x);
		gsl_vector_set(b,i,0.0);
		gsl_matrix_set_col(B,i,x);
	}
	gsl_vector_free(b);
	gsl_vector_free(x);
}
