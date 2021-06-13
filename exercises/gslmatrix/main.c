#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g",gsl_vector_get(v,i));
	printf("\n");
}

int main(){
	int n=3;
TRACE("n=%i\n",n);
TRACE("1\n");
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
TRACE("2\n");
	double vec[3]={6.23, 5.37, 2.29};
	double mat[9]={6.13, -2.90, 5.86, 8.08, -6.31, -3.89, -4.36, 1.00, 0.19};
TRACE("3\n");
	gsl_matrix_set(A,0,0,mat[0]);
	gsl_matrix_set(A,0,1,mat[1]);
	gsl_matrix_set(A,0,2,mat[2]);
	gsl_matrix_set(A,1,0,mat[3]);
	gsl_matrix_set(A,1,1,mat[4]);
	gsl_matrix_set(A,1,2,mat[5]);
	gsl_matrix_set(A,2,0,mat[6]);
	gsl_matrix_set(A,2,1,mat[7]);
	gsl_matrix_set(A,2,2,mat[8]);
TRACE("4\n");
	gsl_matrix_memcpy(Acopy,A);
TRACE("5\n");
	gsl_vector_set(b,0,vec[0]);
	gsl_vector_set(b,1,vec[1]);
	gsl_vector_set(b,2,vec[2]);
TRACE("6\n");
	gsl_linalg_HH_solve(Acopy,b,x);
TRACE("7\n");
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
TRACE("8\n");
	vector_print("Solution to system of linear equations:",x);
	vector_print("check: A*x should be equal b:",y);
TRACE("9\n");
return 0;
}
