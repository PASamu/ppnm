#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#define RND ((double)rand()/RAND_MAX-0.5)*2

int jacobi_diag(gsl_matrix* A, gsl_matrix* V);

void matrix_print(char s[],gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("%10g",gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}

void matrix_sym(gsl_matrix* A){
	for(int i=0;i<A->size1;i++){
		for(int j=i;j<A->size2;j++){
			gsl_matrix_set(A,i,j,RND*20);
			gsl_matrix_set(A,j,i,gsl_matrix_get(A,i,j));
		}
	}
}

int main(){
	int n=4;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* C = gsl_matrix_alloc(n,n);
	matrix_sym(A);
	matrix_print("Randomly generated symmetric matrix A=",A);
	printf("\n");
	gsl_matrix_memcpy(Acopy,A);
	jacobi_diag(A,V);
	matrix_print("matrix A after Jacobi diagonalization:",A);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,Acopy,0,B);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,C);
	matrix_print("V^T*A*V (should equal A after Jacobi diag)=",C);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,A,0,B);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,B,V,0,C);
	matrix_print("V*D*V^T (should equal A before Jacobi diag)=",C);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,B);
	matrix_print("V^T*V (should be identity)=",B);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(Acopy);
	gsl_matrix_free(B);
	gsl_matrix_free(C);
}
	
	

