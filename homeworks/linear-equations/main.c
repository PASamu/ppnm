#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX

void vector_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++){
		printf("%10g",gsl_vector_get(v,i));
		printf("\n");
	}
}

void matrix_print(char s[],gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("%10g",gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}
			
void matrix_gen(gsl_matrix* A){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			gsl_matrix_set(A,i,j,RND);
		}
	}
}

void vector_gen(gsl_vector* b){
	for(int i=0;i<b->size;i++){
		gsl_vector_set(b,i,RND);
	}
}

void gs_decomp(gsl_matrix* A, gsl_matrix* R){
	assert(A->size2==R->size1 && R->size2==R->size1);
	double normsq;
	double prod;
		for(int i=0;i<A->size2;i++){
		gsl_vector* ai = gsl_vector_alloc(A->size1);	//Taking columns in A as vectors
		gsl_matrix_get_col(ai,A,i); 
		gsl_blas_ddot(ai,ai,&normsq);
		gsl_matrix_set(R,i,i,sqrt(normsq)); 	//Finding Rii
		gsl_vector_scale(ai,(double)1/sqrt(normsq)); //Finding qi's
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<A->size2;j++){ 	//Time to find Rij. Recall that it is upper triangular
			gsl_vector* aj = gsl_vector_alloc(A->size1);
			gsl_matrix_get_col(aj,A,j);
			gsl_blas_ddot(ai,aj,&prod);
			gsl_vector_scale(aj, 1./prod);
			gsl_vector_sub(aj,ai);
			gsl_vector_scale(aj,prod);
			gsl_matrix_set_col(A,j,aj);
			gsl_matrix_set(R,i,j,prod);
			gsl_matrix_set(R,j,i,0);	//Zeros under the diagonal
			gsl_vector_free(aj);
		}
	gsl_vector_free(ai);
	}
}	

void backsub(gsl_matrix* U, gsl_vector* c){
	for(int i=c->size-1;i>=0;i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1; k<c->size; k++) s-=gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
		gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
	}
}

void gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	backsub(R,x);
}

int main(){
	int n=6,m=4;

	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	

	matrix_gen(A);
	matrix_print("Randomly generated matrix A = ",A);
	

	gs_decomp(A,R);
	matrix_print("Orthogonal matrix from decomposition Q =",A);
	printf("\n");
	matrix_print("Upper triangular matrix from decomposition R=",R);
	printf("\n");

	gsl_matrix* identity = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,identity);
	matrix_print("Q^T*Q (should be identity) = ", identity);
	printf("\n");
	
	gsl_matrix* QR = gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,QR);
	matrix_print("QR (should be A) = ",QR);
	printf("\n");
	
	

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(identity);
	gsl_matrix_free(QR);
	
	gsl_matrix* P = gsl_matrix_alloc(n,n);
	gsl_matrix* R2 = gsl_matrix_alloc(n,n);
	gsl_matrix* Pcopy = gsl_matrix_alloc(n,n);
	
	printf("We shall now solve Px=b for square P\n");
	matrix_gen(P);
	gsl_matrix_memcpy(Pcopy,P);
	matrix_print("Randomly generated square matrix P=",P);
	printf("\n");
	gs_decomp(P,R2);

	gsl_vector* b = gsl_vector_alloc(n);
	vector_gen(b);
	printf("\n");
	
	gsl_vector* x = gsl_vector_alloc(n);
	gs_solve(P,R2,b,x);
	gsl_vector* Px = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasNoTrans,1.0,Pcopy,x,0.0,Px);
	vector_print("x=",x);
	printf("\n");
	vector_print("b=",b);
	printf("\n");
	vector_print("Checking if Px=b, Px=",Px);

	gsl_matrix_free(P);
	gsl_matrix_free(R2);
	gsl_matrix_free(Pcopy);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(Px);
return 0;
}
