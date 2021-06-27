#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix.h>
#include"qr.h"

int lim=5000;

int newtonroots(void f(gsl_vector* x, gsl_vector* fx),gsl_vector* x, double eps){

	int n=x->size;
	int steps=0;
	double delta=sqrt(DBL_EPSILON);

	gsl_matrix* J  = gsl_matrix_alloc(n,n);
	gsl_matrix* R  = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* y  = gsl_vector_alloc(n);
	gsl_vector* fy = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	while(steps<lim){
		f(x,fx);
		for(int j=0;j<n;j++){
			double xj=gsl_vector_get(x,j);
			gsl_vector_set(x,j,xj+delta);
			f(x,df);
			gsl_vector_sub(df,fx);
			for(int i=0;i<n;i++){
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/delta);
			}
			gsl_vector_set(x,j,xj);
		}
		GS(J,R);
		GS_solve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1.0);
		double s=1.0;
		while(s>1./64){
			gsl_vector_memcpy(y,x);
			gsl_vector_add(y,Dx);
			f(y,fy);
			if(gsl_blas_dnrm2(fy)<(1-s/2)*gsl_blas_dnrm2(fx)) break;
			s*=0.5;
			gsl_vector_scale(Dx,0.5);
		}
		gsl_vector_memcpy(x,y);
		gsl_vector_memcpy(fx,fy);
		if(gsl_blas_dnrm2(Dx)<delta || gsl_blas_dnrm2(fx)<eps){
		       	break;
		}
		steps++;
       	}
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(y);
	gsl_vector_free(fy);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
return steps;
}
