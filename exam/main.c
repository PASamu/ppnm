#include<stdio.h>
#include<math.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_blas.h>

int driver12(int m, void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt),
		gsl_vector_complex* tvals,gsl_matrix_complex* yvals,
		gsl_complex b, gsl_complex h,double acc, double eps);


//Differential equation with complex-valued function, that takes complex variable
void diffeq(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt){
	gsl_vector_complex_set(dydt,0,gsl_complex_mul_real(gsl_vector_complex_get(y,0),-1));
}

int main(){
	gsl_complex a,b,h,y0 = gsl_complex_rect(1,0);
	int n2=100;

	gsl_vector_complex* tvals=gsl_vector_complex_alloc(n2);
	gsl_matrix_complex* yvals=gsl_matrix_complex_alloc(1,n2);

	GSL_SET_COMPLEX(&a,0,0); 		//Start point of straight path	
	GSL_SET_COMPLEX(&b,3,0);		//End point of straight path
	GSL_SET_COMPLEX(&h,0.1,0.1);		//Step
	double acc=1e-2,eps=1e-2;
	gsl_vector_complex_set(tvals,0,a);
	gsl_matrix_complex_set(yvals,0,0,y0);

	int jj=driver12(n2,diffeq,tvals,yvals,b,h,acc,eps);
	FILE* out=fopen("out.txt","w");
	
	for(int i=0;i<jj+1;i++){
		fprintf(out,"%g %g %g %g\n",
				GSL_REAL(gsl_vector_complex_get(tvals,i)),
				GSL_IMAG(gsl_vector_complex_get(tvals,i)),
				GSL_REAL(gsl_matrix_complex_get(yvals,0,i)),
				GSL_IMAG(gsl_matrix_complex_get(yvals,0,i)));
	}
	fclose(out);
}
