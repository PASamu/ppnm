#include<math.h>
#include<stdio.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#define FORLOOP(i) for(int i=0;i<n;i++)

void rkstep12(void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt),
		gsl_complex t, gsl_vector_complex* y, gsl_complex h, gsl_vector_complex* yh, gsl_vector_complex* dy){
	int n=y->size;
	gsl_complex p,q,y1;
	gsl_vector_complex* k0=gsl_vector_complex_alloc(n);
	gsl_vector_complex* k12=gsl_vector_complex_alloc(n);
	gsl_vector_complex* vy=gsl_vector_complex_alloc(n);
	double habs=gsl_complex_abs(h);

	f(t,y,k0); 
	FORLOOP(i){p= gsl_complex_mul_real(gsl_vector_complex_get(k0,i),habs*1./2);
		y1=gsl_complex_add(gsl_vector_complex_get(y,i),p);
		gsl_vector_complex_set(vy,i,y1);}

        f(gsl_complex_add(t,gsl_complex_mul_real(h,1./2)),vy,k12);	
	FORLOOP(i){q=gsl_complex_mul_real(gsl_vector_complex_get(k12,i),habs);
		y1=gsl_complex_add(gsl_vector_complex_get(y,i),q);
		gsl_vector_complex_set(yh,i,y1);}


	FORLOOP(i){p=gsl_complex_sub(gsl_vector_complex_get(k12,i),gsl_vector_complex_get(k0,i)); 
		y1=gsl_complex_mul_real(p,1./2*habs);
		gsl_vector_complex_set(dy,i,y1);}

}

int driver12(int m,void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt),
		gsl_vector_complex* tvals,gsl_matrix_complex* yvals, gsl_complex b,
		gsl_complex h, double acc, double eps){

	const int n=yvals->size1;
	gsl_complex t,a=gsl_vector_complex_get(tvals,0);
	gsl_vector_complex* y=gsl_vector_complex_alloc(n);
	gsl_vector_complex* yh=gsl_vector_complex_alloc(n);
	gsl_vector_complex* dy=gsl_vector_complex_alloc(n);

	//Following part is to ensure a straight line
	double distAB=gsl_complex_abs(gsl_complex_sub(b,a));
	double realL=GSL_REAL(b)-GSL_REAL(a);
	double imagL=GSL_REAL(b)-GSL_REAL(a);
	double p=GSL_REAL(h)/realL,q=GSL_IMAG(h)/imagL;
	double r;

	if(p<q) r=p; else r=q;
	h=gsl_complex_rect(r*realL,r*imagL);
	double dist=0;
	int j=0;
	while(dist<distAB){	//In the real case it was while(t<b). Adapted to complex by using absolute values.
		t=gsl_vector_complex_get(tvals,j);
		gsl_matrix_complex_get_col(y,yvals,j);
		if((dist+gsl_complex_abs(h))>distAB){
			h=gsl_complex_sub(b,t);
		}
		
		rkstep12(f,t,y,h,yh,dy);

	
		double norm_y=gsl_blas_dznrm2(yh);
		double error=gsl_blas_dznrm2(dy);
		
		double tol=(acc+eps*norm_y)*sqrt(gsl_complex_abs(h)/distAB);
		
		if(error<tol){
			j++;
			dist=dist+gsl_complex_abs(h);
			gsl_vector_complex_set(tvals,j,gsl_complex_add(t,h));
			gsl_matrix_complex_set_col(yvals,j,yh);
		}

		h=gsl_complex_mul_real(h,0.95*pow(tol/error,0.25));
		
		
	}
	gsl_vector_complex_free(y);
	gsl_vector_complex_free(yh);
	gsl_vector_complex_free(dy);
return j;
}
