#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"QR_decomp.h"

void lsfit(int m, double f(int i, double x),
		gsl_vector* x, gsl_vector* y, gsl_vector* dy,
		gsl_vector* c, gsl_matrix* S);

int main(){
	double x[]={1,2,3,4,6,9,10,13,15};
	double y[]={117,100,88,72,53,29.5,25.2,15.2,11.1};

	int n=sizeof(x)/sizeof(x[0]);
	double deltay[n];
	for(int i=0;i<n;i++){
		deltay[i]=0.05*y[i];
		}

	for(int i=0;i<n;i++){ //Finding log of y and rewriting the uncertainty accordingly
		deltay[i]/=y[i];
		y[i]=log(y[i]);
	}

	gsl_vector* vx = gsl_vector_alloc(n);
	gsl_vector* vy = gsl_vector_alloc(n);
	gsl_vector* vdeltay = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){ //Putting data into vectors
		gsl_vector_set(vx,i,x[i]);
		gsl_vector_set(vy,i,y[i]);
		gsl_vector_set(vdeltay,i,deltay[i]);
	}

	double func(int i, double t){
		switch(i){
			case 0: return 1; break;
			case 1: return t; break;
			default: return NAN;
		}
	}

	int m=2;
	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	lsfit(m,func,vx,vy,vdeltay,c,S);


	double fit(double x){
		double s=0;
		for(int k=0;k<m;k++)s+=gsl_vector_get(c,k)*func(k,x);
		return s;
	}

	gsl_vector* dc = gsl_vector_alloc(m);
	for(int k=0;k<m;k++){
		double skk=gsl_matrix_get(S,k,k);
		gsl_vector_set(dc,k,sqrt(skk));
	}

	double fit_plus(int i, double x){
		return fit(x)+gsl_vector_get(dc,i)*func(i,x);
	}
	
	double fit_minus(int i, double x){
		return fit(x)-gsl_vector_get(dc,i)*func(i,x);
	}

	double c1 = gsl_vector_get(c,1);
	double dc1 = gsl_vector_get(dc,1);
	double T=-1/c1*log(2.0); //Halftime
	double dT = dc1/c1/c1;
	
	
	printf("# Half-life (fit) = %.3g +- %.2g days\n",T,dT);
	printf("# Half-life (table value) is approx 3.6 days.\n");
	printf("# Time Log(activity) delta(Log(activity))\n");
	


	for(int i=0;i<n;i++)printf("%g %g %g\n",x[i], y[i],deltay[i]);
	printf("\n\n");


	
	for(int i=0;i<m;i++){
		printf("# Time Fit Fit_plus Fit_minus;k=%i\n",i);
		for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz)
			printf("%g %g %g %g\n",z,fit(z),fit_plus(i,z),fit_minus(i,z));
		printf("\n\n");
	}
	
	
	gsl_vector_free(vx);
	gsl_vector_free(vy);
	gsl_vector_free(vdeltay);
	gsl_vector_free(c);
	gsl_matrix_free(S);
	gsl_vector_free(dc);

return 0;
}
