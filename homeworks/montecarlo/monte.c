#include<stdlib.h>
#include<math.h>
#include<complex.h>
#define RND (double)rand()/RAND_MAX

complex plainmonte(int dim, double f(double* x),
	       	double* a, double* b, int N){
	double V=1; for(int i=0;i<dim;i++) V*=b[i]-a[i];
	double sum=0,sum2=0,x[dim];
	for(int i=0;i<N;i++){
		for(int k=0;k<dim;k++) x[k]=a[k]+RND*(b[k]-a[k]);
		double fx=f(x);sum+=fx;sum2+=fx*fx;
	}
	double mean=sum/N,sigma=sqrt(sum2/N-mean*mean);
	complex result=mean*V+I*sigma*V/sqrt(N);
	return result;
}
