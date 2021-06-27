#include<math.h>
#include<stdio.h>
#define FORLOOP(i) for(int i=0;i<n;i++)

void rkstep12(int n,void f(double t, double y[],double dydt[]),
		double t, double y[], double h,
		double yh[], double dy[]){
	double k0[n],y1[n],k1[n];
	
	f(t,y,k0); FORLOOP(i) y1[i]=y[i]+(0.5*h)*k0[i];
	f(t+0.5*h,y1,k1); FORLOOP(i) yh[i]=y[i]+h*k1[i];

       FORLOOP(i) dy[i]=(k1[i]-k0[i])*h*0.5;       
}

void driver12(int n,void f(double t, double y[],double dydt[]),
		double a, double y[], double b, 
		double h, double acc, double eps){
	double t=a;
	while(t<b){
		if(t+h>b){
			h=b-t;
		}
		double yh[n],dy[n];
		rkstep12(n,f,t,y,h,yh,dy);

		double sum=0;FORLOOP(i) sum+=y[i]*y[i];
		double norm_y=sqrt(sum);

		sum=0;FORLOOP(i)sum+=dy[i]*dy[i];
		double error=sqrt(sum);
		double tol=(acc+eps*norm_y)*sqrt(h/(b-a));

		if(error<tol){
			t=t+h;
			FORLOOP(i) y[i]=yh[i];
		}

		if(error>0){ 
			h*=0.95*pow(tol/error,0.25);
		}
		else{ 
			h*=2;
		}
	}
}

