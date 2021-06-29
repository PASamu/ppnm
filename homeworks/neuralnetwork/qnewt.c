#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>


double delta=1.0/524288;

void num_grad (double f(gsl_vector*),gsl_vector* x, gsl_vector* grad){
	double fx=f(x);
	for(int i=0;i<x->size;i++){
		double xi=gsl_vector_get(x,i);
		double dx=fabs(xi)*delta;
		if(fabs(xi)<sqrt(delta))dx=delta;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(f(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int qnewt(double f(gsl_vector*),gsl_vector* x, double acc){
	int n=x->size,nsteps=0,bsteps=0,gsteps=0;
	
	gsl_matrix* B=gsl_matrix_alloc(n,n);
	gsl_vector* gx=gsl_vector_alloc(n);
	gsl_vector* Dx=gsl_vector_alloc(n);
	gsl_vector* z=gsl_vector_alloc(n);
	gsl_vector* gz=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
	gsl_vector* u=gsl_vector_alloc(n);

	gsl_matrix_set_identity(B);
	num_grad(f,x,gx);
	double fx=f(x),fz;
	while(nsteps<2000){
		nsteps++;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gx,0,Dx);
		if(gsl_blas_dnrm2(Dx)<delta*gsl_blas_dnrm2(x)){
			fprintf(stderr,"qnewton: |Dx|<delta*|x|\n");break;}
		if(gsl_blas_dnrm2(gx)<acc){
			fprintf(stderr,"qnewton: |grad|<acc\n");break;}
		double L=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			fz=f(z);
			double sTg; gsl_blas_ddot(Dx,gx,&sTg);
			if(fz<fx+0.01*sTg){
				gsteps++;
				break;
			}
			if(L<delta){
				bsteps++;
				break;
			}
			L*=0.5;
			gsl_vector_scale(Dx,0.5);
		}
		num_grad(f,z,gz);
		gsl_vector_memcpy(y,gz);
		gsl_blas_daxpy(-1,gx,y);
		gsl_vector_memcpy(u,Dx);
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
		double sTy,uTy;
		gsl_blas_ddot(Dx,y,&sTy);
		if(fabs(sTy)>1e-6){
			gsl_blas_ddot(u,y,&uTy);
			double gamma=uTy/2/sTy;
			gsl_blas_daxpy(-gamma,Dx,u);
			gsl_blas_dger(1.0/sTy,u,Dx,B);
			gsl_blas_dger(1.0/sTy,Dx,u,B);
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gx,gz);
		fx=fz;
	}
	gsl_matrix_free(B);
	gsl_vector_free(gx);
	gsl_vector_free(Dx);
	gsl_vector_free(z);
	gsl_vector_free(gz);
	gsl_vector_free(y);
	gsl_vector_free(u);
	fprintf(stderr,"qnewton: nsteps=%i gsteps=%i bsteps=%i fx=%.1e\n",nsteps,gsteps,bsteps,fx);
	return nsteps;
}
		
