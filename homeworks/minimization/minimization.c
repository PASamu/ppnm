#include<math.h>
#include<stdio.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

double delta=sqrt(DBL_EPSILON);

void num_grad(double f(gsl_vector*),gsl_vector* x, gsl_vector* grad){
	double fx=f(x);
	for(int i=0;i<x->size;i++){
		double dx,xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(delta)){
			dx=delta;
		}
		else{
			dx=fabs(xi)*delta;
		}
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(f(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int qnewt(double f(gsl_vector*),gsl_vector* x, double acc){
	int n=x->size,nsteps=0,bsteps=0,gsteps=0;

	gsl_matrix* B 	 = gsl_matrix_alloc(n,n);
	gsl_vector* grad = gsl_vector_alloc(n);
	gsl_vector* Dx 	 = gsl_vector_alloc(n);
	gsl_vector* y 	 = gsl_vector_alloc(n);
	gsl_vector* gy 	 = gsl_vector_alloc(n);
	gsl_vector* z 	 = gsl_vector_alloc(n);
	gsl_vector* u 	 = gsl_vector_alloc(n);
	gsl_vector* a 	 = gsl_vector_alloc(n);

	gsl_matrix_set_identity(B);
	num_grad(f,x,grad);
	double fx=f(x),fy;
	while(nsteps<1500){
		nsteps++;
		gsl_blas_dgemv(CblasNoTrans,-1.0,B,grad,0,Dx);
		if(gsl_blas_dnrm2(Dx)<delta*gsl_blas_dnrm2(x)){
			fprintf(stderr,"qnewton: |Dx|<delta*|x|\n"); break;
		}
		if(gsl_blas_dnrm2(grad)<acc){
			fprintf(stderr,"qnewton: |grad|<acc\n"); break;
		}
		double L=1;
		while(1){
			gsl_vector_memcpy(y,x);
			gsl_vector_add(y,Dx);
			fy=f(y);
			double sTg; gsl_blas_ddot(Dx,grad,&sTg);
			if(fy<fx+0.01*sTg){
				gsteps++;break;
			}
			if(L<delta){
				bsteps++;
				gsl_matrix_set_identity(B);
				break;
			}
			L*=0.5;
			gsl_vector_scale(Dx,0.5);
		}
		num_grad(f,y,gy);
		gsl_vector_memcpy(z,gy);
		gsl_blas_daxpy(-1.0,grad,z);
		gsl_vector_memcpy(u,Dx);
		gsl_blas_dgemv(CblasNoTrans,-1.0,B,z,1.0,u);
		double sTy,uTy;
		gsl_blas_ddot(Dx,z,&sTy);
		if(fabs(sTy)>1e-10){
			gsl_blas_ddot(u,z,&uTy);
			double gam=uTy/2/sTy;
			gsl_blas_daxpy(-gam,Dx,u);
			gsl_blas_dger(1.0/sTy,u,Dx,B);
			gsl_blas_dger(1.0/sTy,Dx,u,B);
		}
		gsl_vector_memcpy(x,y);
		gsl_vector_memcpy(grad,gy);
		fx=fy;
	}
	gsl_matrix_free(B);
	gsl_vector_free(grad);
	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(gy);
	gsl_vector_free(z);
	gsl_vector_free(a);
	gsl_vector_free(u);
	fprintf(stderr,"qnewton: nsteps=%i, gsteps=%i, bsteps=%i fx=%.1e\n",nsteps,gsteps,bsteps,fx);
	return nsteps;
}

