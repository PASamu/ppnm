#include<math.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

void num_grad(double f(gsl_vector*),gsl_vector* x, gsl_vector* grad);
int qnewt(double f(gsl_vector*),gsl_vector* x,double acc);

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

double rosenbrock(gsl_vector* q){
	double x=gsl_vector_get(q,0),y=gsl_vector_get(q,1);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

double himmelb(gsl_vector* q){
	double x=gsl_vector_get(q,0),y=gsl_vector_get(q,1);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

int main(){
	int n=2,nsteps;
	double rsb,him;
	double acc=1e-4;

	gsl_vector* q=gsl_vector_alloc(n);
	//gsl_vector* p=gsl_vector_alloc(n);
	
	printf("Minimization of Rosenbrock's function\n");
	//Initial values
	gsl_vector_set(q,0,2);
	gsl_vector_set(q,1,5);
	vector_print("Init point q1=",q);
	printf("Rsb(q1)=%g\n",rosenbrock(q));
	nsteps=qnewt(rosenbrock,q,acc);
	vector_print("Minimum point:",q);
	printf("Number of steps nsteps=%i\n",nsteps);
	printf("Value at minimum=%g\n",rosenbrock(q));

	printf("\nMinimization of Himmelblau's function\n");
	//Initial values
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,4);
	vector_print("Init point q1=",q);
	printf("Him(q1)=%g\n",himmelb(q));
	nsteps=qnewt(himmelb,q,acc);
	vector_print("Minimum point:",q);
	printf("Number of steps nsteps=%i\n",nsteps);
	printf("Value at minimum=%g\n",himmelb(q));

return 0;
}


