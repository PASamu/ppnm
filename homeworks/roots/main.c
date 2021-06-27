#include<stdio.h>
#include<float.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

int newtonroots(void f(gsl_vector* x, gsl_vector* fx),
		gsl_vector* x, double eps);

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

void rosenbrock(gsl_vector* v, gsl_vector* fx){
	double x=gsl_vector_get(v,0), y=gsl_vector_get(v,1);
	
	//Gradient (partial derivates)
	gsl_vector_set(fx,0,-2.*(1-x)+4.*100.*x*(x*x-y));
	gsl_vector_set(fx,1,200.*(y-x*x));
}

int main(){
	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector_set(x,0,-1);
	gsl_vector_set(x,1,5);

	gsl_vector* fx = gsl_vector_alloc(2);
	printf("Finding extremum of Rosenbrock's function by finding roots of gradient\n");

	vector_print("Initial value of x=",x);
	rosenbrock(x,fx);
	vector_print("Initial value of gradient fx=",fx);

	double eps=1e-3;
	int steps=newtonroots(rosenbrock,x,eps);
	printf("Steps=%i\n",steps);
	vector_print("Extremum point:",x);
	rosenbrock(x,fx);
	vector_print("Gradient value at x (should be zero):",fx);
}

