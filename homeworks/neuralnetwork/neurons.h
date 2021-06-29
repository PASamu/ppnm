#include"gsl/gsl_vector.h"
#ifndef HAVE_NEURONS
#define HAVE_NEURONS
typedef struct { int n; 
	double(*f)(double); 
	gsl_vector* params; 
	} ann;

ann* ann_alloc  (int n, double(*f)(double));
void ann_free (ann* network);
double response (ann* network, double x);
double cost_func (gsl_vector* p);
void train (ann* network, gsl_vector* xs, gsl_vector* ys);
#endif
