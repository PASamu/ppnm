#ifndef HAVE_QR
#define HAVE_QR
void qr_decomp (gsl_matrix* A, gsl_matrix* R);
void qr_back (gsl_matrix* A, gsl_vector* R);
void qr_solve (gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector *x);
void qr_inv (gsl_matrix* Q, gsl_matrix* R, gsl_matrix *B);
#endif
