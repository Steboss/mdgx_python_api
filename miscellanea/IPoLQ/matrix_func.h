#ifndef MATRIXFUNC
#define MATRIXFUNC
#include "matrix.h"

#define SIGN2(a,b)         ((b) >= 0.0 ? fabs(a) : -fabs(a))

double pythag(double a, double b);

void TRED2(double** A, int n, double* d, double* e);

void TQLI(double* d, double* e, int n, double** z);

#endif
