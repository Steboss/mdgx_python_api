#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdgx.h"
#include "time.h"
#include "math.h"

//-----------------------------------------------------------------------------
// CpyDVec: copies a double-precision real vector V of length n into C; C is
//          allocated and returned.                                
//-----------------------------------------------------------------------------
double* CpyDVec(double* V, int n, double *C)
{
  int i;

  printf("Copying double precision vector V to vector C\n");

  //C = (double*)malloc(n*sizeof(double));
  for (i = 0; i < n; i++) {
    //printf("%f\n", V[i]);
    C[i] = V[i];
  }

  return C;
}