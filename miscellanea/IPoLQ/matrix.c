#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

dmat CreateDmat(int M, int N)
{
  int i;
  dmat A;

  // Allocate memory
  A.data = (double*)calloc(M*N, sizeof(double));

  // Make the map
  A.map = (double**)malloc(M*sizeof(double*));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[i*N];
  }

  // Store the dimensionality
  A.row = M;
  A.col = N;

  return A;
}
