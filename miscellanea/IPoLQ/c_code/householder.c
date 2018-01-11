#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "matrix_func.h"




int main() {

  int i,j, maxeigloc;
  double a,x,y,z, maxeig;
  double diag[4], sdiag[4];
  dmat R, U ;
  int N = 4;
  R = CreateDmat(N,N);
  U = CreateDmat(3,3);
  //define the matrix values
  R.data[0] = 4;
  R.data[1] = 1;
  R.data[2] = -2;
  R.data[3] = 2;

  R.data[4] = 1;
  R.data[5] = 2;
  R.data[6] = 0;
  R.data[7] = 1;

  R.data[8] = -2;
  R.data[9] = 0;
  R.data[10] = 3;
  R.data[11] = -2;

  R.data[12] = 2;
  R.data[13] = 1;
  R.data[14] = -2;
  R.data[15] = -1;


  //static float R[4][4] = {4,1,-2,2,1,2,0,1,-2,0,3,-2,2,1,-2,1};

  printf("Tridiagonalization of matrix\n");
  //tridiagonalize the matrix
  TRED2(R.map, N, diag, sdiag);
  printf("Diagonal elements after tridiagonalization\n");
  for (i=0;i<=N;i++) {
    printf("%.4f\n" , diag[i]);
  }
  printf("Subdiagonal elements after tridiagonalization\n");
  for (i=0; i<=N; i++) {
    printf("%.4f\n", sdiag[i]);
  }
  /*printf("Matrix after tridiagonalization\n");
  for (i=0; i<=4; i++) {
    for (j=0;j<=4; j++) {
      printf("A[%d][%d] : %.4f\n", i,j,R.map[i][j]);}
    }*/
  printf("TQLI\n");
  //find eienvalues and vectors
  TQLI(diag, sdiag,N, R.map);

  printf("Eigenvalues\n");
  for (i=0; i<=N; i++) {

    printf("%.8f\n", diag[i]);
  }

  //find the highest eigenvalue
  maxeig = diag[0];
  maxeigloc = 0;
  for (i = 1; i < N; i++) {
    if (diag[i] > maxeig) {
      maxeig = diag[i];
      maxeigloc = i;
    }
  }
  //now create the rotation matrix
  a = R.data[maxeigloc];
  x = R.data[maxeigloc+4];
  y = R.data[maxeigloc+8];
  z = R.data[maxeigloc+12];

  // Construct the rotation matrix
  U.data[0] = a*a + x*x -y*y - z*z;
  U.data[1] = 2.0*(x*y + a*z);
  U.data[2] = 2.0*(z*x - a*y);
  U.data[3] = 2.0*(x*y - a*z);
  U.data[4] = a*a - x*x + y*y - z*z;
  U.data[5] = 2.0*(y*z + a*x);
  U.data[6] = 2.0*(z*x + a*y);
  U.data[7] = 2.0*(y*z - a*x);
  U.data[8] = a*a - x*x - y*y + z*z;
  printf("Rotation matrix U\n");
  for (i=0; i<9; i ++) {
    printf("%.4f\n", U.data[i]);
  }



}
