#ifndef MATRIXSTRUC
#define MATRIXSTRUC


struct DMatrix {
  int row;
  int col;
  double* data;
  double** map;
};
typedef struct DMatrix dmat;

dmat CreateDmat(int M, int N);

#endif
