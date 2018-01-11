#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_func.h"


//accessories
double pythag(double a, double b)
{
  double absa, absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
    return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  }
  else {
    return absb*sqrt(1.0+(absa/absb)*(absa/absb));
  }
}
//
void TRED2(double** A, int n, double* d, double* e)
{
  int i, j, k, l;
  double scale, hh, h, g, f;
  double* tmp;
  double* tm2p;

  for (i = n-1; i >= 1; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    tmp = A[i];
    if (l > 0) {
      for (k = 0; k <= l; k++) {
	scale += fabs(tmp[k]);
      }
      if (scale == 0.0) {
	e[i] = tmp[l];
      }
      else {
	for (k = 0; k <= l; k++) {
	  tmp[k] /= scale;
	  h += tmp[k]*tmp[k];
	}
	f = tmp[l];
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i] = scale*g;
	h -= f*g;
	tmp[l] = f - g;
	f = 0.0;
	for (j = 0; j <= l; j++) {
	  tm2p = A[j];
	  tm2p[i] = tmp[j]/h;
	  g = 0.0;
	  for (k = 0; k <= j; k++) {
	    g += tm2p[k]*tmp[k];
	  }
	  for (k = j+1; k <=l; k++) {
	    g += A[k][j]*tmp[k];
	  }
	  e[j] = g/h;
	  f += e[j]*tmp[j];
	}
	hh = f/(h + h);
	for (j = 0; j <= l; j++) {
	  f = tmp[j];
	  e[j] = g = e[j] - hh*f;
	  tm2p = A[j];
	  for (k = 0; k <= j; k++) {
	    tm2p[k] -= (f*e[k] + g*tmp[k]);
	  }
	}
      }
    }
    else {
      e[i] = tmp[l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    tmp = A[i];
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) {
          g += tmp[k]*A[k][j];
        }
        for (k = 0; k <= l; k++) {
          A[k][j] -= g*A[k][i];
        }
      }
    }
    d[i] = tmp[i];
    tmp[i] = 1.0;
    for (j = 0; j <= l; j++) {
      A[j][i] = tmp[j] = 0.0;
    }
  }
}

//
void TQLI(double* d, double* e, int n, double** z)
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  double* tmp;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;

    // Added to initialize m
    m = l - 1;

    while (m != l) {
      for (m = l; m < n - 1; m++) {
	dd = fabs(d[m]) + fabs(d[m+1]);
	if (fabs(e[m]+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  printf("TQLI >> Error: Too many iterations in tqli\n");
	}
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = pythag(g,1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN2(r,g));
	c = 1.0;
	s = 1.0;
	p = 0.0;
	for (i = m-1; i >= l; i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i] - g)*s + 2.0*c*b;
	  d[i+1] = g + (p = s*r);
	  g = c*r - b;
	  for (k = 0; k < n; k++) {
	    tmp = z[k];
	    f = tmp[i+1];
	    tmp[i+1] = s*tmp[i] + c*f;
	    tmp[i] = c*tmp[i] - s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    }
  }
}
