#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "algo.h"
#include "math.h"
#include <ctime>
//#include <omp.h>
#define byteSize 32
#define nsize 10000000

REAL penalty(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3], REAL  P[3], REAL  Q[3]);

REAL pt(REAL TP1[3], REAL TP2[3], REAL TP3[3], REAL cPoint[3], REAL tq[3]);
REAL segseg(REAL p1[3], REAL p2[3], REAL p3[3], REAL p4[3], REAL P[3], REAL Q[3]);
REAL bf(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3], REAL  P[3], REAL  Q[3]);

REAL rndtriangles(unsigned n, REAL sidelgth, REAL a[][nsize], REAL b[][nsize], REAL c[][nsize], REAL d[][nsize], REAL e[][nsize], REAL f[][nsize])
{
  REAL inlo = 0;
  REAL inhi = sidelgth;
  REAL indiff = inhi-inlo;

  REAL inlo2 = sidelgth+1;
  REAL inhi2 = 2*sidelgth+1;
  REAL indiff2 = inhi2-inlo2;
  
  #pragma simd
  #pragma ivdep  
  for(unsigned i=0;i<n;i++)
  {
    a[0][i] = (indiff * ((REAL) (drand48())) + inlo);
    a[1][i] = (indiff * ((REAL) (drand48())) + inlo);
    a[2][i] = (indiff * ((REAL) (drand48())) + inlo);

    b[0][i] = (indiff * ((REAL) (drand48())) + inlo);
    b[1][i] = (indiff * ((REAL) (drand48())) + inlo);
    b[2][i] = (indiff * ((REAL) (drand48())) + inlo);

    c[0][i] = (indiff * ((REAL) (drand48())) + inlo);
    c[1][i] = (indiff * ((REAL) (drand48())) + inlo);
    c[2][i] = (indiff * ((REAL) (drand48())) + inlo);

    d[0][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    d[1][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    d[2][i] = (indiff2 * ((REAL) (drand48())) + inlo2);

    e[0][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    e[1][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    e[2][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    
    f[0][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    f[1][i] = (indiff2 * ((REAL) (drand48())) + inlo2);
    f[2][i] = (indiff2 * ((REAL) (drand48())) + inlo2); 
  }

  return 0.;
}

int main (int argc, char *argv[])
{
  unsigned n = nsize;

  static REAL a[3][nsize] __attribute__ ((aligned(byteSize)));  
  static REAL b[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL c[3][nsize]  __attribute__ ((aligned(byteSize)));
  static REAL d[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL e[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL f[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL p[3][nsize] __attribute__ ((aligned(byteSize))); 
  static REAL q[3][nsize] __attribute__ ((aligned(byteSize)));

  static REAL distance[nsize] __attribute__ ((aligned(byteSize))); 
   
  rndtriangles(n, 10, a, b, c, d, e, f);
  
  //omp_set_num_threads(8);

  REAL A[3] __attribute__ ((aligned(byteSize)));
  REAL B[3] __attribute__ ((aligned(byteSize)));
  REAL C[3] __attribute__ ((aligned(byteSize)));    
  REAL D[3] __attribute__ ((aligned(byteSize)));    
  REAL E[3] __attribute__ ((aligned(byteSize)));    
  REAL F[3] __attribute__ ((aligned(byteSize)));    
  REAL Q[3] __attribute__ ((aligned(byteSize)));    
  REAL P[3] __attribute__ ((aligned(byteSize)));    
  
  #pragma forceinline recursive
  #pragma simd
  //#pragma omp parallel for simd
  for(unsigned i=0;i<n;i++)
  {	
    A[0]=a[0][i];
    A[1]=a[1][i];
    A[2]=a[2][i];

    B[0]=b[0][i];
    B[1]=b[1][i];
    B[2]=b[2][i];

    C[0]=c[0][i];
    C[1]=c[1][i];
    C[2]=c[2][i];

    D[0]=d[0][i];
    D[1]=d[1][i];
    D[2]=d[2][i];

    E[0]=e[0][i];
    E[1]=e[1][i];
    E[2]=e[2][i];

    F[0]=f[0][i];
    F[1]=f[1][i];
    F[2]=f[2][i];
      
        
    p[0][i] = P[0];
    p[1][i] = P[1];
    p[2][i] = P[2];

    q[0][i] = Q[0];
    q[1][i] = Q[1];
    q[2][i] = Q[2];
  }

  return 0;
}
 
