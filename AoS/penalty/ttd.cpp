#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "algo.h"
#include "math.h"
#include <ctime>
#include <omp.h>
#define byteSize 32
#define nsize 10000

void penalty(
             double   xCoordinatesOfPointsOfGeometryA[3],
             double   yCoordinatesOfPointsOfGeometryA[3],
             double   zCoordinatesOfPointsOfGeometryA[3],
             double   xCoordinatesOfPointsOfGeometryB[3],
             double   yCoordinatesOfPointsOfGeometryB[3],
             double   zCoordinatesOfPointsOfGeometryB[3],
             double&  xPA,
             double&  yPA,
             double&  zPA,
             double&  xPB,
             double&  yPB,
             double&  zPB
             );

iREAL rndtriangles(unsigned n, iREAL length, iREAL a[][nsize], iREAL b[][nsize], iREAL c[][nsize], iREAL d[][nsize], iREAL e[][nsize], iREAL f[][nsize])
{
  int padding = 2*length;
  #pragma simd
  #pragma ivdep  
  for(unsigned i=0;i<n;i++)
  {
    a[0][i] = (length * drand48());
    a[1][i] = (length * drand48());
    a[2][i] = (length * drand48());

    b[0][i] = (length * drand48());
    b[1][i] = (length * drand48());
    b[2][i] = (length * drand48());

    c[0][i] = (length * drand48());
    c[1][i] = (length * drand48());
    c[2][i] = (length * drand48());

    d[0][i] = (length * drand48() + padding);
    d[1][i] = (length * drand48() + padding);
    d[2][i] = (length * drand48() + padding);

    e[0][i] = (length * drand48() + padding);
    e[1][i] = (length * drand48() + padding);
    e[2][i] = (length * drand48() + padding);

    f[0][i] = (length * drand48() + padding);
    f[1][i] = (length * drand48() + padding);
    f[2][i] = (length * drand48() + padding);
  }

  return 0.;
}


int main (int argc, char *argv[])
{
  int n = nsize;
  
  static iREAL a[3][nsize] __attribute__ ((aligned(byteSize)));
  static iREAL b[3][nsize] __attribute__ ((aligned(byteSize)));
  static iREAL c[3][nsize] __attribute__ ((aligned(byteSize)));
  static iREAL d[3][nsize] __attribute__ ((aligned(byteSize)));
  static iREAL e[3][nsize] __attribute__ ((aligned(byteSize)));
  static iREAL f[3][nsize] __attribute__ ((aligned(byteSize)));
  
  rndtriangles(n, 10, a, b, c, d, e, f);
  
  iREAL xCoordinatesOfPointsOfGeometryA[3] __attribute__ ((aligned(byteSize)));
  iREAL yCoordinatesOfPointsOfGeometryA[3] __attribute__ ((aligned(byteSize)));
  iREAL zCoordinatesOfPointsOfGeometryA[3] __attribute__ ((aligned(byteSize)));
  iREAL xCoordinatesOfPointsOfGeometryB[3] __attribute__ ((aligned(byteSize)));
  iREAL yCoordinatesOfPointsOfGeometryB[3] __attribute__ ((aligned(byteSize)));
  iREAL zCoordinatesOfPointsOfGeometryB[3] __attribute__ ((aligned(byteSize)));
  
  iREAL xPA __attribute__ ((aligned(byteSize)));
  iREAL yPA __attribute__ ((aligned(byteSize)));
  iREAL zPA __attribute__ ((aligned(byteSize)));
  iREAL xPB __attribute__ ((aligned(byteSize)));
  iREAL yPB __attribute__ ((aligned(byteSize)));
  iREAL zPB __attribute__ ((aligned(byteSize)));
  
  for (int i=0; i<n; i++)
  {
    xCoordinatesOfPointsOfGeometryA[0] = a[0][i];
    yCoordinatesOfPointsOfGeometryA[0] = a[1][i];
    zCoordinatesOfPointsOfGeometryA[0] = a[2][i];
    
    xCoordinatesOfPointsOfGeometryA[1] = b[0][i];
    yCoordinatesOfPointsOfGeometryA[1] = b[1][i];
    zCoordinatesOfPointsOfGeometryA[1] = b[2][i];
    
    xCoordinatesOfPointsOfGeometryA[2] = c[0][i];
    yCoordinatesOfPointsOfGeometryA[2] = c[1][i];
    zCoordinatesOfPointsOfGeometryA[2] = c[2][i];
    #pragma simd
    for (unsigned j=i+1; j<n; j++)
    {
      
      xCoordinatesOfPointsOfGeometryB[0] = d[0][j];
      yCoordinatesOfPointsOfGeometryB[0] = d[1][j];
      zCoordinatesOfPointsOfGeometryB[0] = d[2][j];
      
      xCoordinatesOfPointsOfGeometryB[1] = e[0][j];
      yCoordinatesOfPointsOfGeometryB[1] = e[1][j];
      zCoordinatesOfPointsOfGeometryB[1] = e[2][j];
      
      xCoordinatesOfPointsOfGeometryB[2] = f[0][j];
      yCoordinatesOfPointsOfGeometryB[2] = f[1][j];
      zCoordinatesOfPointsOfGeometryB[2] = f[2][j];
      
      penalty(\
         xCoordinatesOfPointsOfGeometryA,
         yCoordinatesOfPointsOfGeometryA,
         zCoordinatesOfPointsOfGeometryA,
         xCoordinatesOfPointsOfGeometryB,
         yCoordinatesOfPointsOfGeometryB,
         zCoordinatesOfPointsOfGeometryB,
         xPA, yPA, zPA, xPB, yPB, zPB);
    }
  }

  return 0;
}


/*A[0] = xCoordinatesOfPointsOfGeometryA[0]
 A[1] = yCoordinatesOfPointsOfGeometryA[0]
 A[2] = zCoordinatesOfPointsOfGeometryA[0]
 
 B[0] = xCoordinatesOfPointsOfGeometryA[1]
 B[1] = yCoordinatesOfPointsOfGeometryA[1]
 B[2] = zCoordinatesOfPointsOfGeometryA[1]
 
 C[0] = xCoordinatesOfPointsOfGeometryA[2]
 C[1] = yCoordinatesOfPointsOfGeometryA[2]
 C[2] = zCoordinatesOfPointsOfGeometryA[2]
 
 D[0] = xCoordinatesOfPointsOfGeometryB[0]
 D[1] = yCoordinatesOfPointsOfGeometryB[0]
 D[2] = zCoordinatesOfPointsOfGeometryB[0]
 
 E[0] = xCoordinatesOfPointsOfGeometryB[1]
 E[0] = yCoordinatesOfPointsOfGeometryB[1]
 E[0] = zCoordinatesOfPointsOfGeometryB[1]
 
 F[0] = xCoordinatesOfPointsOfGeometryB[2]
 F[0] = yCoordinatesOfPointsOfGeometryB[2]
 F[0] = zCoordinatesOfPointsOfGeometryB[2]*/

void penalty(
        double   xCoordinatesOfPointsOfGeometryA[3],
        double   yCoordinatesOfPointsOfGeometryA[3],
        double   zCoordinatesOfPointsOfGeometryA[3],
        double   xCoordinatesOfPointsOfGeometryB[3],
        double   yCoordinatesOfPointsOfGeometryB[3],
        double   zCoordinatesOfPointsOfGeometryB[3],
        double&  xPA,
        double&  yPA,
        double&  zPA,
        double&  xPB,
        double&  yPB,
        double&  zPB
        )
{
  iREAL BA[3], CA[3], ED[3], FD[3], hessian[16], x[4];

  BA[0] = xCoordinatesOfPointsOfGeometryA[1] - xCoordinatesOfPointsOfGeometryA[0];
  BA[1] = yCoordinatesOfPointsOfGeometryA[1] - yCoordinatesOfPointsOfGeometryA[0];
  BA[2] = zCoordinatesOfPointsOfGeometryA[1] - zCoordinatesOfPointsOfGeometryA[0];
  
  CA[0] = xCoordinatesOfPointsOfGeometryA[2] - xCoordinatesOfPointsOfGeometryA[0];
  CA[1] = yCoordinatesOfPointsOfGeometryA[2] - yCoordinatesOfPointsOfGeometryA[0];
  CA[2] = zCoordinatesOfPointsOfGeometryA[2] - zCoordinatesOfPointsOfGeometryA[0];
  
  ED[0] = xCoordinatesOfPointsOfGeometryB[1] - xCoordinatesOfPointsOfGeometryB[0];
  ED[1] = yCoordinatesOfPointsOfGeometryB[1] - yCoordinatesOfPointsOfGeometryB[0];
  ED[2] = zCoordinatesOfPointsOfGeometryB[1] - zCoordinatesOfPointsOfGeometryB[0];
  
  FD[0] = xCoordinatesOfPointsOfGeometryB[2] - xCoordinatesOfPointsOfGeometryB[0];
  FD[1] = yCoordinatesOfPointsOfGeometryB[2] - yCoordinatesOfPointsOfGeometryB[0];
  FD[2] = zCoordinatesOfPointsOfGeometryB[2] - zCoordinatesOfPointsOfGeometryB[0];
  
  hessian[0] = 2.*DOT(BA,BA);
  hessian[1] = 2.*DOT(CA,BA);
  hessian[2] = -2.*DOT(ED,BA);
  hessian[3] = -2.*DOT(FD,BA);
  
  hessian[4] = hessian[1]; //use symmetry
  hessian[5] = 2.*DOT(CA,CA);
  hessian[6] = -2.*DOT(ED,CA);
  hessian[7] = -2.*DOT(FD,CA);
  
  hessian[8] = hessian[2];
  hessian[9] = hessian[6];
  hessian[10] = 2.*DOT(ED,ED);
  hessian[11] = 2.*DOT(FD,ED);
  
  hessian[12] = hessian[3];
  hessian[13] = hessian[7];
  hessian[14] = hessian[11];
  hessian[15] = 2.*DOT(FD,FD);
  
  iREAL eps = 1E-16;
  iREAL delta = (hessian[0]+hessian[5]+hessian[10]+hessian[15]) * eps;
  iREAL lambda = sqrt(0.0125*(hessian[0]+hessian[5]+hessian[10]+hessian[15]));
  iREAL r = lambda*pow(10,log10(lambda)+8);
  
  #if iREAL==double
  iREAL tol = DBL_EPSILON;
  #else
  iREAL tol = FLT_EPSILON;
  #endif
  //initial guess
  x[0] = 0.33;
  x[1] = 0.33;
  x[2] = 0.33;
  x[3] = 0.33;
  
  //Declare loop variables;
  iREAL dx[4], a[16], SUBXY[3], b[4], mx[6], dh[8], tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  
  //Newton loop
  for(int i=0;i<4;i++)
  {
    dh[0] = dh[2] = dh[4] = dh[6] = -1;
    dh[1] = dh[3] = dh[5] = dh[7] = 1;
    
    if(-x[0] <= 0)
    {
        dh[0] = mx[0] = 0;
    }else
    {
        mx[0] = -x[0];
    }
    if(-x[1] <= 0)
    {
        dh[2] = mx[1] = 0;
    }else
    {
        mx[1] = -x[1];
    }
    if(x[0]+x[1]-1 <= 0)
    {
        dh[1] =  dh[3] = mx[2] = 0;
    }else
    {
        mx[2] = x[0]+x[1]-1;
    }
    if(-x[2] <= 0)
    {
        dh[4] = mx[3] = 0;
    }else
    {
        mx[3] = -x[2];
    }
    if(-x[3] <= 0)
    {
        dh[6] = mx[4] = 0;
    }else
    {
        mx[4] = -x[3];
    }
    if(x[2]+x[3]-1 <= 0)
    {
        dh[5] = dh[7] = 0;
        
        mx[5] = 0;
    }else
    {
        mx[5] = x[2]+x[3]-1;
    }
    
    delta = i < 4 ? delta : 2000*delta;
    
    SUBXY[0] = (xCoordinatesOfPointsOfGeometryA[0]+(BA[0] * x[0])+(CA[0] * x[1])) - (xCoordinatesOfPointsOfGeometryB[0]+(ED[0] * x[2])+(FD[0] * x[3]));
    SUBXY[1] = (yCoordinatesOfPointsOfGeometryA[0]+(BA[1] * x[0])+(CA[1] * x[1])) - (yCoordinatesOfPointsOfGeometryB[0]+(ED[1] * x[2])+(FD[1] * x[3]));
    SUBXY[2] = (zCoordinatesOfPointsOfGeometryA[0]+(BA[2] * x[0])+(CA[2] * x[1])) - (zCoordinatesOfPointsOfGeometryB[0]+(ED[2] * x[2])+(FD[2] * x[3]));
    
    b[0] = 2*DOT(SUBXY,BA) + r * (dh[0] * mx[0] + dh[1] * mx[2]);
    a[0] = hessian[0] + r * (dh[0] * dh[0] + dh[1] * dh[1]) + delta;
    a[4] = hessian[4] + r * (dh[3] * dh[1]);
    tmp1 = (hessian[1] + r * (dh[1] * dh[3]))/a[0];
    a[13] = hessian[13] - hessian[12] * tmp1;
    a[9] = hessian[9] - hessian[8] * tmp1;
    a[5] = (hessian[5] + r * (dh[2] * dh[2] + dh[3] * dh[3]) + delta) - a[4] * tmp1;
    b[1] = (2*DOT(SUBXY,CA) + r * (dh[2] * mx[1] + dh[3] * mx[2])) - b[0] * tmp1;
    tmp2 = hessian[2]/a[0];
    tmp3 = hessian[3]/a[0];
    tmp4 = ((hessian[6]) - a[4] * tmp2)/a[5];
    a[14] = ((hessian[14] + r * (dh[7] * dh[5])) - hessian[12] * tmp2) - a[13] * tmp4;
    a[10] = ((hessian[10] + r * (dh[4] * dh[4] + dh[5] * dh[5]) + delta) - hessian[8] * tmp2) - a[9] * tmp4;
    b[2] = ((-2*DOT(SUBXY,ED) + r * (dh[4] * mx[3] + dh[5] * mx[5])) - b[0] * tmp2) - b[1] * tmp4;
    tmp5 = (hessian[7] - a[4] * tmp3)/a[5];
    tmp6 = (((hessian[11] + r * (dh[5] * dh[7])) - hessian[8] * tmp3) - a[9] * tmp5)/a[10];
    
    dx[3] = ((((-2*DOT(SUBXY,FD) + r * (dh[6] * mx[4] + dh[7] * mx[5])) - b[2] * tmp6) - b[0] * tmp3) - b[1] * tmp5) / ((((hessian[15] + r * (dh[6] * dh[6] + dh[7] * dh[7]) + delta) - hessian[12] * tmp3) - a[13] * tmp5) - a[14] * tmp6);
    dx[2] = (b[2] - (a[14] * dx[3])) / a[10];
    dx[1] = (b[1] - (a[9] * dx[2] + a[13] * dx[3])) / a[5];
    dx[0] = (b[0] - (a[4] * dx[1] + hessian[8] * dx[2] + hessian[12] * dx[3])) / a[0];
    
    x[0] = x[0] - dx[0];
    x[1] = x[1] - dx[1];
    x[2] = x[2] - dx[2];
    x[3] = x[3] - dx[3];
  }
  
  xPA = xCoordinatesOfPointsOfGeometryA[0]+(BA[0] * x[0])+(CA[0] * x[1]);
  yPA = yCoordinatesOfPointsOfGeometryA[0]+(BA[1] * x[0])+(CA[1] * x[1]);
  zPA = zCoordinatesOfPointsOfGeometryA[0]+(BA[2] * x[0])+(CA[2] * x[1]);
  
  xPB = xCoordinatesOfPointsOfGeometryB[0]+(ED[0] * x[2])+(FD[0] * x[3]);
  yPB = yCoordinatesOfPointsOfGeometryB[0]+(ED[1] * x[2])+(FD[1] * x[3]);
  zPB = zCoordinatesOfPointsOfGeometryB[0]+(ED[2] * x[2])+(FD[2] * x[3]);
}
