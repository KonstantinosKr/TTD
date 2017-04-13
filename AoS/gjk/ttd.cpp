#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "algo.h"
#include "math.h"
#include "gjk.h"
#include <ctime>
#define byteSize 32
#define nsize 10000


iREAL  rndtriangles(int n, iREAL  length, iREAL  a[nsize], iREAL  b[nsize], iREAL  c[nsize], iREAL  d[nsize], iREAL  e[nsize], iREAL  f[nsize])
{
  int padding = 2*length;
  #pragma simd
  #pragma ivdep
  for(int i=0;i<n;i++)
  {
    a[i] = (length * drand48());

    b[i] = (length * drand48());

    c[i] = (length * drand48());

    d[i] = (length * drand48() + padding);

    e[i] = (length * drand48() + padding);

    f[i] = (length * drand48() + padding);
  }

  return 0.;
}


int main (int argc, char *argv[])
{
  int n = nsize;

  static iREAL a[nsize] __attribute__ ((aligned(byteSize)));
  static iREAL b[nsize] __attribute__ ((aligned(byteSize)));
  static iREAL c[nsize] __attribute__ ((aligned(byteSize)));
  static iREAL d[nsize] __attribute__ ((aligned(byteSize)));
  static iREAL e[nsize] __attribute__ ((aligned(byteSize)));
  static iREAL f[nsize] __attribute__ ((aligned(byteSize)));
  
  static iREAL p[10] __attribute__ ((aligned(byteSize)));
  static iREAL q[10] __attribute__ ((aligned(byteSize)));
  
  rndtriangles(n, 10, a, b, c, d, e, f);

  for (int i=0; i<n; i++)
  { 
    #pragma simd
    for (int j=i+1; j<n; j++)
    {      
        gjk (d, 4, e, 4, p, q);        
    }
  }
  return 0;
}
