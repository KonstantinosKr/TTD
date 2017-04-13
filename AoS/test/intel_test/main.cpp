#include <stdio.h>
#include "timing.h"
#define nsize 10000000
#define byteSize 32

inline void test (int n, REAL *restrict a[nsize], REAL *restrict b[nsize], REAL *restrict p[nsize], REAL *restrict q[nsize], REAL distance[]); 

int main (int argc, char *argv[])
{
  const static int n= nsize;

  static REAL distance[nsize] __attribute__ ((aligned(byteSize)));

  REAL *a[3] __attribute__ ((aligned(byteSize)));
  REAL *b[3] __attribute__ ((aligned(byteSize)));
  REAL *p[3] __attribute__ ((aligned(byteSize)));
  REAL *q[3] __attribute__ ((aligned(byteSize)));
  
  a[0]= (REAL*) _mm_malloc (sizeof(REAL)*nsize, byteSize);
  b[0]= (REAL*) _mm_malloc (sizeof(REAL)*nsize, byteSize);
  p[0]= (REAL*) _mm_malloc (sizeof(REAL)*nsize, byteSize);
  q[0]= (REAL*) _mm_malloc (sizeof(REAL)*nsize, byteSize);
  
  for(int i=0;i<n;i=i+2)
  {
    a[0][i] = 99;

    b[0][i] = 99;
    
    p[0][i] = 99;
    q[0][i] = 99;
    distance[i] = 1;
  }

  test (n, a, b, p, q, distance); 

  return 0;
}

//__attribute__((vector))
inline void test (int n, REAL *restrict a[nsize], REAL *restrict b[nsize], REAL *restrict p[nsize], REAL *restrict q[nsize], REAL distance[]) 
{
  REAL aa __attribute__((aligned(byteSize)));
  REAL bb __attribute__((aligned(byteSize)));
  REAL out __attribute__((aligned(byteSize)));

  __assume_aligned(a, byteSize);
  __assume_aligned(b, byteSize);
  __assume_aligned(p, byteSize);
  __assume_aligned(q, byteSize);
  __assume_aligned(distance, byteSize);
 
  //#pragma forceinline recursive
  #pragma simd
  for (int i = 0; i < n; i ++)
  {
    aa = a[0][i];
   
    bb = b[0][i];
    
    out = distance[i] * (aa+aa+aa)*(bb+bb+bb);

    p[0][i] = out;

    q[0][i] = out;

    distance[i] = out;
  }
}
