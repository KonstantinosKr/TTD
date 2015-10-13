#include <stdio.h>
#include "timing.h"
#include "test_ispc.h"
#define nsize 100000
using namespace ispc;

void tuneVectors(int n, REAL *a[3], REAL *b[3], REAL *c[3], REAL *d[3], REAL *e[3], REAL *f[3], REAL *p[3], REAL *q[3], REAL distance[]);

int main (int argc, char *argv[])
{
  int n= nsize;

  static REAL *a[3];
  static REAL *b[3];
  static REAL *c[3];
  static REAL *d[3];
  static REAL *e[3];
  static REAL *f[3];
  static REAL *p[3];
  static REAL *q[3];
  REAL distance[n];
  
  for(int i=0;i<3;i++)
  {
      a[i]= ispc_allocate (n);
      b[i]= ispc_allocate (n);
      c[i]= ispc_allocate (n);
      d[i]= ispc_allocate (n);
      e[i]= ispc_allocate (n);
      f[i]= ispc_allocate (n);
      p[i]= ispc_allocate (n);
      q[i]= ispc_allocate (n); 
  }
  
  for(int i=0;i<n;i++)
  {
    a[0][i] = 99;
    a[1][i] = 99;
    a[2][i] = 99; 

    b[0][i] = 99;
    b[1][i] = 99;
    b[2][i] = 99;

    c[0][i] = 99;
    c[1][i] = 99;
    c[2][i] = 99;

    d[0][i] = 99;
    d[1][i] = 99;
    d[2][i] = 99;

    e[0][i] = 99;
    e[1][i] = 99; 
    e[2][i] = 99;

    f[0][i] = 99;
    f[1][i] = 99;
    f[2][i] = 99; 
  }

  tuneVectors(n,a,b,c,d,e,f,p,q, distance);

  return 0;
}

void serial_test (int n, int m, REAL *a[3], REAL *b[3], REAL *c[3], REAL *d[3], REAL *e[3], REAL *f[3], REAL *p[3], REAL *q[3], REAL distance[])
{
  for (int i = 0; i < n; i ++)
  {
    REAL aa[3], bb[3], cc[3], dd[3], ee[3], ff[3];

    aa[0] = a[0][i];
    aa[1] = a[1][i];
    aa[2] = a[2][i];

    bb[0] = b[0][i];
    bb[1] = b[1][i];
    bb[2] = b[2][i];

    cc[0] = c[0][i];
    cc[1] = c[1][i];
    cc[2] = c[2][i];

    dd[0] = d[0][i];
    dd[1] = d[1][i];
    dd[2] = d[2][i];

    ee[0] = e[0][i];
    ee[1] = e[1][i];
    ee[2] = e[2][i];

    ff[0] = f[0][i];
    ff[1] = f[1][i];
    ff[2] = f[2][i];

    REAL out;
    for (int j = 0; j < m; j ++)
    {
      out = (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2])+
            (aa[0]+aa[1]+aa[2])*(bb[0]+bb[1]+bb[2])*(cc[0]+cc[1]+cc[2])*(dd[0]+dd[1]+dd[2])*(ee[0]+ee[1]+ee[2])*(ff[0]+ff[1]+ff[2]);
    }

    p[0][i] = out;
    p[1][i] = out;
    p[2][i] = out;

    q[0][i] = out;
    q[1][i] = out;
    q[2][i] = out;

    distance[i] = out;
  }
}

void tuneVectors(int n, REAL a[3], REAL *b[3], REAL *c[3], REAL *d[3], REAL *e[3], REAL f[3], REAL *p[3], REAL *q[3], REAL distance[])
{
  int m = 5, l=100;

  REAL tSerial = 0.0;
  reset_and_start_timer();
  for(int j=0;j<l;j++)
  {
    serial_test (n,m,a,b,c,d,e,f,p,q,distance);
  }
  tSerial += get_elapsed_mcycles();

  REAL tISPC = 0.0;
  reset_and_start_timer();
  for(int j=0;j<l;j++)
  {
    ispc_test(n,m,a,b,c,d,e,f,p,q,distance);
  }
  tISPC += get_elapsed_mcycles();

  printf("Test, m = %d, speedup %g\n", m, tSerial/tISPC);
}
