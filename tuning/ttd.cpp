#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "algo.h"
#include "math.h"
#include <ctime>
#include <string.h>

extern REAL serial_pen_test(unsigned n, REAL r, REAL eps, REAL* a[3], REAL* b[3], REAL*  c[3], REAL*  d[3], REAL*  e[3], REAL*  f[3], REAL*  p[3], REAL*  q[3], REAL distance[], REAL fails[], REAL *failrate, unsigned totit[]);
extern REAL serial_pen_test2(unsigned n, REAL eps, REAL setDelta, REAL* a[3], REAL* b[3], REAL*  c[3], REAL*  d[3], REAL*  e[3], REAL*  f[3], REAL*  p[3], REAL*  q[3], REAL distance[], REAL fails[], REAL *failrate, unsigned totit[]);
extern void serial_bf(unsigned n, REAL* a[3], REAL* b[3], REAL* c[3], REAL* d[3], REAL* e[3], REAL* f[3], REAL* p[3], REAL* q[3], REAL distance[]); 

void tPen(unsigned n, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3], REAL * p[3], REAL * q[3]);
void rndtriangles(unsigned n, REAL sidelgth, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3]);

void tune(int lgth_lo_pow, int lgth_hi_pow, int r_lo_pow, int r_hi_pow, int eps_lo_pow, int eps_hi_pow,  unsigned n, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3], REAL * p[3], REAL * q[3]);

extern void loadTriangles(unsigned n, REAL* a[3], REAL* b[3], REAL*  c[3], REAL*  d[3], REAL*  e[3], REAL*  f[3], REAL length[]);

void rndtriangles(unsigned n, REAL length, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3])
{
  int padding = 2*length;
  
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
}

REAL getTriangleLength(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3])
{
    REAL BA[3], CA[3], ED[3], FD[3], hessian[4];

    SUB(B,A, BA);
    SUB(C,A, CA);
    SUB(E,D, ED);
    SUB(F,D, FD);//12 FLOPS

    hessian[0] = 2.*DOT(BA,BA);
    hessian[1] = 2.*DOT(CA,CA);
    hessian[2] = 2.*DOT(ED,ED);
    hessian[3] = 2.*DOT(FD,FD);//10 FLOPS + 5 FLOPS

    return sqrt(0.0125*(hessian[0]+hessian[1]+hessian[2]+hessian[3])); 
}

void loadTriangles(unsigned n, REAL* a[3], REAL* b[3], REAL*  c[3], REAL*  d[3], REAL*  e[3], REAL*  f[3], REAL length[])
{
    REAL A[3],B[3], C[3], D[3], E[3], F[3];
        
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

        length[i] = getTriangleLength(A, B, C, D, E, F);
    }
}
        
int main (int argc, char *argv[])
{
  static unsigned n= 100000;

  REAL *a[3];
  REAL *b[3];
  REAL *c[3];
  REAL *d[3];
  REAL *e[3];
  REAL *f[3];
  REAL *p[3];
  REAL *q[3];
  REAL distance[n];
  REAL length[n];

  for(unsigned i=0;i<3;i++)
  {
    a[i]= new REAL [n];
    b[i]= new REAL [n];
    c[i]= new REAL [n];
    d[i]= new REAL [n];
    e[i]= new REAL [n];
    f[i]= new REAL [n];
    p[i]= new REAL [n];
    q[i]= new REAL [n]; 
  }
 
  rndtriangles(n, 10, a, b, c, d, e, f);

  tPen(n, a, b, c, d, e, f, p, q);
  tune(1, 2, -6, 10, -15, 0,  n, a, b, c, d, e, f, p, q); printf("finished penalty tuning!\n");
  return 0;
}

void tune(int lgth_lo_pow, int lgth_hi_pow, int r_lo_pow, int r_hi_pow, int eps_lo_pow, int eps_hi_pow, unsigned n, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3], REAL * p[3], REAL * q[3])
{
    int u, i, j, ii, jj, kk, zz, uu, xx, t, tt;
    REAL ll;
    int nr = r_hi_pow-r_lo_pow+1,
    neps = eps_hi_pow-eps_lo_pow+1,
    nlgth = lgth_hi_pow-lgth_lo_pow+1;
    
    REAL *lgth = new REAL[nr*neps*nlgth];
    REAL *r = new REAL[nr*neps*nlgth];
    REAL *eps = new REAL[nr*neps*nlgth];
    REAL *failrate = new REAL[nr*neps*nlgth];
    unsigned *totit = new unsigned[n];
    unsigned *stotit = new unsigned[nr*neps*nlgth];
    REAL *fails = new REAL[n];
    REAL *fails_opt = new REAL[n];
    REAL *distance = new REAL[n];
    REAL *distance2 = new REAL[n];
    REAL *p1[3], *p2[3], *q1[3], *q2[3];
    REAL *error = new REAL[nr*neps*nlgth];

    for(unsigned i=0;i<3;i++)
    {
        p1[i] = new REAL[n];
        p2[i] = new REAL[n];
        q1[i] = new REAL[n];
        q2[i] = new REAL[n];
        totit[i] = 0;
    }
    
    for (u = lgth_lo_pow; u <= lgth_hi_pow; u++)
    {
        uu = u - lgth_lo_pow;
        ll = pow(10.,u);
        rndtriangles(n, ll, a, b, c, d, e, f);
        
        REAL r_opt = 0, eps_opt = 0, maxit_opt = 1E99, error_opt=100, failrate_opt=100, fail_opt = -99;
        int swiched = 0; 
        for (i = r_lo_pow; i <= r_hi_pow; i ++)
        {
            ii = i-r_lo_pow;
            zz = uu*nr+ii;
            
            for (j = eps_lo_pow; j <= eps_hi_pow; j ++)
            {
                jj = j-eps_lo_pow;
                kk = zz*neps+jj;
                
                lgth[kk] = ll;
                r[kk] = pow(10.,i);
                eps[kk] = pow(10.,j);

                error[kk] = serial_pen_test (n, r[kk], eps[kk], a, b, c, d, e, f, p, q, distance, fails, &failrate[kk], totit);
                    
                stotit[kk] = 0;
                for(unsigned z=0;z<n;z++) stotit[kk] = stotit[kk] + totit[z];
                    
                if(error[kk] < 1E-2 && swiched != 1)
                { 
                    swiched = 1;
                    r_opt = r[kk];
                    eps_opt = eps[kk];
                    failrate_opt = failrate[kk];
                    error_opt = error[kk];
                    maxit_opt = stotit[kk];
                }
            }
        }
        printf ("len = %g, r_opt = %g, eps_opt = %g, maxit_opt = %g, error = %g, failure_rate = %g\n", ll, r_opt, eps_opt, maxit_opt, error_opt, failrate_opt);
    }
    
    FILE *fp = fopen("ptune.data", "w+");
    fprintf(fp,"ID, lgth, r, eps, failrate%%, totalit, error, neps, nr, nlgth\n");
    for(i = 0; i < nr*neps*nlgth; i++)
    {
        fprintf(fp,"%i, %g, %g, %g, %g, %i, %g, %i, %i, %i\n", i, lgth[i], r[i], eps[i], failrate[i], stotit[i], error[i], neps, nr, nlgth); 
    }
    fclose(fp); 
    
    delete r;
    delete eps;
    delete failrate;
    delete fails;
    delete totit;
}

void tPen(unsigned n, REAL * a[3], REAL * b[3], REAL * c[3], REAL * d[3], REAL * e[3], REAL * f[3], REAL * p[3], REAL * q[3])
{
  REAL *pedserial = new REAL[n];
  REAL *bfdserial = new REAL[n];
  REAL *fails = new REAL[n];
  unsigned *totit = new unsigned[n];  
  REAL *p1[3], *p2[3], *q1[3], *q2[3];
  for(unsigned i=0;i<3;i++)
  {
    p1[i] = new REAL[n];
    p2[i] = new REAL[n];
    q1[i] = new REAL[n];
    q2[i] = new REAL[n];
  }
  
  //for(int o=1;o<4;o++)
  {
      REAL failrate = -1;

      REAL lambda = 10;
      rndtriangles(n, lambda, a,b,c,d,e,f);  

      REAL eps = 1E-15;
      REAL r = lambda*pow(10,log10(lambda)+8);
      REAL error = serial_pen_test (n, r, eps, a, b, c, d, e, f, p1, q1, pedserial, fails, &failrate, totit);
      
      unsigned stotit=0;
      for(unsigned i=0;i<n;i++) stotit = stotit + totit[i];
      printf("Penalty Serial. Failrate:%f. Error:%f. Total it:%i\n", failrate, error, stotit);
      
      serial_bf (n, a, b, c, d, e, f, p2, q2, bfdserial);

      FILE *fp1 = fopen("output.data", "w+");
      FILE *fp2 = fopen("input.data", "w+");

      fprintf(fp1,"ID, p1[0], p2[0], p1[1], p2[1], p1[2], p2[2], q1[0], q2[0], q1[1], q2[1], q1[2], q2[2], penSerial, bfSerial, failed, totalit\n");
      fprintf(fp2,"ID, a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], e[0], e[1], e[2], f[0], f[1], f[2], penSerial, bfSerial, failed, totalit\n");
      
      unsigned counter = 0;
      for(unsigned i = 0; i < n; i++)
      {
          {
                fprintf(fp1,"%i, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %i\n", \
                counter, p1[0][i], p2[0][i], p1[1][i], p2[1][i], p1[2][i], p2[2][i], q1[0][i], q2[0][i], q1[1][i], q2[1][i], q1[2][i], q2[2][i], pedserial[i], bfdserial[i], fails[i], totit[i]);

                fprintf(fp2,"%i, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.2f, %.2f, %.2f, %i\n", \
                counter, a[0][i], a[1][i], a[2][i], b[0][i], b[1][i], b[2][i], c[0][i], c[1][i], c[2][i], d[0][i], d[1][i], d[2][i], e[0][i], e[1][i], e[2][i], f[0][i], f[1][i], f[2][i], pedserial[i], bfdserial[i], fails[i], totit[i]); 
            
                counter++;
          }
      }
      fclose(fp1);
      fclose(fp2); 
  }
}
