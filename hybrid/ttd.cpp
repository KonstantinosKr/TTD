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
#define nsize 10000000

REAL penalty(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3], REAL  P[3], REAL  Q[3], REAL *error);

REAL pt(REAL TP1[3], REAL TP2[3], REAL TP3[3], REAL cPoint[3], REAL tq[3]);
REAL segseg(REAL p1[3], REAL p2[3], REAL p3[3], REAL p4[3], REAL P[3], REAL Q[3]);
REAL bf(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3], REAL  P[3], REAL  Q[3]);

REAL rndtriangles(unsigned n, REAL length, REAL a[][nsize], REAL b[][nsize], REAL c[][nsize], REAL d[][nsize], REAL e[][nsize], REAL f[][nsize])
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
  unsigned n = nsize;

  static REAL a[3][nsize] __attribute__ ((aligned(byteSize)));  
  static REAL b[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL c[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL d[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL e[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL f[3][nsize] __attribute__ ((aligned(byteSize)));
  static REAL p[3][nsize] __attribute__ ((aligned(byteSize))); 
  static REAL q[3][nsize] __attribute__ ((aligned(byteSize)));

  static REAL distance[nsize] __attribute__ ((aligned(byteSize))); 
   
  rndtriangles(n, 10, a, b, c, d, e, f);
  
  omp_set_num_threads(2);

  REAL A[3] __attribute__ ((aligned(byteSize)));
  REAL B[3] __attribute__ ((aligned(byteSize)));
  REAL C[3] __attribute__ ((aligned(byteSize)));    
  REAL D[3] __attribute__ ((aligned(byteSize)));    
  REAL E[3] __attribute__ ((aligned(byteSize)));    
  REAL F[3] __attribute__ ((aligned(byteSize)));    
  REAL Q[3] __attribute__ ((aligned(byteSize)));    
  REAL P[3] __attribute__ ((aligned(byteSize)));
      
  const unsigned batchSize = 4;
  //unsigned counter = 0;
  REAL error[batchSize] __attribute__ ((aligned(byteSize)));
  REAL batchError = 0.0;
  unsigned counter = 0;
  for (unsigned k=0; k<n/batchSize; k++) 
  { // works only for multiples of 16. Adopt      
        
        batchError = 0.0;
	#pragma forceinline recursive
//        #pragma simd reduction (+:batchError)
        #pragma omp parallel for simd reduction (+:batchError)
        for (unsigned l=0; l<batchSize; l++) 
        { 
            unsigned i=k*batchSize + l;
            
            A[0] = a[0][i];
            A[1] = a[1][i];
            A[2] = a[2][i];

            B[0] = b[0][i];
            B[1] = b[1][i];
            B[2] = b[2][i];
            
            C[0] = c[0][i];
            C[1] = c[1][i];
            C[2] = c[2][i];
            
            D[0] = d[0][i];
            D[1] = d[1][i];
            D[2] = d[2][i];
            
            E[0] = e[0][i];
            E[1] = e[1][i];
            E[2] = e[2][i];
            
            F[0] = f[0][i];
            F[1] = f[1][i];
            F[2] = f[2][i];
            
            distance[i] = penalty(A, B, C, D, E, F, P, Q, &error[l]);  // only four Newton steps
            
            p[0][i] = P[0];
            p[1][i] = P[1];
            p[2][i] = P[2];
            
            q[0][i] = Q[0];
            q[1][i] = Q[1];
            q[2][i] = Q[2];
        } // end of batch
        batchError = 0.0;
        for (unsigned l=0; l<batchSize; l++) 
        {
            batchError += error[l];
        } // end of batch
    if(batchError/batchSize>1E-8)
    {
        counter++;
        #pragma forceinline recursive
	//#pragma simd 
	#pragma omp parallel for simd 
	for (unsigned l=0; l<batchSize; l++)
        {
            unsigned i = k*batchSize + l;
            distance[i] = bf(A,B,C, D,E,F, P,Q);

            p[0][i] = P[0];
            p[1][i] = P[1];
            p[2][i] = P[2];
            
            q[0][i] = Q[0];
            q[1][i] = Q[1];
            q[2][i] = Q[2];
        }
    }
  }


  printf("counter: %i",counter);
  return 0;
}
    
    REAL penalty(REAL A[3], REAL B[3], REAL  C[3], REAL  D[3], REAL  E[3], REAL  F[3], REAL  P[3], REAL  Q[3], REAL error[])
    {
        REAL BA[3], CA[3], ED[3], FD[3], hessian[16], x[4];
        
        SUB(B,A, BA);
        SUB(C,A, CA);
        SUB(E,D, ED);
        SUB(F,D, FD);
        
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
        
        REAL eps = 1E-16;
        REAL delta = (hessian[0]+hessian[5]+hessian[10]+hessian[15]) * eps;
        REAL lambda = sqrt(0.0125*(hessian[0]+hessian[5]+hessian[10]+hessian[15]));
        REAL r = lambda*pow(10,log10(lambda)+8);
        
#if REAL==double
        REAL tol = DBL_EPSILON;
#else
        REAL tol = FLT_EPSILON;
#endif
        //initial guess
        x[0] = 0.33;
        x[1] = 0.33;
        x[2] = 0.33;
        x[3] = 0.33;
        
        //Declare loop variables;
        REAL dx[4], a[16], SUBXY[3], b[4], mx[6], dh[8], tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
        
        //Newton loop
        for(int i=0;i<7;i++)
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
            
            SUBXY[0] = (A[0]+(BA[0] * x[0])+(CA[0] * x[1])) - (D[0]+(ED[0] * x[2])+(FD[0] * x[3]));
            SUBXY[1] = (A[1]+(BA[1] * x[0])+(CA[1] * x[1])) - (D[1]+(ED[1] * x[2])+(FD[1] * x[3]));
            SUBXY[2] = (A[2]+(BA[2] * x[0])+(CA[2] * x[1])) - (D[2]+(ED[2] * x[2])+(FD[2] * x[3]));
            
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

        *error = sqrt(DOT4(dx,dx)/DOT4(x,x)); 
        
        P[0] = A[0]+(BA[0] * x[0])+(CA[0] * x[1]);
        P[1] = A[1]+(BA[1] * x[0])+(CA[1] * x[1]);
        P[2] = A[2]+(BA[2] * x[0])+(CA[2] * x[1]);
        
        Q[0] = D[0]+(ED[0] * x[2])+(FD[0] * x[3]);
        Q[1] = D[1]+(ED[1] * x[2])+(FD[1] * x[3]);
        Q[2] = D[2]+(ED[2] * x[2])+(FD[2] * x[3]);
        SUBXY[0] = (P[0]) - (Q[0]);
        SUBXY[1] = (P[1]) - (Q[1]);
        SUBXY[2] = (P[2]) - (Q[2]);
        return sqrt(DOT(SUBXY,SUBXY));
    }

int segt(REAL p1[3], REAL p2[3], REAL A[3], REAL B[3], REAL C[3], REAL P[3]){
  //77 flops
    REAL u[3], v[3], dir[3], w0[3], n[3], a, b, r, uu, uv, vv, w[3], wu, wv, s, D, t;
    
    u[0] = B[0] - A[0];
    u[1] = B[1] - A[1];
    u[2] = B[2] - A[2];
    
    v[0] = C[0] - A[0];
    v[1] = C[1] - A[1];
    v[2] = C[2] - A[2];
    
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];//15 flops
    
    if(n[0]==0 && n[1]==0 && n[2]==0)
        return -1;
    
    dir[0] = p2[0] - p1[0];
    dir[1] = p2[1] - p1[1];
    dir[2] = p2[2] - p1[2];
    
    w0[0] = p1[0] - A[0];
    w0[1] = p1[1] - A[1];
    w0[2] = p1[2] - A[2];//6 flops
    
    a = -DOT(n, w0);
    b = DOT(n,dir);//10 flops
    if (fabs(b) < 1E-99) {
        if (a==0) //segment in triangle plane
            return 2;
    }
    
    r = a / b;
    if (r < 0.0 || r > 1.0 || r == NAN)
        return 0;
    
    P[0] = p1[0] + r * dir[0];
    P[1] = p1[1] + r * dir[1];
    P[2] = p1[2] + r * dir[2];//7 flops
    
    uu = DOT(u,u);
    uv = DOT(u,v);
    vv = DOT(v,v);
    w[0] = P[0] - A[0];
    w[1] = P[1] - A[1];
    w[2] = P[2] - A[2];
    
    wu = DOT(w,u);
    wv = DOT(w,v);
    D = uv*uv - uu*vv;
    
    s = (uv*wv - vv*wu)/D;
    if (s<0.0 || s>1.0)
        return 0;
    
    t = (uv*wu - uu*wv)/D;
    if (t<0.0 || (s+t)>1.0)
        return 0;//39 flops
    
    return 1; //seg in triangle
}

inline REAL segseg(REAL p1[3], REAL p2[3], REAL p3[3], REAL p4[3], REAL P[3], REAL Q[3])
{
    //69 flops
    REAL u[3];
    u[0] = p1[0] - p2[0];
    u[1] = p1[1] - p2[1];
    u[2] = p1[2] - p2[2];
    
    REAL v[3];
    v[0] = p3[0] - p4[0];
    v[1] = p3[1] - p4[1];
    v[2] = p3[2] - p4[2];
    
    REAL w[3];
    w[0] = p2[0] - p4[0];
    w[1] = p2[1] - p4[1];
    w[2] = p2[2] - p4[2];
    
    REAL a = DOT(u,u);
    REAL b = DOT(u,v);
    REAL c = DOT(v,v);
    REAL d = DOT(u,w);
    REAL e = DOT(v,w);
    REAL D = a*c - b*b;
    REAL sD = D;
    REAL tD = D;
    REAL sN=0, tN=0;
    
    REAL SMALL_NUM = 0.00000001;
    
    // compute the line parameters of the two closest points
    if (D < SMALL_NUM){//  //% the lines are almost parallel
        sN = 0.0;       //% force using point P0 on segment S1
        sD = 1.0;       //% to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else {               //% get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0){   //% sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }else if (sN > sD){//% sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
    if (tN < 0.0){     //% tc < 0 => the t=0 edge is visible
        tN = 0.0;
        //% recompute sc for this edge
        if (-d < 0.0){
            sN = 0.0;
        }else if (-d > a){
            sN = sD;
        } else {
            sN = -d;
            sD = a;
        }
        
    } else if (tN > tD){       //% tc > 1 => the t=1 edge is visible
        tN = tD;
        //% recompute sc for this edge
        if ((-d + b) < 0.0){
            sN = 0;
        }else if ((-d + b) > a){
            sN = sD;
        } else {
            sN = (-d + b);
            sD = a;
        }
    }
    REAL sc;
    //% finally do the division to get sc and tc
    if(fabs(sN) < SMALL_NUM){
        sc = 0.0;
    }else{
        sc = sN / sD;
    }
    REAL tc;
    if(fabs(tN) < SMALL_NUM){
        tc = 0.0;
    }else{
        tc = tN / tD;
    }
    
    P[0] = p2[0] + u[0] * sc;
    P[1] = p2[1] + u[1] * sc;
    P[2] = p2[2] + u[2] * sc;
    
    Q[0] = p4[0] + v[0] * tc;
    Q[1] = p4[0] + v[1] * tc;
    Q[2] = p4[0] + v[2] * tc;
    
    REAL dP[3];
    dP[0] = w[0] + (sc*u[0]) - (tc*v[0]);
    dP[1] = w[1] + (sc*u[1]) - (tc*v[1]);
    dP[2] = w[2] + (sc*u[2]) - (tc*v[2]);
    
    return sqrt(DOT(dP,dP));
}

REAL pt(REAL TP1[3], REAL TP2[3], REAL TP3[3], REAL cPoint[3], REAL tq[3])
{
    //191 flops
    REAL E0[3];
    E0[0] = TP2[0] - TP1[0];
    E0[1] = TP2[1] - TP1[1];
    E0[2] = TP2[2] - TP1[2];
    
    REAL E1[3];
    E1[0] = TP3[0] - TP1[0];
    E1[1] = TP3[1] - TP1[1];
    E1[2] = TP3[2] - TP1[2];
    
    REAL D[3];
    D[0] = TP1[0] - cPoint[0];
    D[1] = TP1[1] - cPoint[1];
    D[2] = TP1[2] - cPoint[2];
    
    REAL a = DOT(E0,E0);
    REAL b = DOT(E0,E1);
    REAL c = DOT(E1,E1);
    REAL d = DOT(E0,D);
    REAL e = DOT(E1,D);
    REAL f = DOT(D,D);
    
    REAL det = a*c - b*b; //% do we have to use abs here?
    REAL s   = b*e - c*d;
    REAL t   = b*d - a*e;
    
    REAL sqrDistance=0;
    
    if ((s+t) <= det){
        if (s < 0){
            if (t < 0){
                //region4
                if (d < 0){
                    t = 0;
                    if (-d >= a){
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    }else {
                        s = -d/a;
                        sqrDistance = d*s + f;
                    }
                }else {
                    s = 0;
                    if (e >= 0){
                        t = 0;
                        sqrDistance = f;
                    }else{
                        if (-e >= c){
                            t = 1;
                            sqrDistance = c + 2*e + f;
                        } else {
                            t = -e/c;
                            sqrDistance = e*t + f;
                        }
                    }
                } //end of region 4
            }else {
                // region 3
                s = 0;
                if (e >= 0){
                    t = 0;
                    sqrDistance = f;
                }else {
                    if (-e >= c){
                        t = 1;
                        sqrDistance = c + 2*e +f;
                    }else {
                        t = -e/c;
                        sqrDistance = e*t + f;
                    }
                }
            } //end of region 3
        }else {
            if (t < 0){
                // region 5
                t = 0;
                if (d >= 0){
                    s = 0;
                    sqrDistance = f;
                }else {
                    if (-d >= a){
                        s = 1;
                        sqrDistance = a + 2*d + f;// GF 20101013 fixed typo d*s ->2*d
                    }else {
                        s = -d/a;
                        sqrDistance = d*s + f;
                    }
                }
            }else {
                // region 0
                REAL invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
            }
        }
    }else {
        if (s < 0){
            // region 2
            REAL tmp0 = b + d;
            REAL tmp1 = c + e;
            if (tmp1 > tmp0){ // minimum on edge s+t=1
                REAL numer = tmp1 - tmp0;
                REAL denom = a - 2*b + c;
                if (numer >= denom){
                    s = 1;
                    t = 0;
                    sqrDistance = a + 2*d + f; // GF 20101014 fixed typo 2*b -> 2*d
                }else {
                    s = numer/denom;
                    t = 1-s;
                    sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                }
            }else {         // minimum on edge s=0
                s = 0;
                if (tmp1 <= 0) {
                    t = 1;
                    sqrDistance = c + 2*e + f;
                }else {
                    if (e >= 0){
                        t = 0;
                        sqrDistance = f;
                    }else {
                        t = -e/c;
                        sqrDistance = e*t + f;
                    }
                }
            } //end of region 2
        }else {
            if (t < 0) {
                //region6
                REAL tmp0 = b + e;
                REAL tmp1 = a + d;
                if (tmp1 > tmp0){
                    REAL numer = tmp1 - tmp0;
                    REAL denom = a-2*b+c;
                    if (numer >= denom){
                        t = 1;
                        s = 0;
                        sqrDistance = c + 2*e + f;
                    }else {
                        t = numer/denom;
                        s = 1 - t;
                        sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                }else {
                    t = 0;
                    if (tmp1 <= 0){
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    }else {
                        if (d >= 0) {
                            s = 0;
                            sqrDistance = f;
                        }else {
                            s = -d/a;
                            sqrDistance = d*s + f;
                        }
                    }
                }
                //end of region 6
            }else {
                // region 1
                REAL numer = c + e - b - d;
                if (numer <= 0){
                    s = 0;
                    t = 1;
                    sqrDistance = c + 2*e + f;
                }else {
                    REAL denom = a - 2*b + c;
                    if (numer >= denom){
                        s = 1;
                        t = 0;
                        sqrDistance = a + 2*d + f;
                    }else {
                        s = numer/denom;
                        t = 1-s;
                        sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                } //end of region 1
            }
        }
    }
    
    // account for numerical round-off error
    if (sqrDistance < 0){
        sqrDistance = 0;
    }
    
    tq[0] = TP1[0] + (E1[0] * t) + (E0[0] * s);
    tq[1] = TP1[1] + (E1[1] * t) + (E0[1] * s);
    tq[2] = TP1[2] + (E1[2] * t) + (E0[2] * s);
    
    return sqrt(sqrDistance);
}

REAL bf(REAL A[3], REAL B[3], REAL C[3], REAL D[3], REAL E[3], REAL F[3], REAL P[3], REAL Q[3]){   
    //seg to triangle intersection
    int intersection = 0;
    
    #if 1
    if (segt(D, E, A, B, C, P) == 1)
    {
        intersection = 1;
    } else if (segt(E, F, A, B, C, P) == 1)
    {
        intersection = 1;
    } else if (segt(F, D, A, B, C, P) == 1)
    {
        intersection = 1;
    } else if (segt(A, B, D, E, F, P) == 1)
    {
        intersection = 1;
    } else if (segt(B, C, D, E, F, P) == 1)
    {
        intersection = 1;
    } else if (segt(C, A, D, E, F, P) == 1)
    {
        intersection = 1;
    }
    #endif 

    //%--------Point to Triangle Minimum Distance
    REAL tp[3], tq[3], tmpP[3], ptmin=1E+30, tmp, ssP[3], ssQ[3], tmpQ[3], ssmin=1E+30;

    if (intersection != 1)
    {
      //%Vertex of T1 to Triangle T2
      ptmin= pt(D,E,F, A, tmpP);
      tp[0] = A[0];
      tp[1] = A[1];
      tp[2] = A[2];
    
      tq[0] = tmpP[0];
      tq[1] = tmpP[1];
      tq[2] = tmpP[2];
    
      tmp= pt(D,E,F, B, tmpP);
      if(ptmin>tmp){
        ptmin = tmp;
        tp[0] = B[0];
        tp[1] = B[1];
        tp[2] = B[2];
        
        tq[0] = tmpP[0];
        tq[1] = tmpP[1];
        tq[2] = tmpP[2];
      }
    
      tmp= pt(D,E,F, C, tmpP);
      if(ptmin>tmp){
        ptmin = tmp;
        tp[0] = C[0];
        tp[1] = C[1];
        tp[2] = C[2];
        
        tq[0] = tmpP[0];
        tq[1] = tmpP[1];
        tq[2] = tmpP[2];
      }
    
      //%Vertex of T2 to Triangle T1
      tmp= pt(A,B,C, D, tmpP);
      if(ptmin>tmp){
        ptmin = tmp;
        tq[0] = D[0];
        tq[1] = D[1];
        tq[2] = D[2];
        
        tp[0] = tmpP[0];
        tp[1] = tmpP[1];
        tp[2] = tmpP[2];
      }
    
      tmp= pt(A,B,C, E, tmpP);
      if(ptmin>tmp){
        ptmin = tmp;
        tq[0] = E[0];
        tq[1] = E[1];
        tq[2] = E[2];
        
        tp[0] = tmpP[0];
        tp[1] = tmpP[1];
        tp[2] = tmpP[2];
      }
    
      tmp= pt(A,B,C, F, tmpP);
      if(ptmin>tmp){
        ptmin = tmp;
        tq[0] = F[0];
        tq[1] = F[1];
        tq[2] = F[2];
        
        tp[0] = tmpP[0];
        tp[1] = tmpP[1];
        tp[2] = tmpP[2];
      }

      //-------Seg-Seg Distance--------
      ssmin= segseg(A,B, D,E, tmpP, tmpQ);
      ssP[0] = tmpP[0];
      ssP[1] = tmpP[1];
      ssP[2] = tmpP[2];
    
      ssQ[0] = tmpQ[0];
      ssQ[1] = tmpQ[1];
      ssQ[2] = tmpQ[2];
      tmp= segseg(A,B,E,F, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(A,B, F,D, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(B,C, D,E, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(B,C, E,F, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(B,C, F,D, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(C,A, D,E, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(C,A, E,F, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
      tmp= segseg(C,A, F,D, tmpP, tmpQ);
      if(ssmin>tmp){
        ssmin = tmp;
        
        ssP[0] = tmpP[0];
        ssP[1] = tmpP[1];
        ssP[2] = tmpP[2];
        
        ssQ[0] = tmpQ[0];
        ssQ[1] = tmpQ[1];
        ssQ[2] = tmpQ[2];
      }
    }
    
    //%---Get min distance from ss and point to triangle sets
    REAL min=-1;
    if(intersection == 1){ 
        min = 0;

        Q[0] = P[0];
        Q[1] = P[1];
        Q[2] = P[2];
    } else if(ssmin < ptmin){
        min = ssmin;
        P[0] = ssP[0];
        P[1] = ssP[1];
        P[2] = ssP[2];
        
        Q[0] = ssQ[0];
        Q[1] = ssQ[1];
        Q[2] = ssQ[2];
    } else if(ssmin > ptmin){
        min = ptmin;
        P[0] = tp[0];
        P[1] = tp[1];
        P[2] = tp[2];
        
        Q[0] = tq[0];
        Q[1] = tq[1];
        Q[2] = tq[2];
    } else if(ssmin == ptmin){
        min = ptmin;
        P[0] = ssP[0];
        P[1] = ssP[1];
        P[2] = ssP[2];
        
        Q[0] = tq[0];
        Q[1] = tq[1];
        Q[2] = tq[2];
    }
    return min;
}

