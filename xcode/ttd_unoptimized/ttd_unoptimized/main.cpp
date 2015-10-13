//
//  main.cpp
//  main
//
//  Created by Konstantinos Krestenitis on 9/5/14.
//  Copyright (c) 2014 Durham University. All rights reserved.
//

#include <algorithm>
#include "math.h"
#include "algo.h"
#include <stdlib.h>

using namespace std;

double penalty(double r, double eps, double tol, double A[3], double B[3], double  C[3], double  D[3], double  E[3], double  F[3], double  P[3], double  Q[3]);
void ini(int n, double r, double eps, double tol, double* a[3], double* b[3], double*  c[3], double*  d[3], double*  e[3], double*  f[3], double*  p[3], double*  q[3], double distance[], double fails[]);
double pt(double* TP1, double* TP2, double* TP3, double* cPoint, double tPoint[3]);
double segseg(double* p1,double* p2,double* p3,double* p4,double P1[3], double P2[3]);
int segt(double p1[3], double p2[3], double A[3], double B[3], double C[3], double P[3]);
double bf(double A[3], double B[3], double C[3], double D[3], double E[3], double F[3], double P[3], double Q[3]);
double failureRate(int n, double distance[], double distcmp[], double P1[], double P2[], double Q1[], double Q2[], double epsilon, double fails[]);
int main()
{
    int n = 1;
    
    double *a[3];
    double *b[3];
    double *c[3];
    double *d[3];
    double *e[3];
    double *f[3];
    double *p[3];
    double *q[3];
    double distance[n];
    double fails= 0;
    
    
    for(unsigned int i=0;i<3;i++)
    {
        a[i] = new double[n];
        b[i] = new double[n];
        c[i] = new double[n];
        d[i] = new double[n];
        e[i] = new double[n];
        f[i] = new double[n];
        p[i] = new double[n];
        q[i] = new double[n];
    }
    
    for(int i=0;i<n;i++)
    {

        a[0][i]=83.7429;
        a[1][i]=67.2279;
        a[2][i]=98.6072;
        
        b[0][i]=91.7404;
        b[1][i]=81.3977;
        b[2][i]=50.5941;
        
        c[0][i]=34.5489;
        c[1][i]=63.1802;
        c[2][i]=69.0500;
        
        d[0][i]=122.7448;
        d[1][i]=172.2445;
        d[2][i]=113.2040;
        
        e[0][i]=119.6180;
        e[1][i]=120.0869;
        e[2][i]=100.5043;
        
        f[0][i]=175.4586;
        f[1][i]=132.4482;
        f[2][i]=156.3177;
        
        
        
        a[0][i]=4920.3139;
        a[1][i]=5716.4345;
        a[2][i]=6115.2058;
        
        b[0][i]=8264.0094;
        b[1][i]=3206.4150;
        b[2][i]=217.6251;
        
        c[0][i]=7624.4028;
        c[1][i]=3337.3811;
        c[2][i]=1363.9438;
        
        d[0][i]=13803.9306;
        d[1][i]=15856.1094;
        d[2][i]=16824.0762;
        
        e[0][i]=15441.9754;
        e[1][i]=16475.3678;
        e[2][i]=14700.1829;
        
        f[0][i]=19167.5150;
        f[1][i]=11618.1160;
        f[2][i]=18869.7270;
        
        
        a[0][i]=324.0070;
        a[1][i]=586.0559;
        a[2][i]=840.9284;
        
        b[0][i]=483.7652;
        b[1][i]=641.2078;
        b[2][i]=779.1985;
        
        c[0][i]=989.8194;
        c[1][i]=894.4916;
        c[2][i]=719.9269;
        
        d[0][i]=1812.6839;
        d[1][i]=1972.3481;
        d[2][i]=1448.9996;
        
        e[0][i]=1530.7792;
        e[1][i]=1999.3999;
        e[2][i]=1107.5760;
        
        f[0][i]=1224.1907;
        f[1][i]=1167.7647;
        f[2][i]=1815.9770;
    }
    
    ini(n, 1e+10, 0.01, 1e-3*1000, a, b, c, d, e, f, p, q, distance, &fails);
    printf("fails: %f\n", fails);
    
    printf("\na1:%f, a2:%f, a3:%f, b1:%f, b2:%f, b3:%f, c1:%f, c2:%f, c3:%f, d1:%f, d2:%f, d3:%f, e1:%f, e2:%f, e3:%f, f1:%f, f2:%f, f3:%f\n", \
           a[0][0], a[0][1], a[0][2], b[0][0], b[0][1], b[0][2], c[0][0], c[0][1], c[0][2], d[0][0], d[0][1], d[0][2], e[0][0], e[0][1], e[0][2], f[0][0], f[0][1], f[0][2]);
    
    return 0;
}

double failureRate(int n, double distance[], double distcmp[], double P1[], double P2[], double Q1[], double Q2[], double epsilon, double fails[])
{
    int counter=0;
    double P[3], Q[3];
    for(int i=0;i<n;i++)
    {
        int fail1 = 0, fail2=0;
        if(distance[i] != distance[i]){
            fails[i] = 1;//nan detected;
            counter++;
        }
        else
        {
            if(fabs(distance[i]-distcmp[i]) > epsilon*2)
            {
                fail1 = 1;
            }
            
            P[0] = (P1[0]) - (P2[0]);
            P[1] = (P1[1]) - (P2[1]);
            P[2] = (P1[2]) - (P2[2]);
            
            Q[0] = (Q1[0]) - (Q2[0]);
            Q[1] = (Q1[1]) - (Q2[1]);
            Q[2] = (Q1[2]) - (Q2[2]);
            
            if(sqrt(DOT(P,P) + DOT(Q,Q)) > epsilon*2)
            {
                fail2 = 1;
            }
            
            if(fail1 && fail2){
                counter++;
                fails[i] = 1;
            }
        }
    }
    return (double(double(counter)/double(n))*100.0);
}



void ini(int n, double r, double eps, double tol, double* a[3], double* b[3], double*  c[3], double*  d[3], double*  e[3], double*  f[3], double*  p[3], double*  q[3], double distance[], double fails[])
{
    double *distcmp = new double[n];
    double A[3],B[3], C[3], D[3], E[3], F[3], P1[3], P2[3], Q1[3], Q2[3];
    
    for(int i=0;i<n;i++)
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
        
        distance[i] = penalty(r, eps, tol, A, B, C, D, E, F, P1, Q1);
        distcmp[i] = bf(A, B, C, D, E, F, P2, Q2);
        
        p[0][i] = P1[0];
        p[1][i] = P1[1];
        p[2][i] = P1[2];
        
        q[0][i] = Q1[0];
        q[1][i] = Q1[1];
        q[2][i] = Q1[2];
        
        printf("PENALTY: p[0]:%g, p[1]:%g, p[2]:%g, q[0]:%g, q[1]:%g, q[2]:%g\n", P1[0], P1[1], P1[2], Q1[0], Q1[1], Q1[2]);
        printf("DISTANCE PENALTY:%g\n", distance[0]);
        printf("BF: p[0]:%g, p[1]:%g, p[2]:%g, q[0]:%g, q[1]:%g, q[2]:%g\n", P2[0], P2[1], P2[2], Q2[0], Q2[1], Q2[2]);
        printf("DISTANCE BF:%g\n", distcmp[0]);
    }
    
    double fa = failureRate(n, distance, distcmp, P1, P2, Q1, Q2, tol, fails);
    delete distcmp;
}

double penalty(double r, double eps, double tol, double A[3], double B[3], double  C[3], double  D[3], double  E[3], double  F[3], double  P[3], double  Q[3])
{
    double BA[3];
    double CA[3];
    double ED[3];
    double FD[3];
    double BC[3];
    double EF[3];
    
    SUB(B,A, BA);
    SUB(C,A, CA);
    SUB(E,D, ED);
    SUB(F,D, FD);
    SUB(B,C, BC);
    SUB(E,F, EF);
    
    double hessian[16], vBA, vCA, vED, vFD; /* hessian columnwise */
    
    vBA = DOT(BA,BA);
    vCA = DOT(CA,CA);
    vED = DOT(ED,ED);
    vFD = DOT(FD,FD);
    
    hessian[0] = 2.*vBA;
    hessian[1] = 2.*DOT(CA,BA);
    hessian[2] = -2.*DOT(ED,BA);
    hessian[3] = -2.*DOT(FD,BA);
    
    hessian[4] = hessian[1]; /* use symmetry */
    hessian[5] = 2.*vCA;
    hessian[6] = -2.*DOT(ED,CA);
    hessian[7] = -2.*DOT(FD,CA);
    
    hessian[8] = hessian[2];
    hessian[9] = hessian[6];
    hessian[10] = 2.*vED;
    hessian[11] = 2.*DOT(FD,ED);
    
    hessian[12] = hessian[3];
    hessian[13] = hessian[7];
    hessian[14] = hessian[11];
    hessian[15] = 2.*vFD;
    
    double x[4] = {0.33, 0.33, 0.33, 0.33}; /* initial guess */
    
    //Declare loop variables;
    double error, dx[4], a[16], SUBXY[3], b[4], h[6], mx[6], dh[24], dmx[4]; double dmx44[16];
    
    double T2dd[3], T1dd[3];
    
    //----Newton loop
    for(int i=0;i<99;i++)
    {
        dh[0] = -1;
        dh[1] = 0;
        dh[2] = 1;
        dh[3] = 0;
        dh[4] = 0;
        dh[5] = 0;
        dh[6] = 0;
        dh[7] = -1;
        dh[8] = 1;
        dh[9] = 0;
        dh[10] = 0;
        dh[11] = 0;
        dh[12] = 0;
        dh[13] = 0;
        dh[14] = 0;
        dh[15] = -1;
        dh[16] = 0;
        dh[17] = 1;
        dh[18] = 0;
        dh[19] = 0;
        dh[20] = 0;
        dh[21] = 0;
        dh[22] = -1;
        dh[23] = 1;
        
        h[0] = -x[0];
        h[1] = -x[1];
        h[2] = x[0]+x[1]-1;
        h[3] = -x[2];
        h[4] = -x[3];
        h[5] = x[2]+x[3]-1; //derivative constraint
        
        for(int j=0;j<6;j++)
        {//masking
            if(h[j]<=0)
            {
                dh[j] = 0;
                dh[6+j] = 0;
                dh[12+j] = 0;
                dh[18+j] = 0;
            }
        }
        
        mx[0] = max(0.0,h[0]);
        mx[1] = max(0.0,h[1]);
        mx[2] = max(0.0,h[2]);
        mx[3] = max(0.0,h[3]);
        mx[4] = max(0.0,h[4]);
        mx[5] = max(0.0,h[5]);
        
        dmx[0] = dh[0] * mx[0] + dh[1] * mx[1] + dh[2] * mx[2] + dh[3] * mx[3] + dh[4] * mx[4] + dh[5] * mx[5];
        dmx[1] = dh[6] * mx[0] + dh[7] * mx[1] + dh[8] * mx[2] + dh[9] * mx[3] + dh[10] * mx[4] + dh[11] * mx[5];
        dmx[2] = dh[12] * mx[0] + dh[13] * mx[1] + dh[14] * mx[2] + dh[15] * mx[3] + dh[16] * mx[4] + dh[17] * mx[5];
        dmx[3] = dh[18] * mx[0] + dh[19] * mx[1] + dh[20] * mx[2] + dh[21] * mx[3] + dh[22] * mx[4] + dh[23] * mx[5];
        
        dmx44[0] = dh[0] * dh[0] + dh[1] * dh[1] + dh[2] * dh[2] + dh[3] * dh[3] + dh[4] * dh[4] + dh[5] * dh[5];
        dmx44[1] = dh[0] * dh[6] + dh[1] * dh[7] + dh[2] * dh[8] + dh[3] * dh[9] + dh[4] * dh[10] + dh[5] * dh[11];
        dmx44[2] = dh[0] * dh[12] + dh[1] * dh[13] + dh[2] * dh[14] + dh[3] * dh[15] + dh[4] * dh[16] + dh[5] * dh[17];
        dmx44[3] = dh[0] * dh[18] + dh[1] * dh[19] + dh[2] * dh[20] + dh[3] * dh[21] + dh[4] * dh[22] + dh[5] * dh[23];
        
        dmx44[4] = dh[6] * dh[0] + dh[7] * dh[1] + dh[8] * dh[2] + dh[9] * dh[3] + dh[10] * dh[4] + dh[11] * dh[5];
        dmx44[5] = dh[6] * dh[6] + dh[7] * dh[7] + dh[8] * dh[8] + dh[9] * dh[9] + dh[10] * dh[10] + dh[11] * dh[11];
        dmx44[6] = dh[6] * dh[12] + dh[7] * dh[13] + dh[8] * dh[14] + dh[9] * dh[15] + dh[10] * dh[16] + dh[11] * dh[17];
        dmx44[7] = dh[6] * dh[18] + dh[7] * dh[19] + dh[8] * dh[20] + dh[9] * dh[21] + dh[10] * dh[22] + dh[11] * dh[23];
        
        dmx44[8] = dh[12] * dh[0] + dh[13] * dh[1] + dh[14] * dh[2] + dh[15] * dh[3] + dh[16] * dh[4] + dh[17] * dh[5];
        dmx44[9] = dh[12] * dh[6] + dh[13] * dh[7] + dh[14] * dh[8] + dh[15] * dh[9] + dh[16] * dh[10] + dh[17] * dh[11];
        dmx44[10] = dh[12] * dh[12] + dh[13] * dh[13] + dh[14] * dh[14] + dh[15] * dh[15] + dh[16] * dh[16] + dh[17] * dh[17];
        dmx44[11] = dh[12] * dh[18] + dh[13] * dh[19] + dh[14] * dh[20] + dh[15] * dh[21] + dh[16] * dh[22] + dh[17] * dh[23];
        
        dmx44[12] = dh[18] * dh[0] + dh[19] * dh[1] + dh[20] * dh[2] + dh[21] * dh[3] + dh[22] * dh[4] + dh[23] * dh[5];
        dmx44[13] = dh[18] * dh[6] + dh[19] * dh[7] + dh[20] * dh[8] + dh[21] * dh[9] + dh[22] * dh[10] + dh[23] * dh[11];
        dmx44[14] = dh[18] * dh[12] + dh[19] * dh[13] + dh[20] * dh[14] + dh[21] * dh[15] + dh[22] * dh[16] + dh[23] * dh[17];
        dmx44[15] = dh[18] * dh[18] + dh[19] * dh[19] + dh[20] * dh[20] + dh[21] * dh[21] + dh[22] * dh[22] + dh[23] * dh[23];
        
        double delta = (hessian[0]+hessian[5]+hessian[10]+hessian[15]) *eps;
        
        a[0] = (hessian[0] + r * dmx44[0]) + 1/pow(2,r);
        a[1] = (hessian[1] + r * dmx44[1]);
        a[2] = (hessian[2] + r * dmx44[2]);
        a[3] = (hessian[3] + r * dmx44[3]);
        
        a[4] = (hessian[4] + r * dmx44[4]);
        a[5] = (hessian[5] + r * dmx44[5]) + 1/pow(2,r);
        a[6] = (hessian[6] + r * dmx44[6]);
        a[7] = (hessian[7] + r * dmx44[7]);
        
        a[8] = (hessian[8] + r * dmx44[8]);
        a[9] = (hessian[9] + r * dmx44[9]);
        a[10] = (hessian[10] + r * dmx44[10]) + 1/pow(2,r);
        a[11] = (hessian[11] + r * dmx44[11]);
        
        a[12] = (hessian[12] + r * dmx44[12]);
        a[13] = (hessian[13] + r * dmx44[13]);
        a[14] = (hessian[14] + r * dmx44[14]);
        a[15] = (hessian[15] + r * dmx44[15]) + 1/pow(2,r);
        
        a[0] = (hessian[0] + r * dmx44[0]) + eps;
        a[1] = (hessian[1] + r * dmx44[1]);
        a[2] = (hessian[2] + r * dmx44[2]);
        a[3] = (hessian[3] + r * dmx44[3]);
        
        a[4] = (hessian[4] + r * dmx44[4]);
        a[5] = (hessian[5] + r * dmx44[5]) + eps;
        a[6] = (hessian[6] + r * dmx44[6]);
        a[7] = (hessian[7] + r * dmx44[7]);
        
        a[8] = (hessian[8] + r * dmx44[8]);
        a[9] = (hessian[9] + r * dmx44[9]);
        a[10] = (hessian[10] + r * dmx44[10]) + eps;
        a[11] = (hessian[11] + r * dmx44[11]);
        
        a[12] = (hessian[12] + r * dmx44[12]);
        a[13] = (hessian[13] + r * dmx44[13]);
        a[14] = (hessian[14] + r * dmx44[14]);
        a[15] = (hessian[15] + r * dmx44[15]) + eps;
        
        SUBXY[0] = (A[0]+(BA[0] * x[0])+(CA[0] * x[1])) - (D[0]+(ED[0] * x[2])+(FD[0] * x[3]));
        SUBXY[1] = (A[1]+(BA[1] * x[0])+(CA[1] * x[1])) - (D[1]+(ED[1] * x[2])+(FD[1] * x[3]));
        SUBXY[2] = (A[2]+(BA[2] * x[0])+(CA[2] * x[1])) - (D[2]+(ED[2] * x[2])+(FD[2] * x[3]));
        
        b[0] = 2*DOT(SUBXY,BA) + r * dmx[0];
        b[1] = 2*DOT(SUBXY,CA) + r * dmx[1];
        b[2] = -2*DOT(SUBXY,ED) + r * dmx[2];
        b[3] = -2*DOT(SUBXY,FD) + r * dmx[3];
        
        solve4(a, b, dx);
        
        // differential of x and y at DX
        T1dd[0] = (BA[0] * dx[0]) + (CA[0] * dx[1]);
        T1dd[1] = (BA[1] * dx[0]) + (CA[1] * dx[1]);
        T1dd[2] = (BA[2] * dx[0]) + (CA[2] * dx[1]);
        
        T2dd[0] = (ED[0] * dx[2]) + (FD[0] * dx[3]);
        T2dd[1] = (ED[1] * dx[2]) + (FD[1] * dx[3]);
        T2dd[2] = (ED[2] * dx[2]) + (FD[2] * dx[3]);
        
        error = sqrt(DOT(T1dd, T1dd)+DOT(T2dd,T2dd));
        if (error < tol) {
            break;
        }
        
        x[0] = x[0] - dx[0];
        x[1] = x[1] - dx[1];
        x[2] = x[2] - dx[2];
        x[3] = x[3] - dx[3];
    }
    
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



double bf(double A[3], double B[3], double C[3], double D[3], double E[3], double F[3], double P[3], double Q[3]){
    //seg to triangle intersection
    int intersection = 0;
    
#if 0
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
    double tp[3], tq[3], tmpP[3], ptmin=1E+30, tmp, ssP[3], ssQ[3], tmpQ[3], ssmin=1E+30;
    
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
    double min=-1;
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

int segt(double p1[3], double p2[3], double A[3], double B[3], double C[3], double P[3]){
    //77 flops
    double u[3], v[3], dir[3], w0[3], n[3], a, b, r, uu, uv, vv, w[3], wu, wv, s, D, t;
    
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


double segseg(double p1[3], double p2[3], double p3[3], double p4[3], double P[3], double Q[3]){
    //69 flops
    double u[3];
    u[0] = p1[0] - p2[0];
    u[1] = p1[1] - p2[1];
    u[2] = p1[2] - p2[2];
    
    double v[3];
    v[0] = p3[0] - p4[0];
    v[1] = p3[1] - p4[1];
    v[2] = p3[2] - p4[2];
    
    double w[3];
    w[0] = p2[0] - p4[0];
    w[1] = p2[1] - p4[1];
    w[2] = p2[2] - p4[2];
    
    double a = DOT(u,u);
    double b = DOT(u,v);
    double c = DOT(v,v);
    double d = DOT(u,w);
    double e = DOT(v,w);
    double D = a*c - b*b;
    double sD = D;
    double tD = D;
    double sN=0, tN=0;
    
    double SMALL_NUM = 0.00000001;
    
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
    double sc;
    //% finally do the division to get sc and tc
    if(fabs(sN) < SMALL_NUM){
        sc = 0.0;
    }else{
        sc = sN / sD;
    }
    double tc;
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
    
    double dP[3];
    dP[0] = w[0] + (sc*u[0]) - (tc*v[0]);
    dP[1] = w[1] + (sc*u[1]) - (tc*v[1]);
    dP[2] = w[2] + (sc*u[2]) - (tc*v[2]);
    
    return sqrt(DOT(dP,dP));
}


double pt(double TP1[3], double TP2[3], double TP3[3], double cPoint[3], double tq[3]){
    //191 flops
    double E0[3];
    E0[0] = TP2[0] - TP1[0];
    E0[1] = TP2[1] - TP1[1];
    E0[2] = TP2[2] - TP1[2];
    
    double E1[3];
    E1[0] = TP3[0] - TP1[0];
    E1[1] = TP3[1] - TP1[1];
    E1[2] = TP3[2] - TP1[2];
    
    double D[3];
    D[0] = TP1[0] - cPoint[0];
    D[1] = TP1[1] - cPoint[1];
    D[2] = TP1[2] - cPoint[2];
    
    double a = DOT(E0,E0);
    double b = DOT(E0,E1);
    double c = DOT(E1,E1);
    double d = DOT(E0,D);
    double e = DOT(E1,D);
    double f = DOT(D,D);
    
    double det = a*c - b*b; //% do we have to use abs here?
    double s   = b*e - c*d;
    double t   = b*d - a*e;
    
    double sqrDistance=0;
    
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
                double invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrDistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
            }
        }
    }else {
        if (s < 0){
            // region 2
            double tmp0 = b + d;
            double tmp1 = c + e;
            if (tmp1 > tmp0){ // minimum on edge s+t=1
                double numer = tmp1 - tmp0;
                double denom = a - 2*b + c;
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
                double tmp0 = b + e;
                double tmp1 = a + d;
                if (tmp1 > tmp0){
                    double numer = tmp1 - tmp0;
                    double denom = a-2*b+c;
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
                double numer = c + e - b - d;
                if (numer <= 0){
                    s = 0;
                    t = 1;
                    sqrDistance = c + 2*e + f;
                }else {
                    double denom = a - 2*b + c;
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

