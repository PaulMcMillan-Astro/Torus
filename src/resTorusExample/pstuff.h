// Header required to interface with pstuff

#include "jjb_utils.h"

#ifndef _JJBpress_
#define _JJBpress_ 1


#define SIGN(A,B) ((B)>=0 ?fabs(A):-fabs(A))
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


void spline(double *,double *,int,double,double,double *);
double splint(double *, double *, double *,int,double);
double splintp(double *,double *,double *,int,double);
double splintp2(double *xa, double *ya, double *y2a, int np, double x);


void rk4(double *,double *,int,double,double, void(*)(double,double *,double *));
void rkck(double y[], double dydx[],int N, double x, double h, double yout[],
	  double yerr[], void (*derivs)(double, double [], double []));

void rkqs(double *,double *,int,double *,double,double,
                  double *,double *, double *,void (*)(double, double*, double*));
double zbrent(double *pars,double (*func)(double*,double), double x1, double x2,double fa,double fb,double tol,int itmax);

void four1(float*,unsigned long,int);
void four1(double*,unsigned long,int);
void fourn(float*,unsigned long*,int,int);
//void fourn(double*,unsigned long*,int,int); - now in jjb_utils

#endif
