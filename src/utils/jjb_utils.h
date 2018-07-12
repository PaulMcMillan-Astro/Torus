#ifndef _jjb_utils_
#define _jjb_utils_
#include "../Torus.h"
#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SGN(A)  ((A)>0? 1:-1)
double *dmatrix(int);
double **dmatrix(int,int);
double ***dmatrix(int,int,int);
double ****dmatrix(int,int,int,int);
double ***d3array(int,int,int);
void free_d3array(double***);
void delmatrix(double**,int);
void delmatrix(double***,int,int);
void delmatrix(double****,int,int,int);
void tidy_angles(Angles&);
int quadratic(double,double,double,double*);
double lntp(double*,double*,int,double);
void fourn(double*,unsigned long*,int,int);
void rlft3(double***,double**,unsigned long,unsigned long,unsigned long,int);
#endif
