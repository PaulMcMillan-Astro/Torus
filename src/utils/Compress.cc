/*******************************************************************************
*                                                                              *
* Compress.cc                                                                  *
*                                                                              *
* C++ code created by Walter Dehnen, 1995/96                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  w.dehnen1@physics.ox.ac.uk                                          *
*                                                                              *
*******************************************************************************/

#include "Compress.h"
#include "Numerics.h"
#include "FreeMemory.h"

////////////////////////////////////////////////////////////////////////////////
//
// routines of the OLD VERSION
// 
// NEVER DELETE THEM, because they are used by the routines of the NEW VERSION
//

float uncompress(const char s[5])
// uncompress single float
{
    const int      i1=90;
    const float    at=8388608.f;
    int   m4=s[4];
    if(m4>=124) m4--; if(m4>=95) m4--; m4-=33;
    if(m4==80) return 0.f;
    int   j, nex, m0=s[0], m1=s[1], m2=s[2], m3=s[3];
    if(m0>=124) m0--; if(m0>=95) m0--; m0-=33;
    if(m1>=124) m1--; if(m1>=95) m1--; m1-=33;
    if(m2>=124) m2--; if(m2>=95) m2--; m2-=33;
    if(m3>=124) m3--; if(m3>=95) m3--; m3-=33;
    float x=1.f;
    if(m4>=40) {
        m4-= 40;
        x =-1.f;
    }
    j   = m4/12;
    m4 -= j*12;
    nex = j*i1 + m3 - 129;
    x  *= float(((m4*i1+m2)*i1+m1)*i1+m0)/at + 1.f;
    return x * pow(2.f,nex);
}

void compress(const float x, char s[5])
// compress single float
{
    const int    i1=90, i2=8100, i3=729000;
    const float  at=8388608., small=1.e-28, big=1.e28;
    int m0,m1,m2,m3,m4;
    if(x==0.f) {
	m4=80; m3=m2=m1=m0=0;
    } else {
        int   j,lx,rem;
        float xk=x, ax=fabs(xk), xa;
	if(ax>big) {
	    xk = sign(big,xk);
	    ax = big;
	} else if(ax<small) {
	    xk = sign(small,xk);
	    ax = small;
	}
  	if(ax>=1.) for(xa=ax,lx=128; xa>=1.f; xa*=0.5f,lx++) ;
  	else       for(xa=ax,lx=129; xa< 1.f; xa*=2.f ,lx--) ;
        j   = lx/i1;
	m3  = lx-j*i1;
        lx -= 129;
	rem = int((ax*pow(2.,-lx)-1.)*at);
        m4  = rem/i3;
	rem-= m4*i3;
	m2  = rem/i2;
	rem-= m2*i2;
	m1  = rem/i1;
	m0  = rem-m1*i1;
	m4 += j*12;
	if(xk<0) m4 += 40;
    }
    m0+=33; if(m0>94) m0++; if(m0>123) m0++; s[0]=char(m0);
    m1+=33; if(m1>94) m1++; if(m1>123) m1++; s[1]=char(m1);
    m2+=33; if(m2>94) m2++; if(m2>123) m2++; s[2]=char(m2);
    m3+=33; if(m3>94) m3++; if(m3>123) m3++; s[3]=char(m3);
    m4+=33; if(m4>94) m4++; if(m4>123) m4++; s[4]=char(m4);
}

void put(const float x, ostream& to)
// write single float in compressed format
{
    const int    i1=90, i2=8100, i3=729000;
    const float  at=8388608., small=1.e-28, big=1.e28;
    int m0,m1,m2,m3,m4;
    if(x==0.f) {
	m4=80; m3=m2=m1=m0=0;
    } else {
        int   j,lx,rem;
        float xk=x, ax=fabs(xk), xa;
	if(ax>big) {
	    xk = sign(big,xk);
	    ax = big;
	} else if(ax<small) {
	    xk = sign(small,xk);
	    ax = small;
	}
  	if(ax>=1.) for(xa=ax,lx=128; xa>=1.f; xa*=0.5f,lx++) ;
  	else       for(xa=ax,lx=129; xa< 1.f; xa*=2.f ,lx--) ;
        j   = lx/i1;
	m3  = lx-j*i1;
        lx -= 129;
	rem = int((ax*pow(2.,-lx)-1.)*at);
        m4  = rem/i3;
	rem-= m4*i3;
	m2  = rem/i2;
	rem-= m2*i2;
	m1  = rem/i1;
	m0  = rem-m1*i1;
	m4 += j*12;
	if(xk<0) m4 += 40;
    }
    m0+=33; if(m0>94) m0++; if(m0>123) m0++; to<<char(m0);
    m1+=33; if(m1>94) m1++; if(m1>123) m1++; to<<char(m1);
    m2+=33; if(m2>94) m2++; if(m2>123) m2++; to<<char(m2);
    m3+=33; if(m3>94) m3++; if(m3>123) m3++; to<<char(m3);
    m4+=33; if(m4>94) m4++; if(m4>123) m4++; to<<char(m4);
}

void put(const float* x, const int n, ostream& to)
// put 1D array (rows of 80)
{
    int   k,C=0;
    for(k=0; k<n; C++,k++) {
        put(x[k],to);
        if(C>=15) {
            C=-1;
            to<<'\n';
        }
    }
}

void put(const float** x, int const n[2], ostream& to)
// put 2D array (rows of 80)
{
    int   i,j,C=0;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; C++,j++) {
	    put(x[i][j],to);
            if(C>=15) {
                C=-1;
                to<<'\n';
            }
        }
}

void put(const float*** x, int const n[3], ostream& to)
// put 3D array (rows of 80)
{
    int   i,j,k,C=0;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
            for(k=0; k<n[2]; C++,k++) {
	        put(x[i][j][k],to);
                if(C>=15) {
                    C=-1;
                    to<<'\n';
                }
            }
}

void putFORTRAN(const float** x, int const n[2], ostream& to)
// put 2D array (rows of 75), order of indices as in FORTRAN
{
    int   i,j,C=0;
    for(j=0; j<n[1]; j++)
        for(i=0; i<n[0]; C++,i++) {
	    put(x[i][j],to);
            if(C>=14 && (i<n[0]-1 || j<n[1]-1) ) {
                C=-1;
                to<<'\n';
            }
        }
}

void putFORTRAN(const float*** x, int const n[3], ostream& to)
// put 3D array (rows of 75), order of indices as in FORTRAN
{
    int   i,j,k,C=0;
    for(k=0; k<n[2]; k++)
    	for(j=0; j<n[1]; j++)
    	    for(i=0; i<n[0]; C++,i++) {
	        put(x[i][j][k],to);
                if(C>=14 && (i<n[0]-1 || j<n[1]-1 || k<n[2]-1) ) {
                    C=-1;
                    to<<'\n';
                }
            }
}

void getFORTRAN(float**x, int const n[2], istream& in)
// get 2D array written by FORTRAN routines (eg. of James Binney)
// The difference is in the order of indices
{
    int   m[2]={n[1],n[0]};
    float **y;
    Alloc2D(y,m);
    get(y,m,in);
    int i,j;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
	    x[i][j] = y[j][i];
    Free2D(y);
}

void getFORTRAN(float***x, const int n[3], istream& in)
// get 3D array written by FORTRAN routines (eg. of James Binney)
// The difference is in the order of indices
{
    int   m[3]={n[2],n[1],n[0]};
    float ***y;
    Alloc3D(y,m);
    get(y,m,in);
    int i,j,k;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
            for(k=0; k<n[2]; k++)
	        x[i][j][k] = y[k][j][i];
    Free3D(y);
}

void put(const double* x, const int n, ostream& to)
// put 1D array (rows of 80)
{
    int    k,C=0;
    for(k=0; k<n; C++,k++) {
        put(float(x[k]),to);
        if(C>=15) {
            C=-1;
            to<<'\n';
        }
    }
}

void put(const double** x, int const n[2], ostream& to)
// put 2D array (rows of 75)
{
    int   i,j,C=0;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; C++,j++) {
	    put(float(x[i][j]),to);
            if(C>=15) {
                C=-1;
                to<<'\n';
            }
        }
}

void put(const double*** x, int const n[3], ostream& to)
// put 3D array (rows of 75)
{
    int   i,j,k,C=0;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
            for(k=0; k<n[2]; C++,k++) {
	        put(float(x[i][j][k]),to);
                if(C>=15) {
                    C=-1;
                    to<<'\n';
                }
            }
}

void getFORTRAN(double**x, int const n[2], istream& in)
// get 2D array written by FORTRAN routines (eg. of James Binney)
// The difference is in the order of indices
{
    int    m[2]={n[1],n[0]};
    double **y;
    Alloc2D(y,m);
    get(y,m,in);
    int i,j;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
	    x[i][j] = y[j][i];
    Free2D(y);
}

void getFORTRAN(double***x, int const n[3], istream& in)
// get 3D array written by FORTRAN routines (eg. of James Binney)
// The difference is in the order of indices
{
    int    m[3]={n[2],n[1],n[0]};
    double ***y;
    Alloc3D(y,m);
    get(y,m,in);
    int i,j,k;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
            for(k=0; k<n[2]; k++)
	        x[i][j][k] = y[k][j][i];
    Free3D(y);
}

////////////////////////////////////////////////////////////////////////////////
//
// routines of the NEW VERSION

void Put(const float* x, const int n, ostream& to)
// put 1D array (rows of <= 80)
{
    int   k,C;
    float xk;
    for(C=k=0; k<n; k++) {
	xk = x[k];
	if(xk==0.f) {				// if x == 0
	    if(C>79) { to<<'\n'; C=0; }		//  line already >79 ? -> \n
	    to << char(95);			//  write out '_'
	    C++;				//  line 1 char longer
	} else {				// otherwise: x != 0
	    if(C>75) { to<<'\n'; C=0; }		//  line already >75 ? -> \n
	    put(xk,to);				//  write out compressed
	    C += 5;				//  line 5 chars longer
	}
    }
}

void Put(const float** x, int const n[2], ostream& to)
// put 2D array (rows of 75)
{
    int   i,j,C=0;
    float xk;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++) {
	    xk = x[i][j];
	    if(xk==0.f) {			// if x == 0
	        if(C>79) { to<<'\n'; C=0; }	//  line already >79 ? -> \n
	        to << char(95);			//  write out '_'
	        C++;				//  line 1 char longer
	    } else {				// otherwise: x != 0
	        if(C>75) { to<<'\n'; C=0; }	//  line already >75 ? -> \n
	        put(xk,to);			//  write out compressed
	        C += 5;				//  line 5 chars longer
	    }
        }
}

void Put(const float*** x, int const n[3], ostream& to)
// put 3D array (rows of 75)
{
    int   i,j,k,C=0;
    float xk;
    for(i=0; i<n[0]; i++)
        for(j=0; j<n[1]; j++)
            for(k=0; k<n[2]; k++) {
		xk = x[i][j][k];
		if(xk==0.f) {			// if x == 0
		    if(C>79) { to<<'\n'; C=0; }	//  line already >79 ? -> \n
		    to << char(95);		//  write out '_'
		    C++;			//  line 1 char longer
		} else {			// otherwise: x != 0
		    if(C>75) { to<<'\n'; C=0; }	//  line already >75 ? -> \n
		    put(xk,to);			//  write out compressed
		    C += 5;			//  line 5 chars longer
		}
            }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* Here are the very old (and working) versions, just for reference.
void compress(float* x, char* s, const int npt)
{
    const int      i1=90, i2=8100, i3=729000;
    const float    at=8388608., small=1.e-28, big=1.e28;
    int            m[5];
    int   j,k,lx,rem;
    float xk,ax,xa;
    for(k=0;k<npt;k++) {
	if(x[k]==0.) {
	    m[4]=80; m[3]=0; m[2]=0; m[1]=0; m[0]=0;
	} else {
	    xk = x[k];
	    ax = abs(xk);
	    if(ax>big)            {xk=sign(big,xk);   ax=big;}
	        else if(ax<small) {xk=sign(small,xk); ax=small;}
	    if(ax>=1.) for(xa=ax,lx=128; xa>=1.f; xa*=0.5f,lx++) ;
	        else   for(xa=ax,lx=129; xa< 1.f; xa*=2.f ,lx--) ;
            j    = lx/i1;
	    m[3] = lx-j*i1;
            lx  -= 129;
	    rem  = int((ax*pow(2.,-lx)-1.)*at);
            m[4] = rem/i3;
	    rem -= m[4]*i3;
	    m[2] = rem/i2;
	    rem -= m[2]*i2;
	    m[1] = rem/i1;
	    m[0] = rem - m[1]*i1;
	    m[4]+= j*12;
	    if(xk<0) m[4] += 40;
	}
	for(j=0;j<5;j++) {
	    m[j] += 33;
	    if(m[j]>=94)  m[j]++;
	    if(m[j]>=123) m[j]++;
	    s[5*k+j] = char(m[j]);
	}
    }
}

void uncompress(char* s, float* x, const int npt)
{
    const int    i1=90;
    const float  at=8388608.;
    int          m[5];
    int j,k,nex;
    for(k=0; k<npt; k++) {
	for(j=0;j<5;j++) {
	    m[j] = s[5*k+j];
	    if(m[j]>=124) m[j]--;
	    if(m[j]>=95)  m[j]--;
	    m[j] -= 33;
	}
	if(m[4]==80) x[k]=0.;
	else {
            if(m[4]>=40) { m[4]-= 40; x[k] =-1.; }
	        else       x[k] = 1.;
            j    = m[4]/12;
	    m[4]-= j*12;
	    nex  = j*i1 + m[3] - 129;
	    x[k]*= float(((m[4]*i1+m[2])*i1+m[1])*i1+m[0])/at + 1.;
	    x[k]*= pow(2.f,nex);
	}
    }
}
*/
