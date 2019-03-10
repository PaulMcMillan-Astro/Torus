/***************************************************************************//**
\file Compress.h
\brief Contains functions put get (and variants thereof).
Compresses data and allows reading/writing. I don't think it's used any more.
Delete soon.

*                                                                              *
* Compress.h                                                                   *
*                                                                              *
* C++ code created by Walter Dehnen, 1995/96                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  w.dehnen1@physics.ox.ac.uk                                          *
*                                                                              *
* There are TWO VERSIONS, an old and a new one. The uncompressing routines of  *
* the latter are backward compatible, i.e. they can cope with data written by  *
* the old versions of the compressing routines.                                *
*                                                                              *
*******************************************************************************/

#ifndef _Comp_def_
#define _Comp_def_ 1
#include <iostream>

using std::ostream;
using std::istream;
using std::ios;

////////////////////////////////////////////////////////////////////////////////
//
// OLD VERSIONS
//
// Every float is converted in 5 ASCII characters.
// Arrays (1D,2D, and 3D) are written to ostreams in lines of 80 characters.
// One might skip such arrays in which case nothing is read from the file.
//
////////////////////////////////////////////////////////////////////////////////
//
// 1 pure compressing and uncompressing

void  compress  (const float, char[5]);		// compress single float       S
float uncompress(const char[5]);		// uncompress single float     S

void  compress  (const float*,  char*, const int);
void  compress  (const double*, char*, const int);
void  uncompress(char*, float*,  const int);
void  uncompress(char*, double*, const int);

// 2 compressed in and output 

void  put	(const float,                   ostream&);
float get	(                               istream&); 
void  put	(const float*,    const int,    ostream&);
void  get	(      float*,    const int,    istream&);
void  put	(const float**,   const int[2], ostream&);
void  get	(      float**,   const int[2], istream&);
void  put	(const float***,  const int[3], ostream&);
void  get	(      float***,  const int[3], istream&);
void  putFORTRAN(const float**,   const int[2], ostream&);
void  getFORTRAN(      float**,   const int[2], istream&);
void  putFORTRAN(const float***,  const int[3], ostream&);
void  getFORTRAN(      float***,  const int[3], istream&);

void  put	(const double*,   const int,    ostream&);
void  get	(      double*,   const int,    istream&);
void  put	(const double**,  const int[2], ostream&);
void  get	(      double**,  const int[2], istream&);
void  put	(const double***, const int[3], ostream&);
void  get	(      double***, const int[3], istream&);
void  getFORTRAN(      double**,  const int[2], istream&);
void  getFORTRAN(      double***, const int[3], istream&);

// 3 skip: neither read nor uncompress

void skip	(	       istream&);
void skip1D	(const int,    istream&);
void skip2D	(const int[2], istream&);
void skip3D	(const int[3], istream&);

////////////////////////////////////////////////////////////////////////////////
//
// NEW VERSION (March 1997)
//
// Except for the zero (0) every float if converted in 5 ASCII characters. For
// the zero the special character '_' (underscore, 95th character) is reserved.
// For data with many zeros this yields a substantial improvement over the old
// version. However, there are some shortcomings: the size of converted data
// can no longer be easily evaluated from the amount of data, as it depends on
// the number of zeros contained. In particular, skipping is no longer that 
// easy as it was before.
// The uncompressing routines (Get()) are backward compatible,
// i.e. they also work with input created by the old versions above.
//
// There are no routines for pure compressing and uncompressing.
// The new versions use the old versions, so NEVER DELETE THE OLD VERSIONS.
// 
////////////////////////////////////////////////////////////////////////////////
//
// compressed in and output 

void  Put	(const float,                   ostream&);
float Get	(                               istream&); 
void  Put	(const float*,    const int,    ostream&);
void  Get	(      float*,    const int,    istream&);
void  Put	(const float**,   const int[2], ostream&);
void  Get	(      float**,   const int[2], istream&);
void  Put	(const float***,  const int[3], ostream&);
void  Get	(      float***,  const int[3], istream&);

////////////////////////////////////////////////////////////////////////////////
//******************************************************************************
// INLINE FUNCTIONS
//******************************************************************************
////////////////////////////////////////////////////////////////////////////////
// OLD VERSION

inline void  compress(const float* x, char* s, const int n)
{
    int   k,k5;
    for(k=k5=0; k<n; k++,k5+=5) compress(x[k],s+k5);
}

inline void  uncompress(char* s, float* x, const int n)
{
    int   k5;
    float *xk, *xn=x+n;
    for(xk=x,k5=0; xk<xn; xk++,k5+=5) *xk = uncompress(s+k5);
}

inline void skip(istream& in)
{
    in.seekg(5,ios::cur);
}

inline void skip1D(const int n, istream& in)
{
    int n5 = 5*n;
    in.seekg(n5+n5/80,ios::cur);
}

inline void skip2D(const int n[2], istream& in)
{
    int n5 = 5*n[0]*n[1];
    in.seekg(n5+n5/80,ios::cur);
}

inline void skip3D(const int n[3], istream& in)
{
    int n5 = 5*n[0]*n[1]*n[2];
    in.seekg(n5+n5/80,ios::cur);
}

inline float get(istream& in)
{
    char s[5];
    in>>s[0]>>s[1]>>s[2]>>s[3]>>s[4];
    return uncompress(s);
}

inline void get(float* x, const int n, istream& in)
{
    float *xk, *xn=x+n;
    for(xk=x; xk<xn; xk++) *xk = get(in);
}

inline void get(float** x, int const n[2], istream& in)
{
    float **xk, **xn=x+n[0];
    for(xk=x; xk<xn; xk++) get(*xk,n[1],in);
}

inline void get(float*** x, int const n[3], istream& in)
{
    float ***xk, ***xn=x+n[0];
    for(xk=x; xk<xn; xk++) get(*xk,n+1,in);
}

// double presicion versions

inline void  compress(const double* x, char* s, const int n)
{
    int    k,k5;
    for(k=k5=0; k<n; k++,k5+=5) compress(float(x[k]),s+k5);
}

inline void  uncompress(char* s, double* x, const int n)
{
    int    k5;
    double *xk, *xn=x+n;
    for(xk=x,k5=0; xk<xn; xk++,k5+=5) *xk = double(uncompress(s+k5));
}

inline void get(double* x, const int n, istream& in)
{
    double *xk, *xn=x+n;
    for(xk=x; xk<xn; xk++) *xk = double(get(in));
}

inline void get(double** x, int const n[2], istream& in)
{
    double **xk, **xn=x+n[0];
    for(xk=x; xk<xn; xk++) get(*xk,n[1],in);
}

inline void get(double*** x, int const n[3], istream& in)
{
    double ***xk, ***xn=x+n[0];
    for(xk=x; xk<xn; xk++) get(*xk,n+1,in);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// NEW VERSION
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline void Put(const float x, ostream& to)
{
    if(x) put(x,to);
    else  to<<char(95);
}

inline float Get(istream& in)
{
    char s[5];
    in>>s[0];
    if( int(s[0]) == 95 ) return 0.f;
    in>>s[1]>>s[2]>>s[3]>>s[4];
    return uncompress(s);
}

inline void Get(float* x, const int n, istream& in)
{
    float *xk, *xn=x+n;
    for(xk=x; xk<xn; xk++) *xk = Get(in);
}

inline void Get(float** x, int const n[2], istream& in)
{
    float **xk, **xn=x+n[0];
    for(xk=x; xk<xn; xk++) Get(*xk,n[1],in);
}

inline void Get(float*** x, int const n[3], istream& in)
{
    float ***xk, ***xn=x+n[0];
    for(xk=x; xk<xn; xk++) Get(*xk,n+1,in);
}

#endif
