/***************************************************************************//**
\file Vector.h
\brief Contains templated class Vector.
Useful when working with vectors of known size.

*                                                                              *
* Vector.h                                                                     *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Oxford University, Department of Physics, Theoretical Physics                *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail : dehnen@thphys.ox.ac.uk                                              *
*                                                                              *
*******************************************************************************/

#ifndef Vector__
#define Vector__

#include <iostream>
using std::cerr;
#include "Compress.h"
#include <cstdlib>

template<class T, int N>
class Vector {
protected:
    T a[N];
//  void range_error() const;
static void division_by_zero_error();

public:
    Vector() {}
    Vector(const T);
    Vector(const T*);
    Vector(const Vector&);

    Vector&  operator=  (const Vector&);
    Vector&  operator+= (const Vector&);
    Vector&  operator-= (const Vector&);
    Vector&  operator=  (const T);
    Vector&  operator=  (const T*);
    Vector&  operator+= (const T);
    Vector&  operator-= (const T);
    Vector&  operator*= (const T);
    Vector&  operator/= (const T);
    Vector&  apply	(T(*)(T));

    Vector   operator-  () const;
    Vector   operator+  (const Vector&) const;
    Vector   operator-  (const Vector&) const;
    T        operator*  (const Vector&) const;
    int      operator== (const Vector&) const;
    int      operator!= (const Vector&) const;
    Vector   operator+  (const T) const;
    Vector   operator-  (const T) const;
    Vector   operator*  (const T) const;
    Vector   operator/  (const T) const;

    Vector&  multiply_elements (const Vector&);

    T        norm       () const;
    T        operator() (const int n) const { return  a[n]; }
    T&       operator[] (const int n)       { return  a[n]; }
    int      NumberofTerms() const { return N; }

	     operator T*       ()	 { return a; }
	     operator const T* () const	 { return a; }
};

// Vector cross product for triples

template<class T>
Vector<T,3> operator^ (const Vector<T,3>& A, const Vector<T,3>& B)
{
    Vector<T,3> X;
    X[0] = A(1)*B(2) - A(2)*B(1);
    X[1] = A(2)*B(0) - A(0)*B(2);
    X[2] = A(0)*B(1) - A(1)*B(0);
    return X;
}

// norm defined as non-member function

template<class T, int N>
inline T norm(const Vector<T,N>& V)
    { return V.norm(); }

/*
// operator= between Vectors of same size but different types

template<class T, class S, int N>
inline Vector<T,N>& operator= (Vector<T,N>& Vt, const Vector<S,N>& Vs)
    { for(i=0; i<N; i++) Vt[i] = T(Vs(i));
      return Vt; }
*/

// arithmetic operators with scalar on the left

template<class T, int N>
inline Vector<T,N> operator+ (const T x, const Vector<T,N>& V)
    { return V+x; }
template<class T, int N>
inline Vector<T,N> operator* (const T x, const Vector<T,N>& V)
    { return V*x; }
template<class T, int N>
inline Vector<T,N> operator- (const T x, const Vector<T,N>& V)
    { Vector<T,N> P(x); return P-=V; }

// compressed in and output

template<int N>
inline void put(const Vector<float,N>& V, ostream& to)
    { put((float*)V,N,to); }

template<int N>
inline void put(const Vector<double,N>& V, ostream& to)
    { put((double*)V,N,to); }

template<int N>
inline void put(const Vector<float,N>* V, const int M, ostream& to)
    { float* A = new float[M*N];
      float *Ai=A;
      int   m,n;
      for(m=0; m<M; m++) for(n=0; n<N; n++,Ai++) *Ai = V[m](n);
      ::put(A, M*N, to);
      delete[] A; }

template<int N>
inline void put(const Vector<double,N>* V, const int M, ostream& to)
    { double* A = new double[M*N];
      double *Ai=A;
      int    m,n;
      for(m=0; m<M; m++) for(n=0; n<N; n++,Ai++) *Ai = V[m](n);
      ::put(A, M*N, to);
      delete[] A; }

template<int N>
inline void get(Vector<float,N>& V, istream& from)
    { get((float*)V,N,from); }

template<int N>
inline void get(Vector<double,N>& V, istream& from)
    { get((double*)V,N,from); }

template<int N>
inline void get(Vector<float,N>* V, const int M, istream& from)
    { Vector<float,N> *Vm=V, *VM=V+M;
      for(; Vm<VM; Vm++) get(*Vm,N,from); }

template<int N>
inline void get(Vector<double,N>* V, const int M, istream& from)
    { Vector<double,N> *Vm=V, *VM=V+M;
      for(; Vm<VM; Vm++) get(*Vm,N,from); }

template<class T, int N>
inline Vector<T,N> read_compressed(istream& from)
    { Vector<T,N> x;
      x.get(from);
      return x; }

////////////////////////////////////////////////////////////////////////////////
// protected member functions:

template<class T, int N>
void Vector<T,N>::division_by_zero_error()
{ 
    cerr << " Vector: division by zero \n";
#ifndef ebug
    exit(1); 
#endif
}

////////////////////////////////////////////////////////////////////////////////
// public member functions:
// (A) constructors

template<class T, int N>
Vector<T,N>::Vector(const T fill_value)
    { for(int i=0; i<N; i++) a[i] = fill_value; }

template<class T, int N>
Vector<T,N>::Vector(const T *array)
    { for(int i=0; i<N; i++) a[i] = array[i]; }

template<class T, int N>
Vector<T,N>::Vector(const Vector<T,N>& V)
    { for(int i=0; i<N; i++) a[i] = V.a[i]; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator= (const Vector<T,N>& V)
    { for(int i=0; i<N; i++) a[i] = V.a[i];
      return *this; }

// (B) arithmetic operators with assign: =, +=, -=, *=, /=

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator= (const T fill_value)
    { for(int i=0; i<N; i++) a[i] = fill_value;
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator= (const T* array)
    { for(int i=0; i<N; i++) a[i] = array[i];
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator+= (const Vector<T,N>& V)
    { for(int i=0; i<N; i++) a[i] += V.a[i];
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator-= (const Vector<T,N>& V)
    { for(int i=0; i<N; i++) a[i] -= V.a[i];
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator+= (const T m)
    { for(int i=0; i<N; i++) a[i] += m;
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator-= (const T m)
    { for(int i=0; i<N; i++) a[i] -= m;
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator*= (const T m)
    { for(int i=0; i<N; i++) a[i] *= m;
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::operator/= (const T m)
    { if(m==T(0.)) division_by_zero_error();
      for(int i=0; i<N; i++) a[i] /= m;
      return *this; }

template<class T, int N>
Vector<T,N>& Vector<T,N>::multiply_elements (const Vector& V)
    { for(int i=0; i<N; i++) a[i] *= V.a[i];
      return *this; }

// (C) application of functions T->T onto individual elements

template<class T, int N>
Vector<T,N>& Vector<T,N>::apply	( T(*f)(T) )
    { for(int i=0; i<N; i++) a[i] = f(a[i]);
      return *this; }

// (D) further arithmetic operators, all constant member functions

template<class T, int N>
Vector<T,N> Vector<T,N>::operator- () const
    { Vector<T,N> P(T(0)); return P-=*this; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator+ (const Vector<T,N>& V) const
    { Vector<T,N> P(*this); return P+=V; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator- (const Vector<T,N>& V) const
    { Vector<T,N> P(*this); return P-=V; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator+ (const T x) const
    { Vector<T,N> P(*this); return P+=x; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator- (const T x) const
    { Vector<T,N> P(*this); return P-=x; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator* (const T x) const
    { Vector<T,N> P(*this); return P*=x; }

template<class T, int N>
Vector<T,N> Vector<T,N>::operator/ (const T x) const
    { Vector<T,N> P(*this); return P/=x; }

template<class T, int N>
T Vector<T,N>::operator* (const Vector<T,N>& V) const
    { T x=a[0] * V.a[0];
      for(int i=1; i<N; i++) x += a[i] * V.a[i];
      return x; }

template<class T, int N>
int Vector<T,N>::operator== (const Vector<T,N>& V) const
    { for(int i=0; i<N; i++) if(a[i] != V.a[i]) return 0;
      return 1; }

template<class T, int N>
int Vector<T,N>::operator!= (const Vector<T,N>& V) const
    { for(int i=0; i<N; i++) if(a[i] != V.a[i]) return 1;
      return 0; }


// (E) norm
// For T other than int, long int, float, double, and long double
//	`T norm(const T)'
// must be pre-defined by the user.

template<class T, int N>
T Vector<T,N>::norm() const 
    { T x = a[0]*a[0];
      for(int i=1; i<N; i++) x += a[i]*a[i];
      return x; }

////////////////////////////////////////////////////////////////////////////////
// class Vector: friends for stream I/O

template<class T, int N>
ostream& operator<< (ostream& s, const Vector<T,N>& V)
    { s << V(0);
      for(int i=1; i<N; i++) s << ' ' << V(i);
      return s; }
/* version that writes `(a,b,c)'
    { s << '(' << V(0);
      for(int i=1; i<N; i++) s << ',' << V(i);
      s << ')';
      return s; }
*/

template<class T, int N>
istream& operator>> (istream& s, Vector<T,N>& V)
    { T x[N];
      char c=0;
      int i;
      s >> c;
      if(c == '(') {
          for(i=0; i<N; i++) s >> x[i];
	  s >> c;
	  if(c != ')') s.clear(ios::badbit);
      } else {
          s.putback(c);
          for(i=0; i<N; i++) s >> x[i];
      }
      if(s) for(i=0; i<N; i++) V[i] = x[i];
      return s; }

#endif

#ifdef  __COMPLEX__
#ifndef VectorComplex__
#define VectorComplex__

template<class T, int N>
inline Vector<T,N> real(const Vector<complex<T>,N>& V)
    { Vector<T,N> R;
      for(int i=0; i<N; i++) R[i] = real(V(i)); 
      return R; }
template<class T, int N>
inline Vector<T,N> imag(const Vector<complex<T>,N>& V)
    { Vector<T,N> R;
      for(int i=0; i<N; i++) R[i] = imag(V(i)); 
      return R; }
template<class T, int N>
inline 	Vector<T,N> arg(const Vector<complex<T>,N>& V)
    { Vector<T,N> R;
      for(int i=0; i<N; i++) R[i] = arg(V(i)); 
      return R; }
template<class T, int N>
inline 	Vector<complex<T>,N> conj(const Vector<complex<T>,N>& V)
    { Vector<complex<T>,N> R;
      for(int i=0; i<N; i++) R[i] = conj(V(i)); 
      return R; }

#endif
#endif
