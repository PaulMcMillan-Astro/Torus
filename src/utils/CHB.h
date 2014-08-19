/*******************************************************************************
*                                                                              *
*  CHB.h                                                                       *
*                                                                              *
* C++ code written by Paul McMillan, 2008                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
/*

This is code which defines the class Cheby, which deals with Chebyshev 
polynomials with up to 50 coefficients. It can take the coefficients from 
input, or fit them from data (also input).

Then it can provide the value of the Cheb. polynomial C(x) at any x, also C' 
and C''.

*/

#ifndef _Chebyshev_
#define _Chebyshev_ 1

#include <iostream>

class Cheby {
  void getchebc();
 public:
  int NChb, tabsize;    // Number of coefficients, size of auxiliary table  
  double *s, *e1;       // Auxiliary table, table of coeffs.
  Cheby                 ();             
  Cheby                 (const Cheby&); // copy other one
  ~Cheby                ();
  Cheby&    operator=   (const Cheby&); // copy other one 
  Cheby                 (const int);    // just set up table of coefficients
  Cheby                 (double *, const int); // Already fit
  Cheby                 (double *, double *,const int,const int); // y,x,np,NC
  // fit Chebyshev polynomial with NC coefficients to y(x), defined at np points
 
  void      writecoeffs (std::ostream&) const;            // output coeffs 
  void      setcoeffs   (double *, const int);            // input prefit coeffs
  void      chebyfit    (double*, double*, const int, const int=0);// y,x,np,NC
  // fit Chebyshev polynomial with NC coefficients to y(x), defined at np points
  // if NC undefined on input (i.e.=0), then NChb unchanged.

  void      unfitderiv  (const double, double&, double&, double&) const;
  // For Chebshev polynomial C(x) w. input x, these are: x, C(x), C'(x), C''(x)
  void      unfitderiv  (const double, double&, double&)          const;
  // For Chebshev polynomial C(x) w. input x, these are: x, C(x), C'(x)
  void      unfitn      (double *, double *, const int)           const;
  // Input table x[n], output table C(x)[n], number of terms n.
  double    unfit1      (const double)                            const;
  // given x, returns C(x)
};




#endif
