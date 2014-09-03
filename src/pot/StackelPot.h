/***************************************************************************//**
\file StackelPot.h  
\brief Contains class StackelPotential. 
An axisymmetric Stackel potential.

*                                                                              *
*  StackelPot.h                                                                *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2006/07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#ifndef _Stack_h_
#define _Stack_h_ 1
   
#include <iostream>
#include "Potential.h"
#include "Units.h"

class StackelPotential : public Potential {
private:
  const double rh0, alpha, gamma,Gfac;
  //    double q, qi, q2i, v0, v0sq, v0sqhalf, rc, rc2;
  // double Rei, Rei3;
  // double findG (const double);
  //double findG (const double, double&);
  //void   error(const char*) const;
public:
    StackelPotential(const double=1., const double=1., 
	       const double=1.);
    double operator() (const double, const double) const;
    //double operator() (const double, double&, double&) const;
    double operator() (const double, const double, double&, double&) const;
    //double operator() (const double, const double,
    //		       double&, double&, double&, double&, double&) const;
    double RfromLc(const double, double*) const;
    double LfromRc(const double, double*) const;
    Frequencies KapNuOm(const double) const;
    friend ostream& operator<< (ostream&, const StackelPotential&);
};

inline StackelPotential::StackelPotential(const double rho0, const double Rs,
			      const double q) 
		  : rh0(rho0), alpha(-Rs*Rs), gamma(q*q*alpha), 
		  Gfac(TPi*(Units::G)*rh0*(-alpha))
{ }


#endif
