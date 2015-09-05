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

/** \brief Axisymmetric Stackel potential for an oblate spheroid.

Input parameters are rho_0, Rs, q. 

Potential is for a body with density rho = rho_0/(1+m^2)^2, where: 

m^2 = (R^2 + (z/q)^2) / Rs^2

See De Zeeuw 1985 (MNRAS, 216, 273) for details (Section 3.4.2)

 */

class StackelPotential : public Potential {
private:
  const double rh0, alpha, gamma,Gfac;
public:
    StackelPotential(const double=1., const double=1., 
	       const double=1.);
    double operator() (const double, const double) const;
    double operator() (const double, const double, double&, double&) const;
    //double operator() (const double, const double,
    //		       double&, double&, double&, double&, double&) const;
    //double RfromLc(const double, double*) const;
    //double LfromRc(const double, double*) const;
    Frequencies KapNuOm(const double) const;
    friend ostream& operator<< (ostream&, const StackelPotential&);
};

inline StackelPotential::StackelPotential(const double rho0, const double Rs,
			      const double q) 
		  : rh0(rho0), alpha(-Rs*Rs), gamma(q*q*alpha), 
		  Gfac(TPi*(Units::G)*rh0*(-alpha))
{ }


#endif
