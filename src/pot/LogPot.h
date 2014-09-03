/***************************************************************************//**
\file LogPot.h
\brief Contains class LogPotential. 
Logarithmic potential.
									     
*                                                                              *
* LogPot.h                                                                     *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2007-
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#ifndef LogPotential_def_
#define LogPotential_def_ 1

#include <iostream>
#include <cmath>
#include "Potential.h"

class LogPotential : public Potential {
public:
  double q, qi, q2i, v0, v0sq, v0sqhalf, rc, rc2, plusconst, Rmax;
    double Rei, Rei3;
    void   error(const char*) const;
    //public:
    LogPotential(const double=1., const double=1., const double=1., const double=0.);
    double operator() (const double, const double) const;
    double operator() (const double, double&, double&) const;
    double operator() (const double, const double, double&, double&) const;
    double operator() (const double, const double,
		       double&, double&, double&, double&, double&) const;
    double RfromLc(const double, double* = 0) const;
    double LfromRc(const double, double* = 0) const;
    Frequencies KapNuOm(const double) const;
friend ostream& operator<< (ostream&, const LogPotential&);
};

inline LogPotential::LogPotential(const double V0, const double Q,
		      const double Rc, const double Ri) :
q(Q), qi(1./q), q2i(qi*qi), v0(V0), v0sq(v0*v0), v0sqhalf(0.5*v0sq),
	      rc(Rc), rc2(rc*rc), Rei(Ri), Rei3(Rei*Rei*Rei), Rmax(500.) 
	      
{ 
  plusconst = -v0sqhalf * log(Rmax*Rmax+q2i+rc2);
}

#endif
