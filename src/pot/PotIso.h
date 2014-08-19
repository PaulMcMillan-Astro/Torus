/*******************************************************************************
*                                                                              *
*  IsoPot.h                                                                    *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2006/07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
//
//  Gives the Isochrone potential Phi = -GM/(b+sqrt(b^2+r^2))
//

#ifndef _PotIso_h_
#define _PotIso_h_ 1

#include <iostream>
#include "Potential.h"
#include "Units.h"

class IsoPot : public Potential {
private:
  const double GM, b;
  const double bsq, iGM;
public:
    IsoPot(const double=1., const double=1.);
    double operator() (const double, const double) const;
    double operator() (const double, const double, double&, double&) const;
    double RfromLc(const double, double*) const;
    double LfromRc(const double, double*) const;
    Frequencies KapNuOm(const double) const;
    friend ostream& operator<< (ostream&, const IsoPot&);
};

inline IsoPot::IsoPot(const double M, const double B) 
	      : GM(M*Units::G), b(B), bsq(b*b), iGM(1./GM)
{ }


#endif
