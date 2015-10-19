/*******************************************************************************
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2006-07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
#include "IsochronePot.h"
#include <cmath>

double IsochronePotential::operator() (const double R, const double z) const
{
  double rsq=R*R+z*z;
  return -GM/(b+sqrt(bsq+rsq));
}

double IsochronePotential::operator() (const double R, const double z, double& dPdR,
			       double& dPdz) const
{
  double rsq = R*R+z*z,
    r = sqrt(rsq), ir = 1./r,
    tmp1 = sqrt(bsq+rsq),
    denom = 1./(b+tmp1),
    dPdr = GM*r*denom*denom/tmp1;
  dPdR = dPdr*R/r;
  dPdz = dPdr*z/r;
  return -GM*denom;

}


Frequencies IsochronePotential::KapNuOm(const double R) const
{
  //epi[0] = sqrt(d2PdRR+3*dPdR/R)
  //epi[1] = sqrt(d2Pdzz)
  //epi[2] = sqrt(dPdR/R);

  double rsq = R*R,
    r = R, ir = 1./r,
    tmp1 = sqrt(bsq+rsq),
    denom = 1./(b+tmp1),
    dPdr = GM*r*denom*denom/tmp1,
    dPdR = dPdr,
    d2Pdrr = dPdr* ( ir - 2*r*denom/tmp1 - r/(tmp1*tmp1) ),
    d2PdRR = d2Pdrr,
    d2Pdzz = dPdr/r;

  Frequencies epi;
  epi[0] = sqrt(d2PdRR+3.*dPdR/R);
  epi[1] = sqrt(d2Pdzz);
  epi[2] = sqrt(dPdR/R);

  return epi;

}

ostream& operator<< (ostream& to, const IsochronePotential& P)
{
  to << "Isochrone potential with M=" << P.GM/(Units::G) << " b=" << P.b << '\n';
    return to;
}








