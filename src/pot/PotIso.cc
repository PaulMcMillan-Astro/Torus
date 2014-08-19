/*******************************************************************************
*                                                                              *
* PotIso.cc                                                                    *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2006-07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
#include "PotIso.h"
#include <cmath>

double IsoPot::operator() (const double R, const double z) const
{
  double rsq=R*R+z*z;
  return -GM/(b+sqrt(bsq+rsq));
}

double IsoPot::operator() (const double R, const double z, double& dPdR,
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

//double IsoPot::operator() (const double R, 
//                           double& dPdR, double& d2PdRR) const {}

//double IsoPot::operator() (const double R, const double z,
//			   double& dPdR, double& dPdz,
//			   double& d2PdRR, double& d2Pdzz, double& d2PdRz) const
//{}


double IsoPot::RfromLc(const double L, double* dR) const
{
  bool more=false;
  double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL;
  R=exp(lR);
  P= (*this)(R,0.,dPR,dPz);
  LcR=pow(R*R*R*dPR,0.5);
  if(LcR == L) return R;
  if(L>LcR) more=true;
  oldL=LcR;
  
  for( ; ; ){
    lR += (more)? dlR : -dlR;
    R=exp(lR);
    P= (*this)(R,0.,dPR,dPz);
    LcR=pow(R*R*R*dPR,0.5);
    if(LcR == L) return R;
    if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
	return R;}
    oldL=LcR;
  }
  
}
double IsoPot::LfromRc(const double R, double* dR) const
{
  double dPR,dPz,P;
  P = (*this)(R,0.,dPR,dPz);
  return sqrt(R*R*R*dPR);  
}
Frequencies IsoPot::KapNuOm(const double R) const
{
  //std::cout << "Lies from KapNuOm   ";
  Frequencies output= 100000.*R;
  return output;

}

ostream& operator<< (ostream& to, const IsoPot& P)
{
  to << "Isochrone potential with M=" << P.GM/(Units::G) << " b=" << P.b << '\n';
    return to;
}








