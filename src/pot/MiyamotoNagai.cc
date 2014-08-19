/*******************************************************************************
*                                                                              *
* MiyamotoNagai.cc                                                             *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2006-07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
#include "MiyamotoNagai.h"
#include <cmath>

double MiyamotoNagai::operator() (const double R, const double z) const
{
  double AZB=A+sqrt(z*z+Bq);
  return -GM/sqrt(R*R+AZB*AZB);
}

double MiyamotoNagai::operator() (const double R, const double z, double& dPdR,
			       double& dPdz) const
{
  double  ZB=sqrt(z*z+Bq), AZB = A + ZB,
    F = 1./(R*R+AZB*AZB),rtF=sqrt(F);
  dPdR = GM*R*F*rtF;
  dPdz = ZB? GM*z*F*rtF*AZB/ZB : 0.;
  return -GM*rtF;
}



double MiyamotoNagai::RfromLc(const double L, double* dR) const
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
double MiyamotoNagai::LfromRc(const double R, double* dR) const
{
  double dPR,dPz,P;
  P = (*this)(R,0.,dPR,dPz);
  return sqrt(R*R*R*dPR);  
}
Frequencies MiyamotoNagai::KapNuOm(const double R) const
{
  double dPR,dPz,P;
  P = (*this)(R,0.,dPR,dPz);
  double F = 1./(R*R+ABq),rtF=sqrt(F),
    om2    = dPR/R, 
    nu2    = GM*ABoB*F*rtF,
    kappa2 = GM*(F*rtF - 3.*R*R*F*F*rtF) + 3.*om2;
  Frequencies output= sqrt(kappa2);
  output[1] = sqrt(nu2);
  output[2] = sqrt(om2);

  return output;
}

ostream& operator<< (ostream& to, const MiyamotoNagai& P)
{
  to << "Miyamoto-Nagai potential with mass=" << P.GM/Units::G
     << " a=" << P.A << " b=" << sqrt(P.Bq) << '\n';
    return to;
}








