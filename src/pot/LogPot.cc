/*******************************************************************************     
*                                                                              *
* LogPot.cc                                                                    *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2006-07                                   *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#include <iomanip>
#include "LogPot.h"
#include <cmath>


void  LogPotential::error(const char* msgs) const
{
    cerr << " Error in class LogPotential: " << msgs << '\n';
    exit(1);
}

double LogPotential::operator() (const double R, const double z) const
{
    if(rc2==0. && R==0. && z==0.) error(" (R,z)=0 at zero core radius"); 
   
    return v0sqhalf * log(R*R+z*z*q2i+rc2) + plusconst;
}

double LogPotential::operator() (const double R, const double z, double& dPdR,
			   double& dPdz) const
{
    if(rc2==0. && R==0. && z==0.) error(" (R,z)=0 at zero core radius"); 
    register double Rq= R*R, zq= z*z,
		    mq= Rq + zq * q2i + rc2;
    
    dPdR = v0sq / mq;
    dPdz = dPdR * z * q2i;
    dPdR*= R;
    
    return v0sqhalf * log(mq);
}

double LogPotential::operator() (const double R, double& dPdR, double& d2PdRR) const
{
    if(rc2==0. && R==0.) error(" (R,z)=0 at zero core radius"); 
    register double Rq= R*R,
		    mq= Rq + rc2;
    
    dPdR   = v0sq / mq;
    d2PdRR = dPdR * (1.-2*Rq/mq);
    dPdR  *= R;
    
    return v0sqhalf * log(mq);
}

double LogPotential::operator() (const double R, const double z,
			   double& dPdR, double& dPdz,
			   double& d2PdRR, double& d2Pdzz, double& d2PdRz) const
{
    if(rc2==0. && R==0. && z==0.) error(" (R,z)=0 at zero core radius"); 
    register double Rq = R*R,
		    zq = z*z,
		    mq = Rq + zq * q2i + rc2;
    
    dPdR   = v0sq / mq;
    dPdz   = dPdR * q2i;
    d2PdRR = dPdR * (1.-2*Rq/mq);
    d2Pdzz = dPdz * (1. - 2*zq*q2i/mq);
    d2PdRz =-2 * dPdz * R * z / mq;
    dPdR  *= R;
    dPdz  *= z;
    
    return v0sqhalf * log(mq);
}


//double LogPotential::LfromRc(const double R, double* dR) const
//{
//  double dPR,dPz,P;
//  P = (*this)(R,0.,dPR,dPz);
//  return sqrt(R*R*R*dPR);  
//}



// double LogPotential::RfromLc(const double L_in, double* dR) const
// {
//   bool more=false;
//   double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL,L=fabs(L_in);
//   R=exp(lR);
//   P= (*this)(R,0.,dPR,dPz);
//   LcR=sqrt(R*R*R*dPR);
//   if(LcR == L) return R;
//   if(L>LcR) more=true;
//   oldL=LcR;
  
//   for( ; ; ) {
//     lR += (more)? dlR : -dlR;
//     R=exp(lR);
//     P= (*this)(R,0.,dPR,dPz);
//     LcR=sqrt(R*R*R*dPR);
//     if(LcR == L) return R;
//     if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
// 	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
// 	return R;}
//     oldL=LcR;
//   }
  
// }


Frequencies LogPotential::KapNuOm(const double R) const {

  //if(Rei) cerr << "DANGER: epicycle frequencies will be wrong\n";
  Frequencies epi;
  double dPR,dPz,d2PdRR,d2Pdzz,d2PdRz, tmp;
  (*this)(R,0.,dPR,dPz,d2PdRR,d2Pdzz,d2PdRz);
  tmp = dPR/R;
  epi[2] = sqrt(tmp);
  epi[1] = sqrt(d2Pdzz);
  epi[0] = sqrt(d2PdRR+3*tmp);
  return epi;
}



ostream& operator<< (ostream& to, const LogPotential& P)
{
                 to << "logarithmic Potential Phi = ";
    if(P.v0!=1.) to << P.v0 << "^2";
	else     to << '1';
		 to << "/2 ";
    if(P.q==1.)  to << "ln[R^2 + z^2";
	else     to << "ln[R^2 + (z/" << P.q << ")^2";
    
    if(P.rc!=0.) to << " + " << P.rc << "^2 ";
                 to << ']';
    return to;
}
