/*******************************************************************************
*                                                                              *
* StackelPot.cc                                                                *
*                                                                              *
* C++ code written by Paul McMillan, 2007-                                     *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
#include "StackelPot.h"
#include <cmath>

//double StackelPotential::findG(const double tau)
//{
//  register double x=sqrt(( gamma + tau )/(-gamma)), ix=1./x;  

//  return TPi*(Units::G)*rh0*(-alpha)*ix*atan(x);
//}

//double StackelPotential::findG(const double tau, double& dG)
//{
//  register double x=sqrt((tau+gamma)/(-gamma)), ix=1./x;  
 
//dG =  TPi*(Units::G)*rh0*(-alpha)
//  *0.5*ix*ix*(ix*atan(x)/(gamma) + 1./tau);

//return TPi*(Units::G)*rh0*(-alpha)*ix*atan(x);
//}

double StackelPotential::operator() (const double R, const double z) const
{
  if(alpha==gamma)  {    
    double r=sqrt(R*R + z*z), ia= sqrt(-1./alpha);
    return -Gfac*atan(r*ia)/(r*ia);
  }
  double Rsq=R*R, zsq=z*z, A=-alpha-gamma+Rsq+zsq, // A=lamda+nu, B=lamda*nu 
    B=alpha*gamma-gamma*Rsq-alpha*zsq, temp=sqrt(A*A-4.*B),
    lamda=0.5*(A+temp), nu=0.5*(A-temp); // find de Zeeuw's lamda, nu from R,z

  double Glam,Gnu;
  if(R==0. || lamda <= -alpha){
    double e=sqrt(1-gamma/alpha);
    Glam = Gfac*sqrt(-gamma)*asin(e)/(sqrt(-alpha)*e);
  } else {
    double x=sqrt((gamma+lamda)/(-gamma)),  ix=1./x;
    Glam = Gfac*ix*atan(x);             // de Zeeuw's G(lamda)
  }

  if(z==0. || nu <= -gamma ){ Gnu=Gfac; } else {  
    double x=sqrt((gamma+nu)/(-gamma)), ix=1./x;
    Gnu = Gfac*ix*atan(x);         // G(nu)
  }

  return (-(lamda+gamma)*Glam+(nu+gamma)*Gnu)/(lamda-nu); // Phi
}

double StackelPotential::operator() (const double R, const double z, double& dPdR,
			       double& dPdz) const
{
  if(alpha == gamma)  {
    double r=sqrt(R*R + z*z), ia= sqrt(-1./alpha);
    double P= -Gfac*atan(r*ia)/(r*ia),
      dPdr= (-P/r) -TPi*(Units::G)*rh0*alpha*alpha/(r*(r*r - alpha));
    dPdR=dPdr*R/r;
    dPdz=dPdr*z/r;
    return P;
  }

  double Rsq=R*R, zsq=z*z, A=-alpha-gamma+Rsq+zsq, 
    B=alpha*gamma-gamma*Rsq-alpha*zsq, temp=sqrt(A*A-4.*B),
    lamda=0.5*(A+temp), nu=0.5*(A-temp), 
    dlamdR=R*(1.+(gamma-alpha+Rsq+zsq)/temp),
    dlamdz=z*(1.+(alpha-gamma+Rsq+zsq)/temp),
    dnudR= R*(1.-(gamma-alpha+Rsq+zsq)/temp),
    dnudz= z*(1.-(alpha-gamma+Rsq+zsq)/temp),
    ilammnu=1./(lamda-nu); // lamda, nu, and various derivatives

  double Glam,Gnu,dGlam,dGnu;
  if(R==0. || lamda <= -alpha){
    double e=sqrt(1-gamma/alpha), romesq = sqrt(gamma/alpha);
    Glam=Gfac*romesq*asin(e)/e;
    dGlam = Gfac*romesq*(romesq - asin(e)/e)/(2.*e*e);
  }else{
    double x=sqrt((gamma+lamda)/(-gamma)),  ix=1./x;
    Glam = Gfac*ix*atan(x); // de Zeeuw's G(lamda) & dG/dlamda
    dGlam = Gfac*0.5*ix*ix*(ix*atan(x)/(gamma) + 1./lamda);
  }

  if(z==0. || nu <= -gamma){
    Gnu=Gfac;
    dGnu=-Gfac*alpha/(3.*gamma);
  }else{
    double x=sqrt((gamma+nu)/(-gamma)), ix=1./x;
  Gnu = Gfac*ix*atan(x); // G(nu) and dG/dnu
  dGnu = Gfac*0.5*ix*ix*(ix*atan(x)/(gamma) + 1./nu);
  }

  double P= -((lamda+gamma)*Glam - (nu+gamma)*Gnu)*ilammnu;
  double dPdlam= -ilammnu*(Glam+(lamda+gamma)*dGlam+P),
    dPdnu= -ilammnu*(-Gnu-(nu+gamma)*dGnu-P);

    // Use chain rule to find derivs
  dPdR = dPdlam*dlamdR + dPdnu*dnudR;
  dPdz = dPdlam*dlamdz + dPdnu*dnudz;
  return P;
}

//double StackelPotential::operator() (const double R, 
//                           double& dPdR, double& d2PdRR) const
//{

//}

//double StackelPotential::operator() (const double R, const double z,
//			   double& dPdR, double& dPdz,
//			   double& d2PdRR, double& d2Pdzz, double& d2PdRz) const
//{

//}


// double StackelPotential::RfromLc(const double L, double* dR) const
// {
//   bool more=false;
//   double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL;
//   R=exp(lR);
//   P= (*this)(R,0.,dPR,dPz);
//   LcR=pow(R*R*R*dPR,0.5);
//   if(LcR == L) return R;
//   if(L>LcR) more=true;
//   oldL=LcR;
  
//   for( ; ; ){
//     lR += (more)? dlR : -dlR;
//     R=exp(lR);
//     P= (*this)(R,0.,dPR,dPz);
//     LcR=pow(R*R*R*dPR,0.5);
//     if(LcR == L) return R;
//     if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
// 	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
// 	return R;}
//     oldL=LcR;
//   }
  
// }

//double StackelPotential::LfromRc(const double R, double* dR) const
//{
//  double dPR,dPz,P;
//  P = (*this)(R,0.,dPR,dPz);
//  return sqrt(R*R*R*dPR);  
//}

Frequencies StackelPotential::KapNuOm(const double R) const
{
  std::cout << "Lies from KapNuOm   ";
  Frequencies output= 100000.*R;
  return output;

}

ostream& operator<< (ostream& to, const StackelPotential& P)
{
  to << "Oblate Stackel potential with rho_0=" << P.rh0
     << " a=" << sqrt(-(P.alpha)) << " c=" << sqrt(-(P.gamma)) << '\n';
    return to;
}








