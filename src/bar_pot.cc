#include <math.h>
#include <stdio.h>
#include "bar_pot.h"
void bar_pot::reset_K(double fac){
	K*=fac;
}
void bar_pot::Sormani(void){
	double v0=220./978.,rq=1.5,AS=0.4,Kold=K;
	K=pow(v0,2)*pow(rq,3)*3*exp(2)/20*AS;
	printf("Sormani: Kold, Knew = %f %f\n",Kold,K);
}
double bar_pot::Phi(double R,double z){
	double ms=m2(R,z);
	return -R*R*pow(Rb2+ms,-2.5)*K;
}
double bar_pot::Phi(double R,double z,double &dPdR,double &dPdz){
	dPdR=dPhidR(R,z); dPdz=dPhidz(R,z);
	return Phi(R,z);
}
double bar_pot::dPhidR(double R,double z){
	double ms=m2(R,z);
	return -R*(2*Rb2-3*R*R+2*pow(z/q,2))*pow(Rb2+ms,-3.5)*K;
}
double bar_pot::oRddRRdPhidR(double R,double z){//(1/R)d(R dPhi/dR)
	double ms=m2(R,z);
	return -(2*(2*Rb2-6*R*R+2*pow(z/q,2))*(Rb2+ms)
		 -7*pow(R,2)*(2*Rb2-3*R*R+2*pow(z/q,2)))*pow(Rb2+ms,-4.5)*K;
}
double bar_pot::dPhidz(double R,double z){
	double ms=m2(R,z);
	return 5*R*R*z/pow(q,2)*pow(Rb2+ms,-3.5)*K;
}
double bar_pot::d2Phidz2(double R,double z){
	double ms=m2(R,z);
	return -R*R*pow(q,-2)*(-5*(Rb2+ms)+35*pow(z/q,2))*pow(Rb2+ms,-4.5)*K;
}
double bar_pot::rho(double R,double z){
	return (oRddRRdPhidR(R,z)+d2Phidz2(R,z)-4*Phi(R,z)/pow(R,2))/fourPiG;
/*	return -(rho0_fourPiG*((Rb2+ms)*(4*Rb2-12*R*R+(4*z*z-5*R*R)/pow(q,2))
		      -7*R*R*(2*Rb2-3*R*R+pow(z/q,2)*(2-5/pow(q,2))))*pow(Rb2+ms,-4.5)
			+4*Phi(R,z)/pow(R,2))/fourPiG; */
}
