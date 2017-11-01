/*******************************************************************************
*                                                                              *
* Toy_Isochrone.cc                                                             *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2007
*                     James Binney 2017    *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include "Toy_Isochrone.h"
#include "Numerics.h"
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// class ToyIsochrone


void ToyIsochrone::psifunc(const double p, double& f, double& df) const
{
    psi  = p;
    spsi = sin(psi);
    cpsi = cos(psi);
    f    = psi-eps*spsi-tr;
    df   = 1.-eps*cpsi;
}

void ToyIsochrone::psisolve() const
// solves  tr = psi - eps * sin(psi);    for psi.
{
    const double tol=1.e-8;
    psi = tr;
    if(fabs(eps)<tol || fabs(tr)<tol || fabs(tr-Pi)<tol || fabs(tr-TPi)<tol) {
        spsi = sin(psi);
        cpsi = cos(psi);
        return;
    }
    register double ps=psi;
    if(ps<Pi) ps = rtsafe(this,&ToyIsochrone::psifunc,ps,Pi,tol);
	else  ps = rtsafe(this,&ToyIsochrone::psifunc,Pi,ps,tol);
    if(ps!=psi) {
        psi  = ps;
        spsi = sin(psi);
        cpsi = cos(psi);
    }
}

double ToyIsochrone::catan(const double& aa) const
// Evaluates the function atan(aa tan(psi/2)).
// using tan(x/2) = sin(x)/[1+cos(x)]
{

// Three ranges.
    if(cpsi>=0. && spsi>=0.)                //      0 > psi > Pi/2
        return atan2(aa*spsi,1.+cpsi); 
    else if(cpsi<0. && spsi>=0)             //   Pi/2 > psi > Pi
	return Pih-atan2(1.+cpsi,spsi*aa);
    else if(cpsi<0.) {                      //     Pi > psi > 3Pi/2
	register double temp=atan2(1.+cpsi,spsi*aa);
	if(temp>0) return Pi3h-temp;
	return Pih-temp;
    } else                                  //  3Pi/2 > psi > 2Pi
        return Pi+atan2(aa*spsi,1.+cpsi); 
}

double ToyIsochrone::wfun(const double& fac) const
{
    if(fac==0.) return Pih; // radial orbit <=> e==1.
    else {
      if((1.-e)==0.) return catan(1.e99) + 
		       catan(sqrt((a*(1.+e)+2.)/(a*(1.-e)+2.)))*fac;
      else return   catan(sqrt((1.+e)/(1.-e)))
	       + catan(sqrt((a*(1.+e)+2.)/(a*(1.-e)+2.)))*fac;
    }
}

double ToyIsochrone::wfundw(const double& fac, double& dwdpsi) const
{
    if(fac==0.) { // radial orbit <=> e==1.
	dwdpsi = 0.;
	return Pih;
    } else {
        double a1   = sqrt( (1.+e) / (1.-e) );
        double a2   = sqrt( (a*(1.+e)+2.) / (a*(1.-e)+2.) );
        double cps1 = 1.+cpsi;
        double sps1 = 1.-cpsi;
        dwdpsi      = 1./(cps1/a1+sps1*a1) + fac/(cps1/a2+sps1*a2);
        return        catan(a1) + catan(a2)*fac;
    }
}

double ToyIsochrone::Fint(double p, double q, double Pc, double Qc) const
// Evaluates int_0^psi (Pc+Qc cos[t])/ (p+q cos[t])^2  dt.
{
    double subint = 2./sqrt(p*p-q*q)*catan(sqrt((p-q)/(p+q)));
    return ((p*Qc-q*Pc)*spsi/(p+q*cpsi) + (p*Pc-q*Qc)*subint)/(p*p-q*q);
}

////////////////////////////////////////////////////////////////////////////////

vec4 ToyIsochrone::lower_bounds(const double x, const double v) const
{
    return -upper_bounds(x,v);
}

vec4 ToyIsochrone::upper_bounds(const double x, const double v) const
{
    double a[4] = {sqrt(1.e5*x)*v, sqrt(1.e3*x), 1.e2*x*v, 1.*x};
    return vec4((const double*)a); 
}

void ToyIsochrone::Set()
{
    M = gamma*gamma;
    b = beta*beta;
    //if(b==0.) TorusError("ToyIsochrone: scale radius = 0 not allowed",2);
    //if(M==0.) TorusError("ToyIsochrone: mass = 0 not allowed",2);
    if(b>0. && M>0.) {
	sMb  = sqrt(M*b);
	sMob = sMb/b;
	jp   = fabs(Lz)/sMb;
    }
    derivs_ok = false;
}

void ToyIsochrone::dump() const
{
    cerr<< "\ntest output from Iso:"
        << "\njr,jt,tr,tt     = " << jr <<','<<   jt<<','<<   tr<<','<< tt
        << "\nr,th,p_r,p_th   = " <<  r <<','<<   th<<','<<   pr<<','<< pt
        << "\na,e,eps,u       = " <<  a <<','<<    e<<','<<  eps<<','<< u
	<< "\nH,wr,wt,at      = " <<  H <<','<<   wr<<','<<   wt<<','<< at
	<< "\nwh,sq,sGam,chi  = " << wh <<','<<   sq<<','<< sGam<<','<< chi
	<< "\npsi,spsi,cpsi   = " << psi<<','<< spsi<<','<< cpsi<<'\n'; 
}

ToyIsochrone::ToyIsochrone()
{  
    IsoPar p=1.f;
    set_parameters(p);
    Set();
}

ToyIsochrone::ToyIsochrone(const IsoPar& p)
{
    set_parameters(p);
    Set();
}

ToyIsochrone::ToyIsochrone(const ToyIsochrone& I): gamma(I.gamma),beta(I.beta),Lz(I.Lz),r0(I.r0)
{ 
    Set();
}
inline double ToyIsochrone::GM(double r,double F,double bi){
	//returns GM that gives dP/dr=F at r
	double s=1+sqrt(1+pow(r/bi,2));
	return s*pow(s*(s-1)*bi,2)*F/(s-2);
}
ToyIsochrone::ToyIsochrone(double r1,double F1,double r2,double F2,double Jphi){
	//Returns isochrone that gives dPdr Fi at ri
	if(r1>r2){//make sure r1 is the smaller radius
		double t=r1,f=F1; r1=r2; F1=F2; r2=t; F2=f;
	}
	double bl=.1,bh=10,eps=1.e-4;
	double dl=GM(r2,F2,bl)-GM(r1,F1,bl),dh=GM(r2,F2,bh)-GM(r1,F1,bh);
	while(dl<0){
		bl*=.5; dl=GM(r2,F2,bl)-GM(r1,F1,bl);
	}
	while(dh>0){
		bh*=1.5; dh=GM(r2,F2,bh)-GM(r1,F1,bh);
	}//now we have bracketed b where masses agree
	while(bh-bl>eps){
		double b=sqrt(bl*bh), d=GM(r2,F2,b)-GM(r1,F1,b);
		if(d<0) bh=b; else bl=b;
	}
	beta=sqrt(.5*(bl+bh));
	gamma=sqrt(GM(r2,F2,beta*beta));
	Lz=Jphi;
	r0=0;
	Set();
}

PSPD ToyIsochrone::Forward(const PSPD& JT) const
// Transforms action-angle variables (JT) to phasespace co-ordinates (r,theta)
// and their conjugated momenta; theta is latitude rather than polar angle.
{
	derivs_ok = true;
	register double e2;
	double fac;
// Extract and scale the actions and angles.
	jr = double(JT(0)) / sMb;
	jt = double(JT(1)) / sMb;
	tr = double(JT(2));
	tt = double(JT(3));
// Make sure that tr is in the correct range
	if(std::isnan(tr) || std::isinf(tr) || fabs(tr)>INT_MAX) 
		tr = 0.; // just in case  
	while(tr<0.)  tr+=TPi;
	while(tr>TPi) tr-=TPi;
	if(jr<0. || jt<0.|| jp<0.) {
		TorusError("ToyIsochrone: negative action(s)",2); 
		return PSPD(0.);
	}
// Set up auxilliary variables independent of the angles tr and tt.
	at  = jp+jt;
	sGam= (jp==0.) ? 1. : sqrt(1.-pow(jp/at,2));
	sq  = hypot(2.,at);
	fac = at/sq;
	HH  = 2./(2.*jr+at+sq);
	H2  = HH*HH;
	H   =-0.5*H2;
	a   = 1./H2-1.;
	if(at==0.) e = 1.;
	else {
		e2    = 1. - pow(at/a,2)/H2;
		e     = (e2>0) ? sqrt(e2) : 0;
		if(e2<-1.e-3) TorusError("ToyIsochrone: bad e in Forward",2);
	}
	ae  = a*e;
	eps = ae*H2;
	wr  = H2*HH;
	wt0r= 0.5*(1.+fac);
	wt  = wt0r*wr;
    //cerr << jt << " " << wt << " " << jt*wt<< "\n";
// Solve for the relevant other variables.
	psisolve(); 
	wh  = wfun(fac);
	chi = tt-wt0r*tr+wh;
	u   = a*(1.-e*cpsi);
// Calculate co-ordinates and the conjugate momenta.
	r   = (u==0.)?  0. : sqrt(u*(u+2.));
	th  = asin(sGam*sin(chi));
	pt  = (at==0.)? 0. : at*sGam*cos(chi)/cos(th);
	if(r==0.) pr = sqrt(1.-H2);
	else  pr = HH*ae*spsi/r;
// re-scale and return
	return PSPD(fmax(b*r+r0,b*1.e-6), th, pr*sMob, pt*sMb); 
}

////////////////////////////////////////////////////////////////////////////////
PSPT ToyIsochrone::Forward3D(const PSPT& JT3) const
{ // For the 3D case
	PSPD QP2, JT2 = JT3.Give_PSPD();
	PSPT QP3;

	QP2 = Forward(JT2); 
	QP3.Take_PSPD(QP2);  // do the 2D part

	QP3[5] = JT3(2); // old:JT3(2)/(QP3(0)*cos(QP3(1))); // Angular velocity/mom
	if(JT3(1) == 0.) {
		QP3[2] = (QP3(5)>0.)? JT3(5)+wh-wt0r*JT3(3) : JT3(5)-wh+wt0r*JT3(3);
		if(std::isnan(QP3[2]) || std::isinf(QP3[2]) || fabs(QP3[2])>INT_MAX) 
			QP3[2] = 0.; // just in case  
		while(QP3(2)<0. ) QP3[2] += TPi; 
		while(QP3(2)>TPi) QP3[2] -= TPi;
		return QP3;
	}

	double sinu = (JT3(2)>0.)? tan(QP3(1))*Lz/sqrt(JT3(1)*(JT3(1)+2.*fabs(Lz))) :
		      -tan(QP3(1))*Lz/sqrt(JT3(1)*(JT3(1)+2.*fabs(Lz)));
	ufn = (sinu>1.)? Pih : (sinu<-1.)? -Pih : asin(sinu);
	if(QP3(4)<0.) ufn= Pi-ufn; 
	QP3[2] = JT3(5)+ufn;

	QP3[2] -= (JT3(2)>0.)? JT3(4) : -JT3(4); // Note opposite sign to Backward
	if(std::isnan(QP3[2]) || std::isinf(QP3[2]) || fabs(QP3[2])>INT_MAX) 
		QP3[2] = 0.; // just in case  
	while(QP3(2)<0. ) QP3[2] += TPi; 
	while(QP3(2)>TPi) QP3[2] -= TPi; 

	return QP3;
}

PSPD ToyIsochrone::ForwardWithDerivs(const PSPD& JT, double dQdT[2][2]) const
// Transforms action-angle variables (JT) to phasespace co-ordinates (r,theta)
// and their conjugated momenta, the derivs of (r,theta) w.r.t. the angle vars
// are also returned; theta is latitude rather than polar angle.
{
	PSPD fail(0.);
#include "toy_isochrone_common2.cpp"
	return PSPD(b*r+r0, th, pr*sMob, pt*sMb);
}

PSPD ToyIsochrone::ForwardWithDerivs(const PSPD& JT, double dQdT[2][2], 
			       double dPdT[2][2] ) const
// Transforms action-angle variables (JT) to phasespace co-ordinates (r,theta)
// and their conjugated momenta, the derivs of (r,theta,pr,ptheta) w.r.t. the 
//angle vars are also returned; theta is latitude rather than polar angle.
{
	PSPD fail(0.);
#include "toy_isochrone_common2.cpp"

    dPdT[0][0] = HH*ae/r*((-1.+1./(1.-eps*cpsi))/eps-spsi/r*dQdT[0][0]/b)*sMob;
    dPdT[0][1] = 0.;
    dPdT[1][1] = (at==0) ? 0. : at*sGam/csth*
                                (sGam*schi/csth*cchi*dQdT[1][1]-schi)*sMb;
    dPdT[1][0] = (at==0) ? 0. : at*sGam/csth*(sGam*schi/csth*cchi*dQdT[1][0] 
					       - schi*dchidtr)*sMb;
    return PSPD(b*r+r0, th, pr*sMob, pt*sMb);
}


void ToyIsochrone::Derivatives(double dQPdJ[4][2]) const
// Calculates the derivatives of the spherical phasespace co-ordinates w.r.t.
// Jr,Jt. A call of one of Forward() or ForwardWithDerivs() has to be preceeded.
{
#include "toy_isochrone_common.cpp"
}

void ToyIsochrone::Derivatives(double dQPdJ[4][2], Pdble dQPdA[4]) const
// Calculates the derivatives of the spherical phasespace co-ordinates w.r.t.
// Jr,Jt and gamma,beta,Lz,r0. A prior call of one of Forward() or ForwardWithDerivs()
// is required.
{
#include "toy_isochrone_common.cpp"
// w.r.t. Jp = fabs(Lz):
	if(jt==0.) {
		dQPdA[1][2] = 0.;
		dQPdA[3][2] = 100.;		// actually, it's infinite
	} else {
		double xGam= (cGam-1.)/(sGam*at);
		dQPdA[1][2] = (sGam*cchi*chix+schi*cGam*xGam)/csth;
		dQPdA[3][2] = pt/at + at *( cchi*cGam/csth*xGam + sGam *
					    (-schi/csth*chix +cchi*sith/pow(csth,2)*dQPdA[1][2]));
	}
    dQPdA[0][2] = dQPdJ[0][1];
    dQPdA[2][2] = dQPdJ[2][1];
    dQPdA[1][2]/= sMb;
// w.r.t. M and b:
    double temp = 0.5*(dQPdJ[0][0]*jr+dQPdJ[0][1]*jt+dQPdA[0][2]*jp);
    dQPdA[0][0] =-temp / sMob;
    dQPdA[0][1] = r - temp * sMob;
    temp = 0.5*(dQPdJ[1][0]*jr+dQPdJ[1][1]*jt+dQPdA[1][2]*jp);
    dQPdA[1][0] =-temp / sMob;
    dQPdA[1][1] =-temp * sMob;
    temp = 0.5*(dQPdJ[2][0]*jr+dQPdJ[2][1]*jt+dQPdA[2][2]*jp);
    dQPdA[2][0] = (0.5*pr/b - temp) / sMob;
    dQPdA[2][1] =-(0.5*pr/b + temp) * sMob;
    temp = 0.5*(dQPdJ[3][0]*jr+dQPdJ[3][1]*jt+dQPdA[3][2]*jp);
    dQPdA[3][0] = (0.5*pt - temp) / sMob;
    dQPdA[3][1] = (0.5*pt - temp) * sMob;
// w.r.t. r0:
    dQPdA[0][3] = 1.;
    dQPdA[1][3] = 0.;
    dQPdA[2][3] = 0.;
    dQPdA[3][3] = 0.;
// transfer to gamma, beta, and Lz
    register double tgamma=2*gamma, tbeta=2*beta;
    for(register short i=0; i<4; i++) {
	dQPdA[i][0] *= tgamma;
	dQPdA[i][1] *= tbeta;
	if(Lz < 0.) dQPdA[i][2] = -dQPdA[i][2];
    }
}

PSPD ToyIsochrone::Backward(const PSPD& QP) const
// Transforms phasespace co-ordinates (r,theta) and their conjugated momenta
// (QP) to action-angle variables (JT);
// theta is latitude rather than polar angle.
{
    derivs_ok = true;
    register double e2,csth;
     double fac;
// extract and scale co-ordinates
    r   = (QP(0)-r0) / b;
    th  = QP(1);
    pr  = QP(2) / sMob;      // pr  = dr/dt        => [pr]  = Sqrt[GM/b]
    pt  = QP(3) / sMb;       // pth = r^2 * dth/dt => [pth] = Sqrt[GMb]
    csth= cos(th);
// Perform some consistency checks
    if(csth==0.&&jp!=0.) {
	TorusError("ToyIsochrone: Jp!=0 && R=0 in Backward",2);
	return PSPD(0.);
    }
    if(r==0. && pt!=0.) {
	TorusError("ToyIsochrone: pt!=0 && r=0 in Backward",2);
	return PSPD(0.);
    }
// Set up auxialiary variables
    at  = (jp==0.)?  fabs(pt) : hypot(pt,jp/csth);
    sq  = hypot(2.,at);
    wt0r= 0.5*(1.+at/sq);
    sGam= (jp==0.)? 1. : sqrt(1.-pow(jp/at,2));
    u   = hypot(r,1.)-1.;
    H2  = (pt==0.)? -pr*pr : -pr*pr-pow(pt/r,2);
    if(jp!=0.) H2 -= pow(jp/(r*csth),2);
    H2 += 2./(u+2.);
    if(H2<=0.) {
	TorusError("ToyIsochrone: H>=0 in Backward",2);
	return PSPD(0.);
    }
    H   =-0.5*H2;
    HH  = sqrt(H2);
    a   = 1./H2-1.;
    if(at==0.) e = 1.;
    else {
        e2 = 1. - pow(at/a,2)/H2;
        e  = (e2>0) ? sqrt(e2) : 0;
    }
    ae  = a*e;
    eps = ae*H2;
    fac = at/sq;
    wr  = H2*HH;
    wt  = wt0r*wr;
    //psisolve();
    //wh  = wfun(fac);
    jt  = at-jp;
    if(e==0.) { // circular orbit
      wh  = wfun(fac);
      psi = wh/wt0r;
	cpsi= cos(psi);
        spsi= sin(psi);
        jr  = 0.;
        tr  = psi;
        if(sGam==0.)      chi = 0.;
        else if(sGam==1.) chi = (pt<0) ? Pi-fabs(th) : fabs(th);
        else              chi = acos(fmax(-1.,fmin(1.,pt*csth/at/sGam)));
        if(th<0.) chi=TPi-chi;
        tt  = chi;
    } else {
	cpsi = fmax(-1.,fmin(1.,(1.-u/a)/e));
        if(cpsi==1.)       psi = 0.;
        else if(cpsi==-1.) psi = Pi;
        else               psi = (pr<0.)? TPi-acos(cpsi) : acos(cpsi);
        spsi= sin(psi);
        jr  = fmax(0., 1./HH-0.5*(at+sq));
        tr  = psi-eps*spsi;
        if(sGam==0.)      chi = 0.;
        else if(sGam==1.) chi = (pt<0) ? Pi-fabs(th) : fabs(th);
        else              chi = acos(fmax(-1.,fmin(1.,pt*csth/at/sGam)));
        if(th<0.) chi=TPi-chi;
	wh  = wfun(fac);
        tt  = wt0r*tr+chi-wh;
    }
    if(std::isnan(tt) || std::isinf(tt) || fabs(tt)>INT_MAX) 
      tt = 0.; // just in case  
    while(tt<0.)  tt+=TPi;
    while(tt>TPi) tt-=TPi;
// re-scale and return
    return PSPD(jr*sMb, jt*sMb, tr, tt);
}

////////////////////////////////////////////////////////////////////////////////
PSPT ToyIsochrone::Backward3D(const PSPT& QP3) const
{ // (r,vartheta,phi,pr,pvartheta,Jphi) -> (Jr,J_vartheta,Jphi,theta_r,theta_vartheta,theta_phi)
  PSPD Jt2, QP2 = QP3.Give_PSPD();
  PSPT Jt3;

  Jt2 = Backward(QP2);
  Jt3.Take_PSPD(Jt2); // do the 2D part

  Jt3[2] = QP3(5); // old: QP3(0)*cos(QP3(1))*QP3(5); // Angular momentum. 
  if(QP3(1) == 0. && QP3(4) == 0.) {//planar case
    // This should be 
    // Jt3(5) = tt  = wt0r*tr+chi-wh; so wt0r = omega_theta/omega_r; 
    // chi= angle in some plane;  wh = mess of arctans.
    Jt3[5] = (QP3(5)>0.)? QP3(2)-wh+wt0r*Jt3(3) : QP3(2)+wh-wt0r*Jt3(3);
    AlignAngles3D(Jt3);
    return Jt3;
  }//Now non-planar case
  double sinu = (Jt3(2)>0.)? tan(QP3(1))*Lz/sqrt(Jt3(1)*(Jt3(1)+2.*fabs(Lz))) :
    -tan(QP3(1))*Lz/sqrt(Jt3(1)*(Jt3(1)+2.*fabs(Lz)));
  ufn = (sinu>1.)? Pih : (sinu<-1.)? -Pih : asin(sinu);
  if(QP3(4)<0.) ufn= Pi-ufn; 
  Jt3[5] = QP3(2)-ufn;
  Jt3[5] += (Jt3(2)>0.)? Jt3(4) : -Jt3(4);
  AlignAngles3D(Jt3);
  return Jt3;
}


PSPT ToyIsochrone::Forward3DwithDerivs(const PSPT& JT3, double dQdT3[3][3]) const
// Transforms action-angle variables (JT) to phasespace co-ordinates (r,theta)
// and their conjugated momenta, the derivs of (r,theta,phi) w.r.t. the angle vars
// are also returned; theta is latitude rather than polar angle.
{
	PSPT QP3=Forward3D(JT3);
	PSPD JT=JT3.Give_PSPD();
	double dQdT[4][2];
	PSPT fail(0.);
#include "toy_isochrone_common2.cpp" //fills up dQdT
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++) dQdT3[i][j]=dQdT[i][j];
	double incl=acos(jp/at), coti=1/tan(incl);
	double tanth=tan(th),secth2=1+pow(tanth,2),secu=1/cos(ufn);
	double fac1=secu*coti*secth2;
	dQdT3[2][0]=fac1*dQdT[1][0];
	dQdT3[2][1]=fac1*dQdT[1][1];
	dQdT3[2][2]=1;
	if(jp>0) dQdT3[2][1]-=1; else dQdT3[2][1]+=1;
	dQdT3[0][2]=0; dQdT3[1][2]=0;
	return QP3;
}

void ToyIsochrone::Derivatives3D(double dQdJ[3][3]) const
// Calculates the derivatives of the spherical  co-ordinates w.r.t.
// Jr,Jt,Jp. A previous call of one of Forward() or ForwardWithDerivs()
// is required.
{
	double dQPdJ[4][2];
#include "toy_isochrone_common.cpp"//fill up dQPdJ
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++) dQdJ[i][j]=dQPdJ[i][j];
	double incl=acos(jp/at), coti=1/tan(incl), sin3i=pow(sin(incl),3);
	double tanth=tan(th),secth2=1+pow(tanth,2),secu=1/cos(ufn);
	dQdJ[2][0]=secu* coti*secth2*dQdJ[1][0];
	dQdJ[2][1]=secu*(coti*secth2*dQdJ[1][1]-tanth/sin3i*jp/pow(at,2)/sMb);
	dQdJ[0][2]= dQdJ[0][1];
	if(jt==0.) {
		dQdJ[1][2] = 0.;
	} else {
		double xGam= (cGam-1.)/(sGam*at);
		dQdJ[1][2] = (sGam*cchi*chix+schi*cGam*xGam)/csth;
		//printf("jt=%f %f\n",at,jp+jt);
	}
	dQdJ[1][2]/=sMb;
	dQdJ[2][2]=secu*(coti*secth2*dQdJ[1][2]+tanth/sin3i*jt/pow(at,2)/sMb);
}

///end of Isochrone.cc//////////////////////////////////////////////////////////
