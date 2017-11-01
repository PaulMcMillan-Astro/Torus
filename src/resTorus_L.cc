/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
// Functionality for Nr,0,Nphi resonantly trapped tori

#include "resTorus_L.h"

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SGN(A)  ((A)>0? 1:-1)
#define SMALL 1e-7

static double PI=acos(-1),TPI=2*PI;

/*
Given an estimate of the actions of the resonant torus, the following creator refines that estimate,
creates the corresponding eTorus, finds the required derivatives of the Fourier terms and computes
the grid of tori needed to produce trapped tori
 */
resTorus_L::resTorus_L(Torus **Tgridin,int nrin,Actions Jresin,Potential *Phi,bar_pot *bar,double Omp,int3 resNin,double tolJ){
	tmax=10;
	Tgrid=Tgridin; nr=nrin; resN=resNin; resJp=prime_it(Jresin,0);
	findJp(Jresin,Phi,bar,Omp,tolJ); resJp=prime_it(Jresin,0);
	int ok=eT.AutoFit(Jresin,Phi,bar,Omp,tolJ);
	nangle=200;//nangle must be even
	t1grd=new double[nangle]; trgrd=new double[nangle]; Dgrd=new double[nangle];
	hres[0]=get_driver(eT);	for(int i=1;i<6;i++) hres[i]=0;
	double m2=get_bar(eT);
	printf("driver, bar: %g %g\n",hres[0],m2);
	Actions J=Jresin;
	double dJ=.05*J[0];
	J[0]+=dJ; J[2]+=dJ*resN[2]/resN[0];//holding J3' const, not quite tracking H0=const
	Torus T1; T1.AutoFit(J,Phi,tolJ);
	Frequencies Om=T1.omega(); Om[2]-=Omp;
	double resOm=resCond(Om);
	dJ*=2; J[0]-=dJ; J[2]-=dJ*resN[2]/resN[0];
	Torus T0; T0.AutoFit(J,Phi,tolJ);
	Om=T0.omega(); Om[2]-=Omp; resOm-=resCond(Om);
	m2*=2/Om(2);//only a rough estimate of non-res change to Jphi
	G=resOm/dJ;
	if(G*hres[0]<0) t1c=0;
	else t1c=PI;
	off=0;
	//printf("G %f\n",G);
	dJ*=.5; J[0]+=dJ; J[2]+=dJ*resN[2]/resN[0];
	getImin_max();
	if(G<0)	setI(.99*Imin);
	else if(Imax>0) setI(.99*Imax); else setI(1.001*Imax);
	Actions J0,J1; getJs(J0,J1,.3);
	printf("getJs (%f %f) (%f %f)\n",J0[0],J0[2],J1[0],J1[2]);
	fix_derivs(J,J0,J1,Phi,bar,Omp,tolJ);
	getImin_max();
	double safety=5;
	dJ=safety*fabs(m2);
	if(G<0)	setI(.99*Imin);
	else if(Imax>0) setI(.99*Imax); else setI(1.001*Imax);
	getJs(J0,J1,1);
	printf("getJs, dJ (%f %f) (%f %f) %f\n",J0[0],J0[2],J1[0],J1[2],dJ);
	J0=prime_it(J0,0); J1=prime_it(J1,0); //restore prime
/*	FILE *dfile=fopen("resTorus_L.out","w");
	fprintf(dfile,"%f %g %g %g %g\n",resJp[0],G,hres[0],hres[1],hres[2]);
	int ng=13; double Hb;
	for(int i=0;i<ng;i++){
		double fac=i/(double)(ng-1);
		Actions Jn=(1-fac)*J0+fac*J1;
		double jr=(1-fac)*sqrt(J0[0])+fac*sqrt(J1[0]); Jn[0]=jr*jr;
		printf("Jn: %f\n",Jn[0]);
		Jn=prime_it(Jn,1);//remove prime
		//Actions Jn=J0+i/(double)(ng-1)*(J1-J0);
		eTorus eTT(Jn,Phi,bar,Omp,tolJ);
		if(i==0) Hb=get_H(eTT);
		fprintf(dfile,"%f %g %g\n",prime_it(Jn,0)[0],get_H(eTT)-Hb,get_driver(eTT));
	}
	fclose(dfile);*/
	J0[0]-=dJ; J1[0]+=dJ;
	J0[0]=MAX(J0[0],.0001);
	dJ_grid[0]=(sqrt(J1[0])-sqrt(J0[0]))/(double)(nr-1);
	dJ_grid[1]=0; dJ_grid[2]=2*dJ;
	Jbar_grid[0]=sqrt(J0[0])+.5*dJ_grid[0];
	Jbar_grid[1]=J0[1]; Jbar_grid[2]=J0[2];
	printf("Grid Js: (%f %f) (%f %f)\n",J0(0),J0(2),J1(0),J1(2));
	vec4 tm;
	for(int k=0;k<2;k++){
		if(k==0){
			J0[2]-=dJ; J1[2]-=dJ;
		}else{
			J0[2]+=2*dJ; J1[2]+=2*dJ;
		}
		for(int i=0;i<nr;i++){
			double fac=i/(double)(nr-1);
			Actions Jn=(1-fac)*J0+fac*J1;
			double jr=(1-fac)*sqrt(J0[0])+fac*sqrt(J1[0]); Jn[0]=jr*jr;
			Jn=prime_it(Jn,1);//remove prime
			if(i>100){//to implement make if(i>0)
				Tgrid[i][k].FitWithFixToyPot(Jn,tm,Phi,tolJ);
			}else{
				Tgrid[i][k].AutoFit(Jn,Phi,tolJ);
				tm=Tgrid[i][k].TP();
			}
		}
	}
}
Actions resTorus_L::prime_it(Actions Jin,bool back){//If back==0 go J -> Jp, otherwise Jp -> J
	Actions Jout;
	if(back){
		Jout[0]=resN[0]*Jin[0];
		Jout[1]=resN[1]*Jin[0]+Jin[1];
		Jout[2]=resN[2]*Jin[0]+Jin[2];
	} else {
		Jout[0]=Jin[0]/(double)resN[0];
		Jout[1]=Jin[1]-Jin[0]*resN[1]/(double)resN[0];
		Jout[2]=Jin[2]-Jin[0]*resN[2]/(double)resN[0];
	}
	return Jout;
}
Matrix<double,3,3> resTorus_L::unprime(Matrix<double,3,3> &dJpdt){
	Matrix<double,3,3> dJdt;
	for(int i=0;i<3;i++){//convert to unprimed actions
		dJdt[0][i]=resN[0]*dJpdt[0][i];
		dJdt[1][i]=resN[1]*dJpdt[0][i]+dJpdt[1][i];
		dJdt[2][i]=resN[2]*dJpdt[0][i]+dJpdt[2][i];
	}
	return dJdt;
}
Vector<double,3> resTorus_L::prime_deriv(Vector<double,3> &Din){
	Vector<double,3> Dout;
	Dout[0]=Din[0]/(double)resN[0];
	Dout[1]=Din[1]-resN[1]*Din[0]/(double)resN[0];
	Dout[2]=Din[2]-resN[2]*Din[0]/(double)resN[0];
	return Dout;
}
Matrix<double,3,3> resTorus_L::prime_deriv(Matrix<double,3,3> &Din){
	Matrix<double,3,3> Dout;
	for(int i=0;i<3;i++){
		Dout[i][0]=Din[i][0]/(double)resN[0];
		Dout[i][1]=Din[i][1]-resN[1]*Din[i][0]/(double)resN[0];
		Dout[i][2]=Din[i][2]-resN[2]*Din[i][0]/(double)resN[0];
	}
	return Dout;
}
Matrix<double,3,3> resTorus_L::unprime_deriv(Matrix<double,3,3> &Din){
	Matrix<double,3,3> Dout;
	for(int i=0;i<3;i++){
		Dout[i][0]=resN[0]*Din[i][0];
		Dout[i][1]=resN[1]*Din[i][0]+Din[i][1];
		Dout[i][2]=resN[2]*Din[i][0]+Din[i][2];
	}
	return Dout;
}
Actions resTorus_L::changeJr(Actions J,double dJ){
	Actions Jout;
	Jout[0]=J[0]+dJ;
	for(int i=1;i<3;i++) Jout[i]=J[i]+dJ*resN[i]/resN[0];
	return Jout;
}
void resTorus_L::findJp(Actions &J,Potential *Phi,bar_pot *bar,double Omp,double tolJ){
	Torus T; T.AutoFit(J,Phi,tolJ);
	Frequencies Om=T.omega(); Om[2]-=Omp;
	double dJ=.05*J[2], Jmax=0, Jmin=0, resOm=resCond(Om);
	while(Jmax*Jmin==0){
		if(resN[2]*resOm>0){
			Jmin=J[2]; J[2]+=dJ;
		}else{
			Jmax=J[2]; J[2]-=dJ;
		}
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); Om[2]-=Omp; resOm=resCond(Om);
	}
	while(fabs(resOm)>10*SMALL && Jmax-Jmin>10*SMALL){
		J[2]=.5*(Jmax+Jmin);
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); Om[2]-=Omp; resOm=resCond(Om);
		if(resOm>0) Jmin=J[2];
		else Jmax=J[2];
	}
}
double resTorus_L::get_driver(eTorus &eT){
	for(int i=0;i<tmax;i++){
		if(eT.i1()(i)==resN(0) && eT.i2()(i)==resN(1) && eT.i3()(i)==resN(2)) return eT.hn()(i);
		if(eT.i1()(i)==-resN(0) && eT.i2()(i)==-resN(1) && eT.i3()(i)==-resN(2)) return eT.hn()(i);
	}
}
double resTorus_L::get_H(eTorus &eT){
	for(int i=0;i<tmax;i++)
		if(eT.i1()(i)==0 && eT.i2()(i)==0 && eT.i3()(i)==0) return eT.hn()(i);
}
double resTorus_L::get_bar(eTorus &eT){
	for(int i=0;i<tmax;i++)
		if(eT.i1()(i)==0 && eT.i2()(i)==0 && eT.i3()(i)==2) return eT.hn()(i);
}
bool resTorus_L::fix_derivs(Actions &Jres,Actions &J0,Actions &J1,
			    Potential *Phi,bar_pot *bar,double Omp,double tolJ){
	int ng=13;
	double Jrgrd[ng],Omgrd[ng],resOmgrd[ng],Hgrd[ng],hNgrd[ng],y2[ng];
	for(int i=0;i<ng;i++){
		double f=i/(float)(ng-1);
		Actions J=(1-f)*J0+f*J1;
		Jrgrd[i]=J[0];
		eTorus eTT(J,Phi,bar,Omp,tolJ);
		Frequencies Om=eTT.omega();
		Omgrd[i]=Om(0); resOmgrd[i]=resCond(Om);
		Hgrd[i]=get_H(eTT); hNgrd[i]=get_driver(eTT);
	}
	double Jrres=Jres[0];
	Frequencies Omres=Omgrd[ng/2];
	int nChb=9;
	Cheby resOmchb(Jrgrd,resOmgrd,ng,nChb);
	Cheby hNchb(Jrgrd,hNgrd,ng,nChb);
	double resOm,resOmp; resOmchb.unfitderiv(Jrres,resOm,resOmp);
	if(resN[2]*resOmp>=0) return false;
	int rep=0,repmax=50;
	while(fabs(resOm)>1e-5*Omres[0] && Jrres>Jrgrd[0]){//Refine Jrres by Newton-Raphson
		double dJ=resOm/resOmp;
		while(dJ>.8*Jrres) dJ*=.6;
		Jrres-=dJ;
		if(Jrres<Jrgrd[0]) Jrres=Jrgrd[0];
		if(Jrres>Jrgrd[ng-1]) Jrres=Jrgrd[ng-1];
		resOmchb.unfitderiv(Jrres,resOm,resOmp);
		if(resOmp==0){
			printf("Derivative of resonant freq vanishes!\n");
			return false;
		}
		rep++; if(rep>repmax-20){
			printf("N-R: %g %g %g %g\n",resOm,resOmp,dJ,Jrres);
		}
		if(rep>repmax){
			printf("Newton-R failed\n");
			break;
		}
	}
	Jres=changeJr(Jres,Jrres-Jres[0]); resJp=prime_it(Jres,0);
	G=resOmp;
//value, gradient and 2nd deriv of hN
	hNchb.unfitderiv(Jrres,hres[0],hres[1],hres[2]);
	double Delta[3];
	getImin_max(); I=Imax; getDelta(t1c,Delta);
	printf("bottom: %g\n",resJp[0]+Delta[0]);
	printf("fix_derivs Om, G, hres: %d %f %f %f %f\n",resN[2],G,hres[0],hres[1],hres[2]);
//value of the non-resonant terms
	eT.AutoFit(Jres,Phi,bar,Omp,tolJ);
	return true;
}
int resTorus_L::setI(double Iin){
	/* This must be called after resTorus_L::resTorus_L() to specify
	 * amplitude of libration */
	circulate=true; setI(Iin,circulate);
	return 1;
}
int resTorus_L::setI(double Iin,int circulatein){
	/* This must be called after resTorus_L::resTorus_L() to specify
	 * amplitude of libration */
	I=Iin; circulate=circulatein; lJ=-1;
	if((G<0 && I>Imax) || (G>0 && I<Imin)){
		printf("I out of range: I, Imin, max %f %f %f\n",I,Imin,Imax);
		return 0;
	} else if((G<0 && I>Imin) || (G>0 && I<Imax)){
		circulate=false; set_Dtheta1();
		double Delta[3]; getDelta(t1c,Delta);
		off=MIN(0,resJp[0]+Delta[0]);
		printf("off %g\n",off);
	} else {
		Dt1=PI;//we are circulating
	}
	store_theta();
	printf("t1c, Dtheta circulate: %f %f %d\n",t1c,Dt1,circulate);
	return 1;
}
int resTorus_L::setID(double Iin,int circulatein){
	/* This must be called after resTorus_L::resTorus_L() to specify
	 * amplitude of libration */
	I=Iin; circulate=circulatein; lJ=-1;
	if((G<0 && I>Imax) || (G>0 && I<Imin)){
		printf("I out of range: I, Imin, max %f %f %f\n",I,Imin,Imax);
		return 0;
	} else if((G<0 && I>Imin) || (G>0 && I<Imax)){
		double Delta[3]; getDelta(t1c,Delta);
		if(Delta[0]<0 && resJp[0]<-Delta[0]){
			Dt1=PI; circulate=true;
		}else{
			set_Dtheta1(); circulate=false;
		}
	} else {
		Dt1=PI;//we are circulating
	}
	store_theta();
	printf("t1c, Dtheta circulate: %f %f %d\n",t1c,Dt1,circulate);
	return 1;
}
void resTorus_L::getImin_max(){
	double I0=-2*(hres[0]-hres[3])-pow(hres[1]-hres[4],2)/(.5*G-hres[2]+hres[5]);
	double I1= 2*(hres[0]+hres[3])-pow(hres[1]+hres[4],2)/(.5*G+hres[2]+hres[5]);
	Imax=MAX(I0,I1); Imin=MIN(I0,I1);
	if(Imax>0) Imax*=.999; else Imax*=1.0001;
	if(Imin<0) Imin*=.99999999; else Imin*=1.00000001;
	double D=-resJp[0],cos_tmin=-.5*G*D/(hres[2]*D+hres[1]);
	if(cos_tmin>0 && cos_tmin<=1)
		Imax=((.5*G+hres[2]*cos_tmin)*D+2*hres[1]*cos_tmin)*D+2*hres[0]*cos_tmin;
	Imax=2.1e-5;
	printf("resTorus_L: Imin Imax cos(tmin): %g %g %f\n",Imin,Imax,cos_tmin);
}
void resTorus_L::getDelta(double t1p,double *Delta){//returns J1p-resJp[0]
	double cost1p=cos(t1p),cos2t1p=cos(2*t1p);
	double a=.5*G+hres[2]*cost1p+hres[5]*cos2t1p;
	double b=2*(hres[1]*cost1p+hres[4]*cos2t1p);
	double c=2*(hres[0]*cost1p+hres[3]*cos2t1p)-I;
	if(quadratic(a,b,c,Delta)){
		printf("%g %g\n",t1p,hres[3]);
		Delta[0]=Delta[1]=Delta[2]=0;
	}
	if(cost1p<0) Delta[0]-=off*pow(cost1p,2);
}
double resTorus_L::dDelta(double t1p,double Delta){//returns dJ1p/dt1p
	double cost1p=cos(t1p),cos2t1p=cos(2*t1p);
	double sint1p=sin(t1p),sin2t1p=sin(2*t1p);
	double a=.5*G+hres[2]*cost1p+hres[5]*cos2t1p;
	double dadtp=-hres[2]*sint1p-2*hres[5]*sin2t1p;
	double b=2*(hres[1]*cost1p+hres[4]*cos2t1p);
	double dbdtp=-2*(hres[1]*sint1p+2*hres[4]*sin2t1p);
//	double c=2*(hres[0]*cost1p+hres[3]*cos2t1p)-I;
	double dcdtp=-2*(hres[0]*sint1p+2*hres[3]*sin2t1p);
	return -(Delta*(dadtp*Delta+dbdtp)+dcdtp)/(2*a*Delta+b);
}
#define TINY 1e-6
void resTorus_L::getJs(Actions &J1,Actions &J2,double f){//for f=1 returns extreme unperturbed actions
	double Delta[3]; getDelta(t1c,Delta);
	printf("getJs: %f\n",resJp[0]+Delta[0]);
	J1[0]=MAX(TINY,resJp[0]+f*Delta[0]); for(int i=1;i<3;i++) J1[i]=resJp[i];
	J1=prime_it(J1,1);//remove prime
	J2[0]=MAX(TINY,resJp[0]+f*Delta[1]); for(int i=1;i<3;i++) J2[i]=resJp[i];
	J2=prime_it(J2,1);
}
void resTorus_L::Dtheta1_fn(double t1p,double &f,double &fp) const{
	double cost1p=cos(t1p),cos2t1p=cos(2*t1p);
	double a=.5*G+hres[2]*cost1p+hres[5]*cos2t1p;
	double b=2*(hres[1]*cost1p+hres[4]*cos2t1p);
	double c=2*(hres[0]*cost1p+hres[3]*cos2t1p)-I;
	f=b*b-4*a*c;
	double sin2t1p=sin(t1p),sin4t1p=sin(2*t1p);
	double ap=-hres[2]*sin2t1p-hres[5]*2*sin4t1p;
	double bp=-2*(hres[1]*sin2t1p+hres[4]*2*sin4t1p);
	double cp=-2*(hres[0]*sin2t1p+hres[3]*2*sin4t1p);
	fp=2*b*bp-4*(ap*c+a*cp);
}
void resTorus_L::set_Dtheta1(){//find amplitude of excursions in theta1
	//first find limit to excursions
	double tol=1.e-12,f,f1,fp,t1;
	Dtheta1_fn(t1c,f,fp);
	if(f<=0){
		printf("resTorus_L: f<0: %f\n",f);
		Dt1=0; return;
	}
	if(t1c==0){
		t1=.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0){
			t1+=.05*(PI-t1); Dtheta1_fn(t1,f1,fp);
		}
	}else{
		t1=t1c-.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0){
			t1-=.1*t1; Dtheta1_fn(t1,f1,fp);
		}
	}			
	if(f1>=0){
		Dt1=PI; return;
	}
	if(fabs(f1)>tol) t1=rtsafe(this,&resTorus_L::Dtheta1_fn,t1c,t1,tol);
	if(t1c==0) Dt1=.999999*t1; else Dt1=.999999*fabs(PI-t1);
}
double resTorus_L::librationAction(){//returns libration J
	if(lJ!=-1) return lJ;
	int np=100,iD;
	double dt=Dt1/(float)(np);
	double Delta[3],t1p=t1c; getDelta(t1p,Delta);
	if(circulate) iD=(1+circulate)/2; else iD=2;
	double J1p=0,A=Delta[iD];
	for(int i=0;i<np;i++){
		t1p+=dt;
		getDelta(t1p,Delta);
		double B=Delta[iD];
		J1p+=.5*(A+B)*dt;
		A=B;
	}
	J1p/=PI;
	if(circulate) return resJp[0]+J1p;
	else return J1p;
}
double resTorus_L::theta_dot(double t1p,double D){//with D given
	return G*D+2*(hres[1]+hres[2]*D)*cos(t1p)+2*(hres[4]+hres[5]*D)*cos(2*t1p);
}
void resTorus_L::store_theta(){//Store what's needed by FullMap
	double D,Delta[3],t1p; 
	double s=0,x=0,dx=PI/(double)(nangle-1);
	if(!circulate){
		t1p=t1c; getDelta(t1p,Delta); D=Delta[1];
		double tdot=theta_dot(t1p,D);
		double A=cos(x)/tdot,B;
		t1grd[0]=t1p; trgrd[0]=s; Dgrd[0]=D;
		for(int i=1;i<nangle;i++){
			x+=dx;
			t1p=t1c+Dt1*sin(x); getDelta(t1p,Delta);
			if(i<nangle/2) D=Delta[1]; else D=Delta[0];
			tdot=theta_dot(t1p,D);
			if(tdot!=0){
				B=cos(x)/tdot;
				s+=.5*(A+B); A=B;
			}
			t1grd[i]=t1p; trgrd[i]=s; Dgrd[i]=D;
		}
	} else {
		int iD=(1+circulate)/2;
		t1p=0; getDelta(t1p,Delta); D=Delta[iD];
		double tdot=fabs(theta_dot(t1p,D));
		double A=1/tdot,B;
		t1grd[0]=t1p; trgrd[0]=s; Dgrd[0]=D;
		for(int i=1;i<nangle;i++){
			t1p+=dx; getDelta(t1p,Delta); D=Delta[iD];
			tdot=fabs(theta_dot(t1p,D));
			if(tdot!=0){
				B=1/tdot;
				s+=.5*(A+B); A=B;
			}
			t1grd[i]=t1p; trgrd[i]=s; Dgrd[i]=D;
		}
	}
	lOmega=PI/trgrd[nangle-1];
	for(int i=0;i<nangle;i++){ trgrd[i]*=lOmega;
	}
	if(circulate) lOmega=fabs(lOmega)/dx; else lOmega=fabs(lOmega)/(Dt1*dx);
}
void resTorus_L::check_angle(void){//printout to check libration_angle is working
	double D;
	printf("         res_theta            theta1p               Delta\n");
	for(int i=0;i<nangle;i+=5){
		double tr=from_librationAngle(trgrd[i],D);
		double lt=librationAngle(tr,D);
		printf("(%f %f) (%f %f) ",trgrd[i],lt,t1grd[i],tr);
		printf("(%f %f)\n",D,Dgrd[i]);
	}
	for(int i=0;i<nangle;i+=5){
		int ip=nangle-1-i;
		double tr=from_librationAngle(TPI-trgrd[ip],D);
		double lt=librationAngle(tr,D);
		printf("(%f %f) (%f %f) ",TPI-trgrd[ip],lt,TPI-t1grd[ip],tr);
		printf("(%f %f)\n",Dgrd[ip],D);
	}
}
double resTorus_L::librationAngle(double t1p,double Delta){//returns libration theta
	if(!circulate){
		while(t1p-t1c>PI) t1p-=TPI; while(t1p-t1c<-PI) t1p+=TPI;//now t1c-PI<t1p<t1c+PI
		double t1=t1c+fabs(t1c-t1p);//now t1c<t1<TPI
		double psi,td=theta_dot(t1,Delta);
		int nh=nangle/2;
		if((G>0 && td>0) || (G<0 && td<0))//heading out from t1c
			psi=lntp(t1grd,trgrd,nh,fabs(t1));//look in table up to PI/2
		else
			psi=lntp(t1grd+nh,trgrd+nh,nh,fabs(t1));//look in 2nd half of table
		if(t1p<t1c) psi=TPI-psi;
		return psi;
	} else {
		while(t1p>TPI) t1p-=TPI; while(t1p<0) t1p+=TPI;//now 0<tp1/TPI
		if(t1p<PI) return lntp(t1grd,trgrd,nangle,t1p);
		else return TPI-lntp(t1grd,trgrd,nangle,TPI-t1p);
	}
}
double resTorus_L::from_librationAngle(double tl,double &D){
	// Return t_p and action offset Jp-resJp[0] given the true libration angle tl
	while(tl>TPI) tl-=TPI; while(tl<0) tl+=TPI;// now 0<= tl <= 2*PI
	bool sn=false;
	if(tl>PI){
		tl=TPI-tl; sn=true;//now 0<= tl <= PI
	}
	int bot=0,top=nangle-1;
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((trgrd[top]-tl)*(tl-trgrd[n])>=0) bot=n;
		else top=n;
	}
	double f=(trgrd[top]-tl)/(trgrd[top]-trgrd[bot]);
	double t1p=(f*t1grd[bot]+(1-f)*t1grd[top]);
	if(circulate){
		if(sn) t1p=TPI-t1p;
	}else{
		if(sn) t1p=2*t1c-t1p;
	}
	D=f*Dgrd[bot]+(1-f)*Dgrd[top];
	return t1p;
}
double resTorus_L::from_librationAngle(double tl,double &dt1dla,double &D){
	// as above but with derivative of t1 wrt libration angle
	while(tl>TPI) tl-=TPI; while(tl<0) tl+=TPI;// now 0<= tl <= 2*PI
	bool sn=false;
	if(tl>PI){
		tl=TPI-tl; sn=true;//now 0<= tl <= PI
	}
	int bot=0,top=nangle-1;
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((trgrd[top]-tl)*(tl-trgrd[n])>=0) bot=n;
		else top=n;
	}
	double f=(trgrd[top]-tl)/(trgrd[top]-trgrd[bot]);
	double t1p=(f*t1grd[bot]+(1-f)*t1grd[top]);
	int r,s,t;
	if((f>.5 && bot>0) || (f<=.5 && top==nangle-1)){//t1 nearer bot
		r=bot-1; s=bot; t=top;
	}else if((f<=.5 && top<nangle-1) || (f>.5 && bot==0)){
		r=bot; s=top; t=top+1;
	}
	double a=(t1grd[t]-t1grd[r])/(trgrd[t]-trgrd[r]);
	double b=(t1grd[s]-t1grd[r]-a*(trgrd[s]-trgrd[r]))/
		  ((trgrd[s]-trgrd[r])*(trgrd[t]-trgrd[s]));
	dt1dla=a+b*(trgrd[t]+trgrd[r]-2*tl);
	if(circulate){
		if(sn) t1p=TPI-t1p;
	}else{
		if(sn) t1p=2*t1c-t1p;
	}
	D=f*Dgrd[bot]+(1-f)*Dgrd[top];
	return t1p;
}
Actions resTorus_L::fixJ(const Actions Jp,const Angles& A){
	//Adds non-resonant terms to primed actions with unprimed angles
	Frequencies Om=eT.omega();
	Actions Jfix=Jp;
	int10 iN=eT.i1(),jN=eT.i2(),kN=eT.i3();
	vec10 hN=eT.hn();
	for(int i=1;i<tmax;i++){
		if((iN[i]==resN[0] && kN[i]==resN[2]) || (iN[i]==-resN[0] && kN[i]==-resN[2])) continue;//exclude resonant line
		double cw=2*hN[i]*cos(iN[i]*A[0]+kN[i]*A[2])/Om(2); 
		Jfix[0]-=iN[i]*cw/(double)(kN[i]*resN[0]-iN[i]*resN[2]);
		Jfix[2]-=cw;
	}
	//if(Jfix[0]<0) printf("J[0]<0 in fixJ %f\n",Jfix[0]);
	return Jfix;
}
Actions resTorus_L::fixJ(const Actions Jp,const Angles& A,
			 double dJ1pdt1p,Matrix<double,3,3> &dJdt){
	//As above but computes derivatives of primed actions wrt
	//unprimed angles
	Frequencies Om=eT.omega();
	Actions Jfix=Jp;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++) dJdt[i][j]=0;
	}
	int10 iN=eT.i1(),jN=eT.i2(),kN=eT.i3();
	vec10 hN=eT.hn();
	for(int i=1;i<tmax;i++){
		if((iN[i]==resN[0] && kN[i]==resN[2]) || (iN[i]==-resN[0] && kN[i]==-resN[2])) continue;//exclude resonant line
		//if(jN(i)) continue;//exclude variation in Jz
		double th=iN[i]*A[0]+kN[i]*A[2];
		double cw=cos(th),sw=sin(th);
		double fac=2*hN(i)*cw/Om(2),fas=2*hN(i)*sw/Om(2);
		Jfix[0]-=iN[i]*fac/(double)(kN[i]*resN[0]-iN[i]*resN[2]);
		dJdt[0][0]+=iN[i]*fas/(double)(kN[i]*resN[0]-iN[i]*resN[2])*iN[0]/(double)resN[0];
		dJdt[0][2]+=iN[i]*fas/(double)(kN[i]*resN[0]-iN[i]*resN[2])*2;
		Jfix[2]-=fac;
		dJdt[2][0]+=fas*iN[0]/(double)resN[0];
		dJdt[2][2]+=fas*2;
	}
	//dJdt[0][0]+=resN[0]*dJ1pdt1p;
	return Jfix;
}
GCY resTorus_L::FullMap(const Angles &Ap){
	/* First angle is libration angle, then theta_z, theta_phi.
	 * The main job is to get theta_r from lib angle.  */
	Actions Jp=resJp; Angles A=Ap;
	double D,t1p=from_librationAngle(Ap[0],D);
	A[0]=(t1p-resN[2]*A[2])/(double)resN[0];//thetar
	Jp[0]+=D;
	Jp=fixJ(Jp,A);
	Torus T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	GCY gcy=T.FullMap(A);
	if(Jp[0]<1e-4) gcy[0]=-1;
	return gcy;
}
double resTorus_L::librationOmega(void){
	return lOmega;
}
void resTorus_L::check_derivs(Angles A){
	double tlib=A[0],dA=.02;
	double dtrdla,D,t1p=from_librationAngle(tlib,dtrdla,D);
	double dJ1pdt1p=dDelta(t1p,D);
	A[0]=(t1p-resN[2]*A[2])/(double)resN[0];//thetar
	Matrix<double,3,3> dQdt,dJdt;
	Actions Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
	Torus T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	dJdt=unprime(dJdt);
	T.get_derivs(A,dJdt,dQdt);
	dQdt=prime_deriv(dQdt);
	for(int i=0;i<3;i++) dQdt[i][0]*=dtrdla;
	A[1]+=dA;
	Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
	T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	GCY gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
	double dphi=gcy[2],dz=gcy[1];
	A[1]-=2*dA;
	Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
	T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
	dphi-=gcy[2]; dz-=gcy[1];
	//printf("(%f %f) (%f %f) ",.5*dz/dA,dQdt[1][1],.5*dphi/dA,dQdt[2][1]);
	A[1]+=dA; A[2]+=dA;
	A[0]=(t1p-resN[2]*A[2])/(double)resN[0];//thetar
	Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
	T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
	dphi=gcy[2]; dz=gcy[1];
	A[2]-=2*dA;
	A[0]=(t1p-resN[2]*A[2])/(double)resN[0];//thetar
	Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
	T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
	dphi-=gcy[2]; dz-=gcy[1];
	printf("(%f %f) (%f %f)\n",.5*dz/dA,dQdt[1][2],.5*dphi/dA,dQdt[2][2]);
}
void resTorus_L::SOS(ostream& sfile,int np){//work round with libration angle
	double tlib,z0=0,phi0=0,oldChi=100;
	Matrix<double,3,3> dQdt,dJdt;
	Angles A;  GCY gcy;
	for(int i=0;i<np;i++){
		if(t1c==0) tlib=i/(double)(np-1)*(PI-1.e-6);
		else tlib=i/(double)(np-1)*(TPI-1.e-6);
		if(oldChi>1e-4){ A[1]=0; A[2]=0;}
		double dtrdla,D,t1p=from_librationAngle(tlib,dtrdla,D);
		double dJ1pdt1p=dDelta(t1p,D);
		oldChi=100;
		bool eflag=false;
		for(int k=0;k<20;k++){
			A[0]=(t1p-resN[2]*A[2])/(double)resN[0];//thetar
			Actions Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,A,dJ1pdt1p,dJdt);
			if(Jp[0]<=1e-4) eflag=true;
			Torus T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
			gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
			double dphi=phi0-gcy[2],dz=z0-gcy[1];
			double newChi=dz*dz+dphi*dphi;
			if(newChi<1.e-8 || (k>15 && newChi>1.1*oldChi)) break;
			dJdt=unprime(dJdt);
			T.get_derivs(A,dJdt,dQdt);
			dQdt=prime_deriv(dQdt);
			double det=dQdt[1][1]*dQdt[2][2]-dQdt[1][2]*dQdt[2][1];
			double dA1=(dz*dQdt[2][2]-dphi*dQdt[1][2])/det;
			double dA2=(dphi*dQdt[1][1]-dz*dQdt[2][1])/det;
			while(MAX(fabs(dA1),fabs(dA2))>.05){dA1*=.5; dA2*=.5;}
			A[1]+=dA1; A[2]+=dA2; tidy_angles(A);
			oldChi=newChi;
		}
		if(eflag) gcy[0]=-1;
		sfile << std::setw(12) << gcy << "\n";
		if(gcy[5]<0) printf("%f %f\n",gcy[0],gcy[5]);
	}
}/*
void resTorus_L::SOS(ostream& sfile,int np){//work round with libration angle
	double Delta[3],tlib;
	for(int i=0;i<np;i++){
		if(t1c==0) tlib=i/(double)(np-1)*(PI-1.e-6);
		else tlib=i/(double)(np-1)*(TPI-1.e-6);
		Actions Jp=resJp;
		double D,t1p=from_librationAngle(tlib,D);
		Angles A; A[1]=0; A[2]=0;
		A[0]=(t1p-resN[2]*A[2])/(double)resN[0];
		Jp[0]+=D;
		Actions Jtry=fixJ(Jp,A);
		Torus T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jtry);
		GCY gcy=T.FullMap(A);
		while(gcy[1]<-0.01){
			A[1]+=.001; gcy=T.FullMap(A);
		}
		while(gcy[1]>0.01){
			A[1]-=.001; gcy=T.FullMap(A);
		}
		if(gcy[2]>PI) gcy[2]-=TPI;
		double Ttop,Tbot,Atop,Abot,dA=.1,eps=.0001;
		if(fabs(gcy[2])<eps){
			sfile << std::setw(12) << gcy <<'\n'; continue;
		}
		if(gcy[2]<0){
			Abot=A[2]; Tbot=gcy[2];
			while(gcy[2]<0){
				A[0]-=dA*resN[2]; A[2]+=dA*resN[0];
				Jtry=fixJ(Jp,A);
				T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jtry);
				gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
				if(gcy[2]<0){
					Abot=A[2]; Tbot=gcy[2];
				}
			}
			Atop=A[2]; Ttop=gcy[2];
		}else{
			Atop=A[2]; Ttop=gcy[2];
			while(gcy[2]>0){
				A[0]+=dA*resN[2]; A[2]-=dA*resN[0];
				Jtry=fixJ(Jp,A);
				T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jtry);
				gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
				if(gcy[2]>0){
					Atop=A[2]; Ttop=gcy[2];
				}
			}
			Abot=A[2]; Tbot=gcy[2];
		}
		while(1){
			double f=Tbot/(Tbot-Ttop);
			dA=(f*Atop+(1-f)*Abot)-A[2];
			A[0]-=(dA*resN[2])/(double)resN[0]; A[2]+=dA;
			Jtry=fixJ(Jp,A);
			T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jtry);
			gcy=T.FullMap(A); if(gcy[2]>PI) gcy[2]-=TPI;
			if(fabs(gcy[2])<eps) break;
			if(gcy[2]<0){
				Abot=A[2]; Tbot=gcy[2];
			}else{
				Atop=A[2]; Ttop=gcy[2];
			}
		}
		while(gcy[1]<-0.01){
			A[1]+=.001; gcy=T.FullMap(A);
		}
		while(gcy[1]>0.01){
			A[1]-=.001; gcy=T.FullMap(A);
		}
		sfile << std::setw(12) << gcy << "\n";
	}
}*/
bool resTorus_L::oldV(Velocity *Vout,Angles *Aout,int nv,const Velocity &V,
		      const Angles &A,double &stat){
	if(nv==0) return false;
	bool res=false; stat=100;
	double tolerance=.05,sq=0; for(int j=0;j<3;j++) sq+=pow(V[j],2);
	for(int i=0;i<nv;i++){
		double dot=0,sqo=0,diff=0;
		for(int j=0;j<3;j++){
			dot+=Vout[i][j]*V[j]; sqo+=pow(Vout[i][j],2);
			double dA=Aout[i][j]-A[j];
			if(dA>Pi) dA-=TPi; if(dA<-Pi) dA+=TPi;
			diff+=fabs(dA);
		}
		double st=sqrt(MAX(0,1-pow(dot,2)/(sq*sqo)))+diff; stat=MIN(stat,st);
		if(st<tolerance) res=true;
	}
	return res;
}
void resTorus_L::LevMstart(Angles theta,const Position Q,PSPT &QP3,Vector<double,3> &B,
			   Matrix<double,3,3> &A,Matrix<double,3,3> &dQdt,double &chisqo){
	//input angles true, with libration angle at position 0
	Matrix<double,3,3> dJdt;
	double dtrdla,D,t1p=from_librationAngle(theta[0],dtrdla,D);
	theta[0]=(t1p-resN[2]*theta[2])/(double)resN[0];//thetar
	double dJ1pdt1p=dDelta(t1p,D);
	Actions Jp=resJp; Jp[0]+=D; Jp=fixJ(Jp,theta,dJ1pdt1p,dJdt);
	dJdt=unprime(dJdt);
	Torus T=InterpTorus_n2(Tgrid,nr,Jbar_grid,dJ_grid,Jp);
	T.LevCof3DTrue(theta,dJdt,Q,QP3,chisqo,B,dQdt);
	B=prime_deriv(B); dQdt=prime_deriv(dQdt);
	B[0]*=dtrdla;//convert derivs to libration angle
	for(int i=0;i<3;i++) dQdt[i][0]*=dtrdla;
	A[0][0] = pow(dQdt[0][0],2) + pow(dQdt[1][0],2)+ pow(Q(0)*dQdt[2][0],2);
	A[0][1] = dQdt[0][0]*dQdt[0][1] + dQdt[1][0]*dQdt[1][1] + pow(Q(0),2)*dQdt[2][0]*dQdt[2][1];
	A[1][0] = A[0][1];
	A[0][2] = dQdt[0][0]*dQdt[0][2] + dQdt[1][0]*dQdt[1][2] + pow(Q(0),2)*dQdt[2][0]*dQdt[2][2];
	A[2][0] = A[0][2];
	A[1][1] = pow(dQdt[0][1],2) + pow(dQdt[1][1],2) + pow(Q(0)*dQdt[2][1],2);
	A[1][2] = dQdt[0][1]*dQdt[0][2] + dQdt[1][1]*dQdt[1][2] + pow(Q(0),2)*dQdt[2][1]*dQdt[2][2];
	A[2][1] = A[1][2];
	A[2][2] = pow(dQdt[0][2],2) + pow(dQdt[1][2],2) + pow(Q(0)*dQdt[2][2],2);
}
void resTorus_L::LevMstep(Angles &theta,const Position Q,PSPT &QP3,Vector<double,3> &B,Matrix<double,3,3> &A,
			  Matrix<double,3,3> &dQdt,double &chisqo,double &lam){
	double chisq,BB[3],**AA; AA=dmatrix(3,3);
	Vector<double,3> Btry;
	Matrix<double,3,3> Atry,dQdtry;
	for(int i=0;i<3;i++){
		BB[i]=B[i];
		for(int j=0;j<3;j++){
			AA[i][j]=A[i][j];
			if(j==i) AA[i][j]*=1+lam;
		}
	}
	GaussJordan(AA,3,BB);
	Angles t_try;
	for(int i=0;i<3;i++){
		t_try[i]=theta[i]+BB[i];
	}
	LevMstart(t_try,Q,QP3,Btry,Atry,dQdtry,chisq);
	if(chisq<chisqo  && !std::isnan(Btry(0))) {
		lam *= 0.125;
		chisqo = chisq; theta = t_try; tidy_angles(theta);
		A = Atry; B = Btry; dQdt = dQdtry;
	} else {
		lam *= 8.;
	}
	delmatrix(AA,3);
}
int resTorus_L::containsPoint(const Position &Q,Velocity *Vout,Angles *Aout,double *densout){
	const    int    maxit1=50, maxtry=100;
	const    double tiny=1.e-7;

	double rq = Q[0]*Q[0]+Q[1]*Q[1],rtin = sqrt(rq)*tiny; // tolerance in position 
	double chisqo,stat;
	Matrix<double,3,3> A,dQdt,dJtry,Atry,dQdtry;;
	Vector<double,3> B,Btry;
	Angles theta;
	PSPT QP3;
	int nfail=0,nv=0,ntry=0;
	Random3 ran(7986535);
	while(nfail<maxtry && ntry<maxtry){
		theta[0] = TPi*(ran.RandomDouble());
		theta[1] = TPi*(ran.RandomDouble());
		theta[2] = TPi*(ran.RandomDouble());
		LevMstart(theta,Q,QP3,B,A,dQdt,chisqo);
		int it=0;
		double lam=0.5;
		while(chisqo>rtin && maxit1>it++ && lam < 1.e10 ) {  // Lev Mar iteration
			LevMstep(theta,Q,QP3,B,A,dQdt,chisqo,lam);
		}
		if(lam<1e-4 && chisqo > 1e-5){
			printf("iTorus::containsPoint no can do: %d %g\n",it,lam);
			return 0;
		}
		if(chisqo<rtin){
			for(int i=0;i<3;i++) Vout[nv][i]=QP3(i+3);
			Aout[nv]=theta;
			if(!oldV(Vout,Aout,nv,Vout[nv],Aout[nv],stat)){
				densout[nv]=fabs(LUDet3(dQdt));
				//densout[nv]=fabs(dQdt[0][0]*dQdt[1][1]-dQdt[0][1]*dQdt[1][0])*Q[0];
				nv++;
			} else nfail++;
		}
		ntry++; nfail++;
	}
	return nv;
}
