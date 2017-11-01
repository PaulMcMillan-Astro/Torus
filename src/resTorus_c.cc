/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
// Functionality for 0,0,Nphi resonantly trapped tori per Binney (2017)
 
#include "resTorus_c.h"

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SGN(A)  ((A)>0? 1:-1)
#define SMALL 1e-7

static double PI=acos(-1),TPI=2*PI,PIH=.5*PI;

resTorus_c::resTorus_c(Torus **Tgridin,int npin,Actions Jresin,Potential *Phi,bar_pot *bar,double Omp,double tolJ){
	Tgrid=Tgridin; np=npin; resN=0; resN[2]=2;
	findJp(Jresin,Phi,bar,Omp,tolJ); resJp=Jresin;
	int ok=eT.AutoFit(resJp,Phi,bar,Omp,tolJ);
	nangle=200;//nangle must be even
	t1grd=new double[nangle]; tpgrd=new double[nangle]; Dgrd=new double[nangle];
	hres[0]=get_driver(eT);	for(int i=1;i<6;i++) hres[i]=0;
	if(hres[0]>0) t1c=0;
	else t1c=PIH;
	double m2=get_OLR(eT);
	Actions J=resJp;
	double dJ=.05*J[2];
	J[2]+=dJ;
	Torus T1; T1.AutoFit(J,Phi,tolJ);
	Frequencies Om=T1.omega(); Om[2]-=Omp;
	double resOm=resCond(Om);
	dJ*=2;
	J[2]-=dJ;
	Torus T0; T0.AutoFit(J,Phi,tolJ);
	Om=T0.omega(); Om[2]-=Omp; resOm-=resCond(Om);
	m2*=2/Om(0);//only a rough estimate of non-res change to Jphi
	G=resOm/dJ;
	getImin_max();	setI(.99*Imin);
	Actions J0,J1; getJs(J0,J1,1);
	fix_derivs(resJp,J0,J1,Phi,bar,Omp,tolJ);
	//printf("Non res dJp, G, hres: %f %f %g %f %f\n",m2,G,hres[0],hres[1],hres[2]);
	getImin_max();
	double safety=10;
	getJs(J0,J1,1);
	J0[0]=MAX(J0[0],.001);
	J0[2]-=safety*fabs(m2); J1[2]+=safety*fabs(m2);
	dJ=.35*resJp[0];
	dJ_grid[2]=(J1[2]-J0[2])/(double)(np-1);
	dJ_grid[1]=0; dJ_grid[0]=2*dJ;
	Jbar_grid[2]=J0[2]+.5*dJ_grid[2];
	Jbar_grid[1]=J0[1]; Jbar_grid[0]=J0[0];
	printf("Grid Js: (%f %f) (%f %f)\n",J0(0),J0(2),J1(0),J1(2));
	vec4 tm;
	for(int k=0;k<2;k++){
		if(k==0){
			J0[0]-=dJ; J1[0]-=dJ;
		}else{
			J0[0]+=2*dJ; J1[0]+=2*dJ;
		}
		for(int i=0;i<np;i++){
			double fac=i/(double)(np-1);
			Actions Jn=(1-fac)*J0+fac*J1;
			if(i>100){//to implement make if(i>0)
				Tgrid[k][i].FitWithFixToyPot(Jn,tm,Phi,tolJ);
			}else{
				Tgrid[k][i].AutoFit(Jn,Phi,tolJ);
				tm=Tgrid[k][i].TP();
			}
		}
	}
}
void resTorus_c::findJp(Actions &J,Potential *Phi,bar_pot *bar,double Omp,double tolJ){
	Torus T; T.AutoFit(J,Phi,tolJ);
	Frequencies Om=T.omega(); Om[2]-=Omp;
	double dJ=.05*J[2], Jmax=0, Jmin=0, resOm=resCond(Om);
	while(Jmax*Jmin==0){
		if(resOm>0){
			Jmin=J[2]; J[2]+=dJ;
		}else{
			Jmax=J[2]; J[2]-=dJ;
		}
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); Om[2]-=Omp; resOm=resCond(Om);
	}
	while(fabs(resOm)>10*SMALL && Jmax-Jmin>10*SMALL){
		J[2]=.5*(Jmax+Jmin);
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); Om[2]-=Omp; resOm=resCond(Om);
		//printf("2: %g %f %f\n",resOm,Jmin,Jmax);
		if(resOm>0) Jmin=J[2];
		else Jmax=J[2];
	}
}
double resTorus_c::get_driver(eTorus &eT){
	for(int i=0;i<10;i++)
		if(eT.i1()(i)==resN(0) && eT.i2()(i)==resN(1) && eT.i3()(i)==resN(2)) return eT.hn()(i);
}
double resTorus_c::get_H(eTorus &eT){
	for(int i=0;i<10;i++)
		if(eT.i1()(i)==0 && eT.i2()(i)==0 && eT.i3()(i)==0) return eT.hn()(i);
}
double resTorus_c::get_OLR(eTorus &eT){
	for(int i=0;i<10;i++)
		if(eT.i1()(i)==1 && eT.i2()(i)==0 && eT.i3()(i)==2) return eT.hn()(i);
}
bool resTorus_c::fix_derivs(Actions &Jres,Actions &J0,Actions &J1,
			    Potential *Phi,bar_pot *bar,double Omp,double tolJ){
	int ng=6;
/*	Torus T,Tgrd[2];
	double Jpgrid[2]={J0[2],J1[2]};
	Tgrd[0].AutoFit(J0,Phi,tolJ); Tgrd[1].AutoFit(J1,Phi,tolJ);
	Frequencies Om=Tgrd[0].omega(); Om[2]-=Omp;
	Frequencies Omres=.5*(Tgrd[0].omega()+Tgrd[1].omega());
	Omres[2]-=Omp;*/
	//printf("Omres: %g\n",Omres(2));
	double Jpgrd[ng],Omgrd[ng],resOmgrd[ng],Hgrd[ng],hNgrd[ng],y2[ng];
	for(int i=0;i<ng;i++){
		double f=i/(float)(ng-1);
		Actions J=(1-f)*J0+f*J1;
		Jpgrd[i]=J[2];
		eTorus eTT(J,Phi,bar,Omp,tolJ);
		Frequencies Om=eTT.omega();
		Omgrd[i]=Om(0); resOmgrd[i]=resCond(Om);
		Hgrd[i]=get_H(eTT); hNgrd[i]=get_driver(eTT);
	}
	double Jpres=Jres[2];
	Frequencies Omres=Omgrd[ng/2];
	int nChb=9;
	Cheby resOmchb(Jpgrd,resOmgrd,ng,nChb);
	Cheby hNchb(Jpgrd,hNgrd,ng,nChb);
	double resOm,resOmp; resOmchb.unfitderiv(Jpres,resOm,resOmp);
	if(resOmp>=0) return false;
	int rep=0,repmax=50;
	//printf("refining...");
	while(fabs(resOm)>1e-5*Omres[0] && Jpres>Jpgrd[0]){//Refine Jpres by Newton-Raphson
		double dJ=resOm/resOmp;
		while(dJ>.8*Jpres) dJ*=.6;
		Jpres-=dJ;
		if(Jpres<Jpgrd[0]) Jpres=Jpgrd[0]; else if(Jpres>Jpgrd[ng-1]) Jpres=Jpgrd[ng-1];
		resOmchb.unfitderiv(Jpres,resOm,resOmp);
		if(resOmp==0){
			printf("Derivative of resonant freq vanishes!\n");
			return false;
		}
		rep++; if(rep>repmax-20){
			printf("N-R: %g %g %g %g\n",resOm,resOmp,dJ,Jpres);
		}
		if(rep>repmax){
			printf("Newton-R failed\n");
			break;
		}
	}
	Jres[2]=Jpres;
	G=resOmp;
//value, gradient and 2nd deriv of hN
	hNchb.unfitderiv(Jpres,hres[0],hres[1],hres[2]);
//value of the non-resonant terms
	eT.AutoFit(Jres,Phi,bar,Omp,tolJ);
	return true;
}
int resTorus_c::setI(double Iin){
	/* This must be called after resTorus_c::resTorus_c() to specify
	 * amplitude of libration */
	inout=1; setI(Iin,inout);
	return 1;
}
int resTorus_c::setI(double Iin,int inoutin){
	/* This must be called after resTorus_c::resTorus_c() to specify
	 * amplitude of libration */
	I=Iin; inout=inoutin; lJ=-1; 
	if(I>Imax){
		printf("I out of range: I, Imin, max %f %f %f\n",I,Imin,Imax);
		return 0;
	} else if(I>Imin){
		//printf("in setI %f %f\n",I,t1c);
		inout=0; set_Dtheta1();
	} else {
		Dt1=PIH;
		//printf("circulating\n");
	}
	//printf("t1c, Dtheta Jl: %f %f %f\n",t1c,Dt1,librationAction());
	store_theta(); return 1;
}
void resTorus_c::getImin_max(){
	double I0=-2*(hres[0]-hres[3])-pow(hres[1]-hres[4],2)/(.5*G-hres[2]+hres[5]);
	double I1= 2*(hres[0]+hres[3])-pow(hres[1]+hres[4],2)/(.5*G+hres[2]+hres[5]);
	Imax=.99999999*MAX(I0,I1); Imin=MIN(I0,I1);
	if(Imin<0) Imin*=.99999999; else Imin*=1.00000001;
//	printf("resTorus_c: Imin Imax: %g %g\n",Imin,Imax);
}
void resTorus_c::getDelta(double t1p,double *Delta){//returns J1p-resJp[0]
	double cos2t1p=cos(2*t1p),cos4t1p=cos(4*t1p);
	double a=.5*G+hres[2]*cos2t1p+hres[5]*cos4t1p;
	double b=2*(hres[1]*cos2t1p+hres[4]*cos4t1p);
	double c=2*(hres[0]*cos2t1p+hres[3]*cos4t1p)-I;
	quadratic(a,b,c,Delta);
}
double resTorus_c::dDelta(double t1p,double Delta){//returns dJ1p/dtp
	double cos2t1p=cos(2*t1p),cos4t1p=cos(4*t1p);
	double sin2t1p=sin(2*t1p),sin4t1p=sin(4*t1p);
	double a=.5*G+hres[2]*cos2t1p+hres[5]*cos4t1p;
	double dadtp=-2*hres[2]*sin2t1p-4*hres[5]*sin4t1p;
	double b=2*(hres[1]*cos2t1p+hres[4]*cos4t1p);
	double dbdtp=-2*(2*hres[1]*sin2t1p+4*hres[4]*sin4t1p);
//	double c=2*(hres[0]*cos2t1p+hres[3]*cos4t1p)-I;
	double dcdtp=-2*(2*hres[0]*sin2t1p+4*hres[3]*sin4t1p);
	return -(Delta*(dadtp*Delta+dbdtp)+dcdtp)/(2*a*Delta+b);
}
#define TINY 1e-6
void resTorus_c::getJs(Actions &J1,Actions &J2,double f){//return f times  extreme unperturbed actions for libration
	double Delta[3]; getDelta(t1c,Delta);
	for(int i=0;i<2;i++) J1[i]=resJp[i]; J1[2]=MAX(TINY,resJp[2]+f*Delta[0]);
	for(int i=0;i<2;i++) J2[i]=resJp[i]; J2[2]=MAX(TINY,resJp[2]+f*Delta[1]);
}
void resTorus_c::Dtheta1_fn(double t1p,double &f,double &fp) const{
	double cos2t1p=cos(2*t1p),cos4t1p=cos(4*t1p);
	double a=.5*G+hres[2]*cos2t1p+hres[5]*cos4t1p;
	double b=2*(hres[1]*cos2t1p+hres[4]*cos4t1p);
	double c=2*(hres[0]*cos2t1p+hres[3]*cos4t1p)-I;
	f=b*b-4*a*c;
	double sin2t1p=sin(2*t1p),sin4t1p=sin(4*t1p);
	double ap=-hres[2]*2*sin2t1p-hres[5]*4*sin4t1p;
	double bp=-2*(hres[1]*2*sin2t1p+hres[4]*4*sin4t1p);
	double cp=-2*(hres[0]*2*sin2t1p+hres[3]*4*sin4t1p);
	fp=2*b*bp-4*(ap*c+a*cp);
}
bool resTorus_c::set_Dtheta1(){//find amplitude of excursions in theta1
	//Max range is +-PIH because n=2. first find limit to excursions
	double tol=1.e-12,f,f1,fp,t1;
	Dtheta1_fn(t1c,f,fp);
	int nc=0;
	if(f<=0){
		printf("resTorus_c: f<0: %f\n",f);
		Dt1=0; return 0;
	}
	if(t1c==0){
		t1=.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0){
			t1+=.05*(PIH-t1); Dtheta1_fn(t1,f1,fp);
		}
	}else{
		t1=t1c-.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0 && nc<100){
			t1-=.1*t1; Dtheta1_fn(t1,f1,fp);
			//if(nc>50) printf("%f %g\n",t1,f1); nc++;
		}
	}
	if(f1>=0){
		Dt1=PIH; return 1;
	}
	if(fabs(f1)>tol) t1=rtsafe(this,&resTorus_c::Dtheta1_fn,t1c,t1,tol);
	if(t1c==0) Dt1=.999999*t1; else Dt1=.999999*fabs(PIH-t1);
	return 1;
}
double resTorus_c::librationAction(){//returns libration J
	if(lJ!=-1) return lJ;
	int npt=100,iD;
	double dt=Dt1/(float)(npt);
	double Delta[3],t1p=t1c; getDelta(t1p,Delta);
	if(inout) iD=(1+inout)/2; else iD=2;
	double J1p=0,A=Delta[iD];
	for(int i=0;i<npt;i++){
		t1p+=dt;
		getDelta(t1p,Delta);
		double B=Delta[iD];
		J1p+=.5*(A+B)*dt;
		A=B;
	}
	J1p/=PI;
	if(inout) return resJp[2]+2*J1p;
	else return J1p;
}
double resTorus_c::theta_dot(double t1p,double D){//with D given
	return G*D+2*(hres[1]+hres[2]*D)*cos(2*t1p)+2*(hres[4]+hres[5]*D)*cos(4*t1p);
}
void resTorus_c::store_theta(){//Store what's needed by FullMap
	double D,Delta[3],t1p; 
	double s=0,x=0,dx=PI/(double)(nangle-1);
	if(!inout){
		t1p=t1c;
		getDelta(t1p,Delta); D=Delta[1];
		double tdot=theta_dot(t1p,D);
		double A=cos(x)/tdot,B;
		t1grd[0]=t1p; tpgrd[0]=s; Dgrd[0]=D;
		for(int i=1;i<nangle;i++){
			x+=dx;
			t1p=t1c+Dt1*sin(x); getDelta(t1p,Delta);
			if(i<nangle/2) D=Delta[1]; else D=Delta[0];
			tdot=theta_dot(t1p,D);
			if(tdot!=0){
				B=cos(x)/tdot;
				s+=.5*(A+B); A=B;
			}
			t1grd[i]=t1p; tpgrd[i]=s; Dgrd[i]=D;
		}
	} else {
		int iD=(1+inout)/2;
		t1p=0;
		getDelta(t1p,Delta); D=Delta[iD];
		double tdot=theta_dot(t1p,D);
		double A=1/tdot,B;
		t1grd[0]=t1p; tpgrd[0]=s; Dgrd[0]=D;
		for(int i=1;i<nangle;i++){
			t1p+=dx; getDelta(t1p,Delta); D=Delta[iD];
			tdot=theta_dot(t1p,D);
			if(tdot!=0){
				B=1/tdot;
				s+=.5*(A+B); A=B;
			}
			t1grd[i]=t1p; tpgrd[i]=s; Dgrd[i]=D;
		}
	}
	lOmega=PI/tpgrd[nangle-1];
	for(int i=0;i<nangle;i++) tpgrd[i]*=lOmega;
	if(inout) lOmega/=dx; else lOmega/=Dt1*dx;
}
void resTorus_c::check_angle(void){//printout to check libration_angle is working
	double D;
	printf("   Libration_theta         res_theta            Delta\n");
	for(int i=0;i<nangle;i+=5){
		double tr=from_librationAngle(tpgrd[i],D);
		double lt=librationAngle(tr,D);
		printf("(%f %f) (%f %f) ",tpgrd[i],lt,t1grd[i],tr);
		printf("(%f %f)\n",D,Dgrd[i]);
	}
	for(int i=0;i<nangle;i+=5){
		int ip=nangle-1-i;
		double tr=from_librationAngle(TPI-tpgrd[ip],D);
		double lt=librationAngle(tr,D);
		if(inout) printf("(%f %f) (%f %f) ",TPI-tpgrd[ip],lt,TPI-t1grd[ip],tr);
		else printf("(%f %f) (%f %f) ",TPI-tpgrd[ip],lt,PI-t1grd[i],tr);
		printf("(%f %f)\n",Dgrd[ip],D);
	}
}
double resTorus_c::librationAngle(double t1p,double Delta){//returns libration theta
	if(!inout){
		while(t1p-t1c>PI) t1p-=TPI; while(t1p-t1c<-PI) t1p+=TPI;//now t1c-PI<t1p<t1c+PI
		double t1=t1c+fabs(t1c-t1p);//now t1c<t1<PI
		double psi,td=theta_dot(t1,Delta);
		int nh=nangle/2;
		if(td<0)//heading out from t1c
			psi=lntp(t1grd,tpgrd,nh,fabs(t1));//look in table up to PI/2
		else
			psi=lntp(t1grd+nh,tpgrd+nh,nh,fabs(t1));//look in 2nd half of table
		if(t1p<t1c) psi=2*PI-psi;
		return psi;
	} else {
		while(t1p>TPI) t1p-=TPI; while(t1p<0) t1p+=TPI;
		if(t1p<PI) return lntp(t1grd,tpgrd,nangle,t1p);
		else return TPI-lntp(t1grd,tpgrd,nangle,TPI-t1p);
	}
}
double resTorus_c::from_librationAngle(double tl,double &D){
	// Return t_p and action offset Jp-resJp[2] given the true libration angle tl
	while(tl>TPI) tl-=TPI; while(tl<0) tl+=TPI;// now 0<= tl <= 2*PI
	bool sn=false;
	if(tl>PI){
		tl=TPI-tl; sn=true;//now 0<= tl <= PI
	}
	int bot=0,top=nangle-1;
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((tpgrd[top]-tl)*(tl-tpgrd[n])>=0) bot=n;
		else top=n;
	}
	double f=(tpgrd[top]-tl)/(tpgrd[top]-tpgrd[bot]);
	double t1p=(f*t1grd[bot]+(1-f)*t1grd[top]);
	if(inout){
		if(sn) t1p=TPI-t1p;
	}else{
		if(sn) t1p=2*t1c-t1p;
	}
	D=f*Dgrd[bot]+(1-f)*Dgrd[top];
	return t1p;
}
double resTorus_c::from_librationAngle(double tl,double &dt2dla,double &D){
	// as above but with derivative of t2 wrt libration angle
	while(tl>TPI) tl-=TPI; while(tl<0) tl+=TPI;// now 0<= tl <= 2*PI
	bool sn=false;
	if(tl>PI){
		tl=TPI-tl; sn=true;//now 0<= tl <= PI
	}
	int bot=0,top=nangle-1;
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((tpgrd[top]-tl)*(tl-tpgrd[n])>=0) bot=n;
		else top=n;
	}
	double f=(tpgrd[top]-tl)/(tpgrd[top]-tpgrd[bot]);
	double t1p=(f*t1grd[bot]+(1-f)*t1grd[top]);
	if(f>.5 && bot>0){//tl nearer bot
		int m=bot-1;
		double a=(t1grd[top]-t1grd[m])/(tpgrd[top]-tpgrd[m]);
		double b=(t1grd[bot]-t1grd[m]-a*(tpgrd[bot]-tpgrd[m]))/
			 ((tpgrd[top]-tpgrd[bot])*(tpgrd[bot]-tpgrd[m]));
		dt2dla=a+b*(tpgrd[top]+tpgrd[m]-2*tl);
	}else if(f<.5 && top<nangle-1){//nearer top
		int m=top+1;
		double a=(t1grd[m]-t1grd[bot])/(tpgrd[m]-tpgrd[bot]);
		double b=(t1grd[top]-t1grd[bot]-a*(tpgrd[top]-tpgrd[bot]))/
			 ((tpgrd[m]-tpgrd[top])*(tpgrd[top]-tpgrd[bot]));
		dt2dla=a+b*(tpgrd[m]+tpgrd[bot]-2*tl);
	} else {
		dt2dla=(t1grd[top]-t1grd[bot])/(tpgrd[top]-tpgrd[bot]);
	}
	if(inout){
		if(sn) t1p=TPI-t1p;
	}else{
		if(sn) t1p=2*t1c-t1p;
	}
	D=f*Dgrd[bot]+(1-f)*Dgrd[top];
	return t1p;
}
Actions resTorus_c::fixJ(const Actions Jp,const Angles& A,double tp){
	//Adds non-resonant terms
	Frequencies Om=eT.omega();
	Actions Jfix=Jp;
	int10 iN=eT.i1(),jN=eT.i2(),kN=eT.i3();
	vec10 hN=eT.hn();
	for(int i=1;i<10;i++){
		if(iN[i]==resN[0] && kN[i]==resN[2]) continue;//exclude resonant line
		double cw=2*hN[i]*cos(iN[i]*A[0]+kN[i]*tp)/Om(0); 
		Jfix[0]-=cw;
		Jfix[2]-=kN[i]*cw/(double)iN[i];
	}
	return Jfix;
}
Actions resTorus_c::fixJ(const Actions Jp,const Angles& A,double tp,
			 double dJ2dtp,Matrix<double,3,3> &dJdt){
	//As above but computes derivatives
	Frequencies Om=eT.omega();
	Actions Jfix=Jp;
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) dJdt[i][j]=0;
	dJdt[2][2]=dJ2dtp;
	int10 iN=eT.i1(),jN=eT.i2(),kN=eT.i3();
	vec10 hN=eT.hn();
	for(int i=1;i<10;i++){
		if(iN[i]==resN[0] && kN[i]==resN[2]) continue;//exclude resonant line
		if(jN(i)) continue;//exclude variation in Jz
		double th=iN(i)*A(0)+kN(i)*tp;
		double cw=cos(th),sw=sin(th);
		double fac=2*hN(i)*cw/Om(0),fas=2*hN(i)*sw/Om(0);
		if(fabs(fac)>0.4) printf("iTorus::fixJ %d %d %f\n",iN(i),kN(i),Om(0));
		Jfix[0]-=fac; Jfix[2]-=kN(i)*fac/(double)iN(i);
		dJdt[0][0]+=iN(i)*fas; dJdt[0][2]+=kN(i)*fas/(double)iN(i);
		dJdt[2][0]+=kN(i)*fas; dJdt[2][2]+=kN(i)*kN(i)*fas/(double)iN(i);
	}
	return Jfix;
}
GCY resTorus_c::FullMap(const Angles& A){
	/* Last angle is libration angle, second theta_z, first theta_r.
	 * The main job is to get theta_p from lib angle.  */
	Actions Jp=resJp;
	double D,t1p=from_librationAngle(A[2],D); Jp[2]+=D;
	Jp=fixJ(Jp,A,t1p);
	Angles theta=A; theta[2]=t1p;
	Torus T=InterpTorus_2n(Tgrid,np,Jbar_grid,dJ_grid,Jp);
	GCY gcy=T.FullMap(theta);
	return gcy;
}
double resTorus_c::librationOmega(void){
	return lOmega;
}
void resTorus_c::SOS(ostream& sfile,int npt){//work round with true angle
	npt/=2;
	double eps=.0001,A_old;
	Angles A;
	for(int krep=-1; krep<2; krep+=2){
		if(inout && krep==-1){npt*=2; continue;}
		A_old=krep*PIH;
		for(int i=0;i<npt;i++){
			//printf("(%d %d)",krep,i);
			A[0]=PIH+i/(double)(npt-1)*(TPI-1.e-6);
			A[1]=0;
			A[2]=A_old;
			GCY gcy=FullMap(A);
			while(gcy[1]<-0.01){
				A[1]+=.001; gcy=FullMap(A);
			}
			while(gcy[1]>0.01){
				A[1]-=.001; gcy=FullMap(A);
			}
			double diff=gcy[2]-PIH,diff_old=diff,dA=.1;
			int n=0;
			while(diff*diff_old>eps*eps){
				if(n>40)printf("(%g %g %g)\n",diff,diff_old,dA);
				if(n>40) exit(0); n++;
				if(n>1 && fabs(diff)>fabs(diff_old)){//we're going the wrong way
					dA=-dA; A[2]+=dA; diff=diff_old;//reverse,undo
				} else {
					diff_old=diff; A_old=A[2];
				}
				A[2]+=dA;
				gcy=FullMap(A);
				diff=gcy[2]-PIH;
			}//we straddle root
			while(fabs(diff)>eps && n<20){
				double A_new=A[2];
				A[2]=.5*(A_new+A_old);
				gcy=FullMap(A);
				double test=gcy[2]-PIH;
				if(diff*test>0){//replace last value A_new
					A_new=A[2]; diff=test;
				}else{//replace A_old
					A_old=A_new; diff_old=diff;
					A_new=A[2]; diff=test;
				}
			}
			sfile << std::setw(12) << gcy <<'\n';
		}
	}
}
bool resTorus_c::oldV(Velocity *Vout,Angles *Aout,int nv,const Velocity &V,
		      const Angles &A,double &stat){
	if(nv==0) return false;
	bool res=false; stat=100;
	double tolerance=.05,sq=0; for(int j=0;j<3;j++) sq+=pow(V[j],2);
	for(int i=0;i<nv;i++){
		double dot=0,sqo=0,diff=0;
		for(int j=0;j<3;j++){
			dot+=Vout[i][j]*V[j]; sqo+=pow(Vout[i][j],2);
			diff+=fabs(Aout[i][j]-A[j]);
		}
		double st=sqrt(MAX(0,1-pow(dot,2)/(sq*sqo)))+diff; stat=MIN(stat,st);
		if(st<tolerance) res=true;
	}
	return res;
}
void resTorus_c::LevMstart(Angles theta,const Position Q,PSPT &QP3,Vector<double,3> &B,
			   Matrix<double,3,3> &A,Matrix<double,3,3> &dQdt,
			   double &dt2dla,double &chisqo){
	Matrix<double,3,3> dJdt;
	double D,tp=from_librationAngle(theta[2],dt2dla,D); theta[2]=tp;
	double dJ2dtp=dDelta(tp,D);
	Actions Jp=resJp; Jp[2]+=D; Jp=fixJ(Jp,theta,tp,dJ2dtp,dJdt);
	Torus T=InterpTorus_2n(Tgrid,np,Jbar_grid,dJ_grid,Jp);
	T.LevCof3DTrue(theta,dJdt,Q,QP3,chisqo,B,dQdt);
	B[2]*=dt2dla;//convert derivs to libration angle
	for(int i=0;i<3;i++) dQdt[i][2]*=dt2dla;
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
void resTorus_c::LevMstep(Angles &theta,const Position Q,PSPT &QP3,Vector<double,3> &B,Matrix<double,3,3> &A,
			  Matrix<double,3,3> &dQdt,double &dt2dla,double &chisqo,double &lam){
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
	LevMstart(t_try,Q,QP3,Btry,Atry,dQdtry,dt2dla,chisq);
	if(chisq<chisqo  && !std::isnan(Btry(0))) {
		lam *= 0.125;
		chisqo = chisq; theta = t_try; tidy_angles(theta);
		A = Atry; B = Btry; dQdt = dQdtry;
	} else {
		lam *= 8.;
	}
	delmatrix(AA,3);
}
int resTorus_c::containsPoint(const Position &Q,Velocity *Vout,Angles *Aout,double *densout){
	const    int    maxit1=50, maxtry=100;
	const    double tiny=1.e-7;

	double rq = Q[0]*Q[0]+Q[1]*Q[1],rtin = sqrt(rq)*tiny; // tolerance in position 
	double dt2dla,chisqo,stat;
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
		LevMstart(theta,Q,QP3,B,A,dQdt,dt2dla,chisqo);
		int it=0;
		double lam=0.5;
		while(chisqo>rtin && maxit1>it++ && lam < 1.e10 ) {  // Lev Mar iteration
			LevMstep(theta,Q,QP3,B,A,dQdt,dt2dla,chisqo,lam);
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
				//fabs((dQdt[0][0]*dQdt[1][1]-dQdt[0][1]*dQdt[1][0])*Q[0]*dt2dla);
				nv++;
			} else nfail++;
		}
		ntry++; nfail++;
	}
	return nv;
}
