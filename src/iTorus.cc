/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
#include "iTorus.h"
#include "Numerics.h"
#include "Random.h"

iTorus::iTorus(Actions Jin,eTorus **Tgridin,int npin,int nrin,Actions Jbarin,Actions dJin){
	Tgrid=Tgridin; np=npin; nr=nrin; Jbar=Jbarin; dJ=dJin;
	setJ(Jin);
}
eTorus iTorus::eT1(Actions Jp) const{
	eTorus T=Interp_eTorus_mn(Tgrid,np,nr,Jbar,dJ,Jp);
	return T;
}
void iTorus::setJ(Actions Jin){//changes J and eT
	J=Jin; eT=Interp_eTorus_mn(Tgrid,np,nr,Jbar,dJ,J);
}
bool iTorus::oldV (Velocity *Vout,Angles *Aout,int nv,const Velocity &V,
		   const Angles &A,double &stat) const{
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
inline void tidy(Angles &theta){
	for(int i=0;i<3;i++){
		if(fabs(theta[i])>100) theta[i] -= TPi*int(theta[i]/TPi);
		while (theta[i]< 0. ) theta[i] += TPi;
		while (theta[i]> TPi) theta[i] -= TPi;
	}
}
void iTorus::LevMstart(Angles theta,const Position Q,PSPT &QP3,Vector<double,3> &B,
			   Matrix<double,3,3> &A,Matrix<double,3,3> &dQdt,double &chio) const{
	Matrix<double,3,3> dJdt;
	eTorus T=fixJ(theta,dJdt);
	T.LevCof3DTrue(theta,dJdt,Q,QP3,chio,B,dQdt);
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
double iTorus::LevMstep(Angles &theta,const Position Q,PSPT &QP3,Vector<double,3> &B,Matrix<double,3,3> &A,
			  Matrix<double,3,3> &dQdt,double &chisqo,double &lam) const{
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
	Angles t_try; double Bmag=0;
	for(int i=0;i<3;i++){
		t_try[i]=theta[i]+BB[i];
		Bmag+=pow(BB[i],2);
	}
	LevMstart(t_try,Q,QP3,Btry,Atry,dQdtry,chisq);
	if(chisq<chisqo  && !std::isnan(Btry(0))) {
		lam *= 0.125;
		chisqo = chisq;	theta = t_try; tidy(theta);
		A = Atry; B = Btry; dQdt = dQdtry;
	} else {
		lam *= 8.;
	}
	delmatrix(AA,3); return Bmag;
}
int iTorus::containsPoint(const Position &Q,Velocity *Vout,Angles *Aout,double *densout,int vmax) const{
	const    int    maxit1=50, maxtry=100;
	const    double tiny=1.e-10;

	double rq = Q[0]*Q[0]+Q[1]*Q[1],rtin = sqrt(rq)*tiny; // tolerance in position
	double chisqo,stat;
	Matrix<double,3,3> A,dQdt,Atry,dQdtry;;
	Vector<double,3> B,Btry;
	Angles theta;
	PSPT QP3;
	int nfail=0,nv=0,ntry=0;
	Random3 ran(7986535);
	while(ntry<maxtry && nv<vmax){
		theta[0] = TPi*(ran.RandomDouble());
		theta[1] = TPi*(ran.RandomDouble());
		theta[2] = TPi*(ran.RandomDouble());
		LevMstart(theta,Q,QP3,B,A,dQdt,chisqo);
		int it=0;
		double lam=0.5,Bmag=1;
		while(maxit1>it++ && lam < 1.e10 && Bmag>1e-12) {  // Lev Mar iteration
			Bmag=LevMstep(theta,Q,QP3,B,A,dQdt,chisqo,lam);
		}
		if(lam<1e-4 && Bmag<1e-11 && chisqo > 1e-8 && chisqo < 10){
			//printf("iTorus::containsPoint no can do: %d %g\n",it,lam);
			return 0;
		}
		if(chisqo<rtin){
			for(int i=0;i<3;i++) Vout[nv][i]=QP3(i+3);
			Aout[nv]=theta;
			if(!oldV(Vout,Aout,nv,Vout[nv],Aout[nv],stat)){
				densout[nv]=fabs(LUDet3(dQdt));
				nv++;
			} else nfail++;
		}
		ntry++;
	}
	return nv;
}
bool iTorus::InOrbit(const Position &Rzphi){
	Velocity Vm[2]; Angles Am[2]; double dens[2];
	if(containsPoint(Rzphi,Vm,Am,dens,1)) return 1; else return 0;
}
double iTorus::refine(const Position &Rzphi,double &J0,double &J2){
	//J0 should always lie outside orbit and J2 in it
	double SMALL=1e-6;
	Actions Jtry=J;
	Jtry[2]=J0; setJ(Jtry);
	if(InOrbit(Rzphi)){
		printf("iTorus::InOrbit J0 doesn't put us outside orbit!\n");
	}
	Jtry[2]=J2; setJ(Jtry);
	if(!InOrbit(Rzphi)){
		printf("iTorus::InOrbit J2 doesn't put us inside orbit!\n");
	}
	while(fabs(J0-J2)>SMALL){
		Jtry[2]=.5*(J0+J2); setJ(Jtry);
		if(!InOrbit(Rzphi)) J0=Jtry[2]; else J2=Jtry[2];
	}
	Jtry[2]=J2; setJ(Jtry);
	//printf("refine: %d %f\n",InOrbit(Rzphi),J2);
	return J2;
}
void iTorus::get_crit_Jp(const Position &Rzphi,double &Japo,double &Jperi,
			 double Jbot,double Jtop){
	//For given Jr,Jz and Jp=Jperi Rzphi is at pericentre of orbit, Jp=Japo & it's
	//at apocentre. When Jperi>Jtop, returns Jtop; when Japo<Jbot
	//returns Jbot
	Actions Jin=J;//actions();
	double Jmid,Jlast;
	int nv=InOrbit(Rzphi);
	if(!nv){
		printf("iTorus::get_crit_Jp called with Jp=%f putting Rzphi outside orbit\n",J[2]);
		exit(0); //return;
	}// We are inside orbit
	Jmid=Jin[2];//so we can restore state
	int k=0;
	while(nv && k<10){//go after Japo
		Jlast=Jin[2]; Jin[2]=.5*(Jin[2]+Jbot); setJ(Jin);
		nv=InOrbit(Rzphi); k++;
	}// now J[2]<Japo<Jlast
	if(!nv) Japo=refine(Rzphi,Jin[2],Jlast);
	else Japo=Jbot;
	Jin[2]=Jmid; setJ(Jin); nv=1; k=0; //restore state
	while(nv && k<10){//go after Jperi
		Jlast=Jin[2]; Jin[2]=.5*(Jin[2]+Jtop); setJ(Jin);
		nv=InOrbit(Rzphi); k++;
	}//Jlast<Jperi<J[2]
	if(!nv) Jperi=refine(Rzphi,Jin[2],Jlast);
	else Jperi=Jtop;
	Jin[2]=Jmid; setJ(Jin);//restore original state
}
eTorus iTorus::fixJ(const Angles& A) const{
	//These methods return the eTorus trough the given angles using
	//the Fourier coefficients of eT. They don't change J or eT
	Actions Jfix=J;
	Frequencies Om=omega();
	int10 iN=i1(),jN=i2(),kN=i3();
	vec10 hN=hn();
	for(int i=1;i<10;i++){
		if(jN(i)) continue;//exclude variation in Jz
		double Om_r=iN(i)*Om(0)+jN(i)*Om(1)+kN(i)*Om(2);
		double cw=cos(iN(i)*A(0)+jN(i)*A(1)+kN(i)*A(2));
		double fac=2*hN(i)*cw/Om_r;
		if(fabs(fac)>0.4) printf("iTorus::fixJ %d %d %f\n",iN(i),kN(i),Om_r);
		Jfix[0]-=iN(i)*fac; Jfix[1]-=jN(i)*fac; Jfix[2]-=kN(i)*fac;
		if(std::isnan(Jfix[0]) || std::isnan(Jfix[2])){
			printf("error in iTorus::fix\n%d (%g,%g) (%d %d) %g %g %g\n",
			       i,Jfix[2],Jfix[0],iN(i),kN(i),hN[i],Om_r,fac);
			printf("(%g %g %g)\n",A(0),A(1),A(2));
		}
	}
	eTorus T=Interp_eTorus_mn(Tgrid,np,nr,Jbar,dJ,Jfix);
	return T;
}
eTorus iTorus::fixJ(const Angles& A,Matrix<double,3,3> &dJdt) const{
	//As above but computes derivatives
	Actions Jfix=J;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) dJdt[i][j]=0;
	Frequencies Om=omega();
	int10 iN=i1(),jN=i2(),kN=i3();
	vec10 hN=hn();
	for(int i=1;i<10;i++){
		if(jN(i)) continue;//exclude variation in Jz
		double Om_r=iN(i)*Om(0)+jN(i)*Om(1)+kN(i)*Om(2);
		double th=iN(i)*A(0)+jN(i)*A(1)+kN(i)*A(2);
		double cw=cos(th),sw=sin(th);
		double fac=2*hN(i)*cw/Om_r,fas=2*hN(i)*sw/Om_r;
		if(fabs(fac)>0.4) printf("iTorus::fixJ %d %d %f\n",iN(i),kN(i),Om_r);
		Jfix[0]-=iN(i)*fac; Jfix[1]-=jN(i)*fac; Jfix[2]-=kN(i)*fac;
		dJdt[0][0]+=iN(i)*iN(i)*fas; dJdt[0][1]+=iN(i)*jN(i)*fas; dJdt[0][2]+=iN(i)*kN(i)*fas;
		dJdt[1][0]+=jN(i)*iN(i)*fas; dJdt[1][1]+=jN(i)*jN(i)*fas; dJdt[1][2]+=jN(i)*kN(i)*fas;
		dJdt[2][0]+=kN(i)*iN(i)*fas; dJdt[2][1]+=kN(i)*jN(i)*fas; dJdt[2][2]+=kN(i)*kN(i)*fas;
		if(std::isnan(Jfix[0]) || std::isnan(Jfix[2])){
			printf("error in iTorus::fix\n%d (%g,%g) (%d %d) %g %g %g\n",
			       i,Jfix[2],Jfix[0],iN(i),kN(i),hN[i],Om_r,fac);
			printf("(%g %g %g)\n",A(0),A(1),A(2));
		}
	}
	eTorus T=Interp_eTorus_mn(Tgrid,np,nr,Jbar,dJ,Jfix);
	return T;
}
GCY iTorus::FullMap(const Angles& A) const{
	eTorus T=fixJ(A);
	GCY gcy=T.FullMap(A);
	return gcy;
}
void iTorus::SOS(ostream& sfile,const int npt){//work round with true angle
	double PiH=.5*Pi,eps=.0001,A_old;
	Angles A;
	A_old=PiH;
	for(int i=0;i<npt;i++){
		A[0]=PiH+i/(double)(npt-1)*(TPi-1.e-6);
		A[1]=0;
		A[2]=A_old;
		GCY gcy=FullMap(A);
		while(gcy[1]<-0.01){
			A[1]+=.001; gcy=FullMap(A);
		}
		while(gcy[1]>0.01){
			A[1]-=.001; gcy=FullMap(A);
		}
		if(gcy[2]>Pi) gcy[2]-=TPi;
		double diff=gcy[2],diff_old=diff,dA=.1;
		int n=0;
		while(diff*diff_old>eps*eps){
			if(n>40)printf("error in iTorus::SOS (%g %g %g)\n",diff,diff_old,dA);
			if(n>45) exit(0); n++;
			if(n>1 && fabs(diff)>fabs(diff_old)){//we're going the wrong way
				dA=-dA; A[2]+=dA; diff=diff_old;//reverse,undo
			} else {
				diff_old=diff; A_old=A[2];
			}
			A[2]+=dA;
			gcy=FullMap(A);	if(gcy[2]>Pi) gcy[2]-=TPi;
			diff=gcy[2];
		}//we straddle root
		while(fabs(diff)>eps && n<20){
			double A_new=A[2];
			A[2]=.5*(A_new+A_old);
			gcy=FullMap(A);	if(gcy[2]>Pi) gcy[2]-=TPi;
			double test=gcy[2];
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
