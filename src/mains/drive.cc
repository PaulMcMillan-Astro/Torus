#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "Torus.h"
//#include "/u/sm/mongo.h"
//#include "/u/c/spress/hacked/press.h"
#include "falPot.h"

#include "pstuff.h"

using namespace std;
using namespace Units;

//#define N  20000 //16384// 8192 //4096
#define NS 16000
#define NV 40
#define NP 140
#define NT 40
//#define SGN(A) ((A)>=0 ? 1:-1)
const int kpl=1,kres=1;
Potential *Phi;
static double PI=acos(-1);

double get_closed(double,double,double&,double&);

void quadratic(double a,double b,double c,double &r1, double &r2){
	double discr=b*b-4*a*c;
	if(discr<0) printf("discr<0 in quadratic: %g %g %g %g\n",a,b,c,discr);
	double q=-.5*(b+SGN(b)*sqrt(discr));
	r1=q/a; r2=c/q;
}
void fitquad(double *Jm,double *hm,double J,double *h){//Interpolates & gets derivs between 3 points
	double dp1m1=Jm[1]-Jm[-1],d0m1=Jm[0]-Jm[-1],dp10=Jm[1]-Jm[0];
	h[2]=-2/(d0m1*dp10)*(hm[0]-(hm[-1]*dp10+hm[1]*d0m1)/dp1m1);
	h[1]=(hm[1]-hm[-1])/dp1m1+(Jm[1]+Jm[-1]-2*J)/
	     (d0m1*dp10)*(hm[0]-(hm[-1]*dp10+hm[1]*d0m1)/dp1m1);
	h[0]=(hm[-1]*(Jm[1]-J)+hm[1]*(J-Jm[-1]))/dp1m1
	     +(J-Jm[-1])*(Jm[1]-J)/(d0m1*dp10)*(hm[0]-(hm[-1]*dp10+hm[1]*d0m1)/dp1m1);
}
void extract(double *data,double **sr,double **st,int nf){
	//pull FFTs of Jr and Jt from data
	int i,j,k,l,n2=nf*nf;
	for(i=1;i<nf/2;i++){//case of 2 non-vanishing frequencies
		for(j=1;j<nf;j++){
			k=nf*i+j; l=nf*nf-(i-1)*nf-j;
			sr[k][0]=.5/n2*(data[2*k]+data[2*l]);
			sr[k][1]=.5/n2*(data[2*k+1]-data[2*l+1]);
			st[k][0]=.5/n2*(data[2*k+1]+data[2*l+1]);
			st[k][1]=-.5/n2*(data[2*k]-data[2*l]);
		}
	}
	for(j=1;j<nf/2;j++){//cases where one frequency vanishes
		k=j; l=nf-j;//line freq0=0
		sr[k][0]=.5/n2*(data[2*k]+data[2*l]);
		sr[k][1]=.5/n2*(data[2*k+1]-data[2*l+1]);
		st[k][0]=.5/n2*(data[2*k+1]+data[2*l+1]);
		st[k][1]=-.5/n2*(data[2*k]-data[2*l]);
		k=nf*j; l=nf*(nf-j);//line freq1=0
		sr[k][0]=.5/n2*(data[2*k]+data[2*l]);
		sr[k][1]=.5/n2*(data[2*k+1]-data[2*l+1]);
		st[k][0]=.5/n2*(data[2*k+1]+data[2*l+1]);
		st[k][1]=-.5/n2*(data[2*k]-data[2*l]);
	}
	sr[0][0]=data[0]/n2; sr[0][1]=0;//zero-freqs
	st[0][0]=data[1]/n2; st[0][1]=0;
}
double H0(Torus *T,Potential *Phi,Angles &thetas){
	GCY w=T->FullMap(thetas);
	return .5*(pow(w[3],2)+pow(w[4],2))+Phi->eff(w[0],w[1]);
}
void getHn(Torus *T,Potential *Phi,double **sr,double **st,int nf){
	Angles thetas; thetas[2]=0;
	double *h = new double[2*nf*nf];
	double dt=2*PI/(double)nf;
	for(int i=0;i<nf;i++)
		for(int j=0;j<nf;j++){ h[2*(i*nf+j)+1]=0;//H is real
	h[2*(i*nf+j)]=1e6;
		}
	for(int i=0;i<=nf/2;i++){
		int i1=nf-i;
		thetas[0]=i*dt;
		for(int j=0;j<=nf/2;j++){
			int j1=nf-j,j2=nf/2-j,j3=nf/2+j;
			thetas[1]=j*dt; double a=H0(T,Phi,thetas);
			h[2*(i*nf+j)]=a;
			if(j3<nf) h[2*(i*nf+j3)]=a;//N-S symmetry
			if(i1!=i && i1<nf){
				if(j1<nf) h[2*(i1*nf+j1)]=a;//time-reverse symmetry
				if(j2<nf) h[2*(i1*nf+j2)]=a;//both symmetries
			}
		}
	}
/*	for(int i=0;i<nf;i++)
		for(int j=0;j<nf;j++){
			if(h[2*(i*nf+j)]>1e5) printf("error %d %d\n",i,j);
		}*/
	unsigned long nn[2]={nf,nf};
	fourn(h-1,nn-1,2,1);
	extract(h,sr,st,nf);
	delete [] h;
}
void old_getHn(Torus *T,Potential *Phi,double **sr,double **st,int nf){
	Angles thetas; thetas[2]=0;
	double *h = new double[2*nf*nf];
	double dt=2*PI/(double)nf;
	for(int i=0;i<nf;i++)
		for(int j=0;j<nf;j++) h[2*(i*nf+j)+1]=0;//H is real
	for(int i=0;i<nf;i++){
		thetas[0]=i*dt;
		for(int j=0;j<nf;j++){
			thetas[1]=j*dt;
			h[2*(i*nf+j)]=H0(T,Phi,thetas);;
		}
	}
	unsigned long nn[2]={nf,nf};
	fourn(h-1,nn-1,2,1);
	extract(h,sr,st,nf);
	delete [] h;
}
/*void plotnet(Torus *T){
	Angles theta(0.);
	PSPD Rz(7,2,0,0);
	float *R = new float[NP*NP];
	float *z = new float[NP*NP];
	int k=0;
	for(int i=0;i<NP;i++){
		theta[0]=i*2*PI/(double)NP;
		for(int j=0;j<NP;j++){
			theta[1]=j*2*PI/(double)NP;
			Rz=T->MapfromToy(theta);
			R[k]=Rz[0]; z[k]=Rz[1]; k++;
		}
	}
	printf("%d points plotted\n",k);
	plots(1,R,z,6,11,-3,3,"R",1,"z",1,-.9,10);
	points(63.05,1,R,z,k);
	grend(1);
	delete[] R; delete[] z;
}*/
Actions changeJr(Actions J,Frequencies Om,double dJr){
	Actions J1=J;
	J1[0]+=dJr; J1[1]-=Om[0]/Om[1]*dJr;
	return J1;
}
Actions changeJz(Actions J,Frequencies Om,double dJz){
	Actions J1=J;
	J1[0]-=Om[1]/Om[0]*dJz; J1[1]+=dJz;
	return J1;
}
void BoxMass(Torus *T,Position Bottom,Position Top,float *vd,int &nvisit){//Bottom and Top edges of box for which mass wanted
	if(Bottom[1]<0) Bottom[1]*=-1;
	if(Top[1]<0) Top[1]*=-1;//get point into upper half plane
	double db1,db2,dt1,dt2;
	Velocity vb1,vb2,vt1,vt2;
	if(1!=T->containsPoint(Bottom,vb1,db1,vb2,db2)) return;//bottom of box not in orbit
	nvisit++;
	while(1!=T->containsPoint(Top,vt1,dt1,vt2,dt2)){
		Top=.5*(Top+Bottom);
	}//now both Top and Bottom in box
	double dbs=db1*db1, dts=dt1*dt1;//squared Jacobians
	double zts=pow(Top[1],2),zbs=pow(Bottom[1],2);
	double fac=(zts-zbs)/(dbs-dts);
	double m1,m2,z0,z0s=((dts*zbs-dbs*zts)/(dts-dbs));
	if(fac<0 || z0s<0){
		printf("In BoxMass z0s, fac: %f %f\n",z0s,fac);
		m1=2*(Top[1]-Bottom[1])/(fabs(dt1)+fabs(db1));
		//m1=0;
	}else{
		z0=sqrt(z0s);
		double zt=MIN(z0,Top[1]);// integrate to either top of box or top of orbit
		//m1=sqrt(z0s-zbs)/fabs(db1)*(asin(zt/z0)-asin(Bottom[1]/z0));
		//m1=sqrt(z0s-zts)/fabs(dt1)*(asin(zt/z0)-asin(Bottom[1]/z0));
		m1=sqrt(fac)*(asin(zt/z0)-asin(Bottom[1]/z0));
	}
	dbs=db2*db2; dts=dt2*dt2;//squared Jacobians
	z0s=((dts*zbs-dbs*zts)/(dts-dbs)); fac=(zbs-zts)/(dts-dbs);
	if(fac<0 || z0s<0){
		printf("In BoxMass z0s, fac: %f %f\n",z0s,fac);
		m2=2*(Top[1]-Bottom[1])/(fabs(dt2)+fabs(db2));
		//m2=0;
	}else{
		z0=sqrt(z0s);
		double zt=MIN(z0,Top[1]);// integrate to either top of box or top of orbit
		//m2=sqrt(z0s-zbs)/fabs(db2)*(asin(zt/z0)-asin(Bottom[1]/z0));
		//m2=sqrt(z0s-zts)/fabs(dt2)*(asin(zt/z0)-asin(Bottom[1]/z0));
		m2=sqrt(fac)*(asin(zt/z0)-asin(Bottom[1]/z0));
	}
	Velocity v1=.5*(vb1+vt1),v2=.5*(vb2+vt2);
	float psi=atan2(v1[1],v1[0]);
	int nv=NV*.5*(psi/PI+1.); if(nv==NV) nv--;
	vd[nv]+=m1;
	psi=atan2(-v1[1],-v1[0]);
	nv=NV*.5*(psi/PI+1.); if(nv==NV) nv--;
	vd[nv]+=m1;
	psi=atan2(v2[1],v2[0]);
	nv=NV*.5*(psi/PI+1.); if(nv==NV) nv--;
	vd[nv]+=m2;
	psi=atan2(-v2[1],-v2[0]);
	nv=NV*.5*(psi/PI+1.); if(nv==NV) nv--;
	vd[nv]+=m2;
}
void find_hres(Torus *T,Potential *Phi,int nf,double *hres,int ifp){// Get resid H on resonant torus
	double **sr,**st; sr=dmatrix(nf*nf,2); st=dmatrix(nf*nf,2);
	getHn(T,Phi,sr,st,nf);
	if(ifp) printf("h0, h2, h4, h6: %f %f %f %f\n",sr[0][0],
	       sr[2*nf+nf-2][0]/pow(Units::kms,2),sr[4*nf+nf-4][0]/pow(Units::kms,2)
	       ,sr[6*nf+nf-6][0]/pow(Units::kms,2));
	hres[0]=sr[0][0]; hres[1]=sr[2*nf+nf-2][0]; hres[2]=sr[4*nf+nf-4][0];
	delmatrix(sr,2); delmatrix(st,2);
}
double plsos(Torus *T){
	float Rm[NS],zm[NS],Rs[NS],pRs[NS];
	ofstream to; to.open("SOS.dat");
	T->SOS(to);
	to.close();
	FILE* ifile; ifile=fopen("SOS.dat","r");
	int i=0; double Rmin=1e3;
	while(!feof(ifile) && i<NS){
		if(5!=fscanf(ifile,"%f %f %f %f %f",Rs+i,zm+i,pRs+i,Rm+i,zm+(NS-1-i))) break;
		Rmin=MIN(Rs[i],Rmin);
		pRs[i]/=Units::kms; i++;
	}
	if(i==NS) printf("NS too small in plSOS\n");
	fclose(ifile);
//	if(kpl==1) connect(Rs,pRs,i);
	for(int j=0;j<i;j++) pRs[j]=-pRs[j];
//	if(kpl==1) connect(Rs,pRs,i);
	return Rmin;
}
/* void define_region2(Torus *Tgrid,double *Jrgrid,int ng,Actions &Jres,Frequencies &Omres,
		    Potential *Phi,double tolJ,int nf,double &omp,double *h){
	resTorus rT(Tgrid,Jrgrid,ng,Jres,omp,h); rT.setI(rT.Imin);
	Actions J,J1,J2; rT.getJs(J1,J2);
	if(ng%2==0){ printf("ng must be odd\n"); return;}
	int j;
	Torus T; Frequencies Om;
	J=J1;
	double Omgrd[ng],dOmgrd[ng],y2[ng];
	double **hresgrd; hresgrd=dmatrix(ng,3);
	T.AutoFit(J,Phi,tolJ); Om=T.omega();
	for(int i=0;i<ng;i+=2){
//		double f=.3+.4*i/(float)(ng-1);//from 0.2 to 0.7
		double f=i/(float)(ng-1);
		Jrgrid[i]=(1-f)*J1[0]+f*J2[0];
		J=changeJr(J,Om,Jrgrid[i]-J[0]);
		Tgrid[i].AutoFit(J,Phi,tolJ);
		printf("PT: %d ",(Tgrid[i].canmap()).NumberofParameters());
		Om=Tgrid[i].omega(); Omgrd[i]=Om(0); dOmgrd[i]=Om(0)-Om(1);
		int ifp=1;// if(i==ng/2) ifp=1; else ifp=0;
		find_hres(&Tgrid[i],Phi,nf,hresgrd[i],ifp);
	}
	J=J1;
	for(int i=1;i<ng;i+=2){
		double f=i/(float)(ng-1);
		Jrgrid[i]=(1-f)*J1[0]+f*J2[0];
		J=changeJr(J,Om,Jrgrid[i]-J[0]);
		Tgrid[i]=InterpTorus2(Tgrid,Jrgrid,ng,Jrgrid[i]);
		Om=Tgrid[i].omega(); Omgrd[i]=Om(0); dOmgrd[i]=Om(0)-Om(1);
		int ifp=1;// if(i==ng/2) ifp=1; else ifp=0;
		find_hres(&Tgrid[i],Phi,nf,hresgrd[i],ifp);
	}
	double x[ng],y[ng],z[ng];
	for(int i=0;i<ng;i++){
		x[i]=Jrgrid[i]/Units::kms; y[i]=hresgrd[i][1]/pow(Units::kms,2);
		z[i]=20*hresgrd[i][2]/pow(Units::kms,2);
	}
	connect(x,y,ng); setltype(2);
	connect(x,z,ng); setltype(0);
	double Jrres=Jres[0];
	spline(Jrgrid,dOmgrd,ng,2e30,2e30,y2);
	double om=splint(Jrgrid,dOmgrd,y2,ng,Jrres);
	omp=splintp(Jrgrid,dOmgrd,y2,ng,Jrres);
	while(fabs(om)>1e-5*Omres[0]){//Refine Jrres by Newton-Raphson
		Jrres-=om/omp;
		om=splint(Jrgrid,dOmgrd,y2,ng,Jrres);
		omp=splintp(Jrgrid,dOmgrd,y2,ng,Jrres);
	}
	Jres=changeJr(Jres,Omres,Jrres-Jres[0]);
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][0];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	double om2=splintp(Jrgrid,y,y2,ng,Jrres);
	double omp2=splintp2(Jrgrid,y,y2,ng,Jrres);
	printf("Jrres, Omega, G: %f %f %f %f\n",Jrres,om2,omp,omp2); //omp=omp2;
	spline(Jrgrid,Omgrd,ng,2e30,2e30,y2);
	Omres[0]=splint(Jrgrid,Omgrd,y2,ng,Jrres);
	Omres[1]=Omres[0];
//value, gradient and 2nd deriv of h2
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][1];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	h[0]=splint(Jrgrid,y,y2,ng,Jrres);
	h[1]=splintp(Jrgrid,y,y2,ng,Jrres);
	h[2]=splintp2(Jrgrid,y,y2,ng,Jrres);
	double dJ=Jrgrid[1]-Jrgrid[0],x0=Jrres-dJ,y0=h[0]-dJ*h[1],x1=Jrres+dJ,y1=h[0]+dJ*h[1];
//	relocate(x0/Units::kms,y0/pow(Units::kms,2));
//	draw(x1/Units::kms,y1/pow(Units::kms,2));
//value, gradient and 2nd deriv of h4
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][2];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	h[3]=splint(Jrgrid,y,y2,ng,Jrres);
	h[4]=splintp(Jrgrid,y,y2,ng,Jrres);
	h[5]=splintp2(Jrgrid,y,y2,ng,Jrres);
	printf("h: %f %f %f %f %f %f\n",h[0]/pow(Units::kms,2),h[1]/Units::kms,h[2],
	       h[3]/pow(Units::kms,2),h[4]/Units::kms,h[5]);
	//h[3]=h[4]=h[5]=0;
	rT.reset(Jres,omp); rT.setI(rT.Imin);
	rT.getJs(J1,J2);
	printf("Refined boundaries: %f %f\n",J1[0]/Units::kms,J2[0]/Units::kms);
	delmatrix(hresgrd,ng);
}*/
void define_region1(Torus *Tgrid,double *Jrgrid,int ng,Actions &Jres,Frequencies &Omres,
		    Potential *Phi,double tolJ,int nf,double &omp,double *h){
	resTorus rT(Tgrid,Jrgrid,ng,Jres,omp,h); rT.setI(rT.Imin);
	Actions J,J1,J2; rT.getJs(J1,J2);
	int j;
	Torus T; Frequencies Om;
	J=J1;
	double Omgrd[ng],dOmgrd[ng],y2[ng];
	double **hresgrd; hresgrd=dmatrix(ng,3);
	T.AutoFit(J,Phi,tolJ); Om=T.omega();
	for(int i=0;i<ng;i++){
//		double f=.3+.4*i/(float)(ng-1);//from 0.2 to 0.7
		double f=i/(float)(ng-1);
		Jrgrid[i]=(1-f)*J1[0]+f*J2[0];
		J=changeJr(J,Om,Jrgrid[i]-J[0]);
		Tgrid[i].AutoFit(J,Phi,tolJ);
		printf("PT: %d ",(Tgrid[i].canmap()).NumberofParameters());
		Om=Tgrid[i].omega(); Omgrd[i]=Om(0); dOmgrd[i]=Om(0)-Om(1);
		int ifp=1;// if(i==ng/2) ifp=1; else ifp=0;
		find_hres(&Tgrid[i],Phi,nf,hresgrd[i],ifp);
	}
	double x[ng],y[ng],z[ng];
	for(int i=0;i<ng;i++){
		x[i]=Jrgrid[i]/Units::kms; y[i]=hresgrd[i][1]/pow(Units::kms,2);
		z[i]=20*hresgrd[i][2]/pow(Units::kms,2);
	}
//	connect(x,y,ng); setltype(2);
//	connect(x,z,ng); setltype(0);
	double Jrres=Jres[0];
	spline(Jrgrid,dOmgrd,ng,2e30,2e30,y2);
	double om=splint(Jrgrid,dOmgrd,y2,ng,Jrres);
	omp=splintp(Jrgrid,dOmgrd,y2,ng,Jrres);
	while(fabs(om)>1e-5*Omres[0]){//Refine Jrres by Newton-Raphson
		Jrres-=om/omp;
		om=splint(Jrgrid,dOmgrd,y2,ng,Jrres);
		omp=splintp(Jrgrid,dOmgrd,y2,ng,Jrres);
	}
	Jres=changeJr(Jres,Omres,Jrres-Jres[0]);
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][0];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	double om2=splintp(Jrgrid,y,y2,ng,Jrres);
	double omp2=splintp2(Jrgrid,y,y2,ng,Jrres);
	printf("Jrres, Omega, G: %f %f %f %f\n",Jrres,om2,omp,omp2); //omp=omp2;
	spline(Jrgrid,Omgrd,ng,2e30,2e30,y2);
	Omres[0]=splint(Jrgrid,Omgrd,y2,ng,Jrres);
	Omres[1]=Omres[0];
//value, gradient and 2nd deriv of h2
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][1];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	h[0]=splint(Jrgrid,y,y2,ng,Jrres);
	h[1]=splintp(Jrgrid,y,y2,ng,Jrres);
	h[2]=splintp2(Jrgrid,y,y2,ng,Jrres);
	double dJ=Jrgrid[1]-Jrgrid[0],x0=Jrres-dJ,y0=h[0]-dJ*h[1],x1=Jrres+dJ,y1=h[0]+dJ*h[1];
//	relocate(x0/Units::kms,y0/pow(Units::kms,2));
//	draw(x1/Units::kms,y1/pow(Units::kms,2));
//value, gradient and 2nd deriv of h4
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][2];
	spline(Jrgrid,y,ng,2e30,2e30,y2);
	h[3]=splint(Jrgrid,y,y2,ng,Jrres);
	h[4]=splintp(Jrgrid,y,y2,ng,Jrres);
	h[5]=splintp2(Jrgrid,y,y2,ng,Jrres);
	printf("h: %f %f %f %f %f %f\n",h[0]/pow(Units::kms,2),h[1]/Units::kms,h[2],
	       h[3]/pow(Units::kms,2),h[4]/Units::kms,h[5]);
	//h[4]=h[5]=0;
	rT.reset(Jres,omp); rT.setI(rT.Imin);
	rT.getJs(J1,J2);
	printf("Refined boundaries: %f %f\n",J1[0]/Units::kms,J2[0]/Units::kms);
	delmatrix(hresgrd,ng);
}
void define_region0(Torus *Tgrid,double *Jrgrid,int ng0,Actions &Jres,Frequencies &Omres,
		   Potential *Phi,double tolJ,int nf,double &omp,double *h){
	resTorus rT(Tgrid,Jrgrid,ng0,Jres,omp,h); rT.setI(rT.Imin);
	Actions J,J1,J2; rT.getJs(J1,J2);
	Jrgrid[0]=J1[0]; Jrgrid[1]=J2[0];
	Tgrid[0].AutoFit(J1,Phi,tolJ); Tgrid[1].AutoFit(J2,Phi,tolJ);
	int ng=9;
	Torus T; Frequencies Om;
	J=J1; Om=Tgrid[0].omega();
	double Jrgrd[ng],Omgrd[ng],dOmgrd[ng],y2[ng];
	double **hresgrd; hresgrd=dmatrix(ng,3);
	for(int i=0;i<ng;i++){
		double f=.01+.98*i/(float)(ng-1);//from 0.2 to 0.8
		Jrgrd[i]=(1-f)*Jrgrid[0]+f*Jrgrid[1];
		T=InterpTorus(Tgrid,Jrgrid,2,Jrgrd[i]);
		Om=T.omega(); Omgrd[i]=Om(0); dOmgrd[i]=Om(0)-Om(1);
		int ifp=1;// if(i==ng/2) ifp=1; else ifp=0;
		find_hres(&T,Phi,nf,hresgrd[i],ifp);
	}
	double x[ng],y[ng],z[ng];
	for(int i=0;i<ng;i++){
		x[i]=Jrgrd[i]/Units::kms; y[i]=hresgrd[i][1]/pow(Units::kms,2);
		z[i]=20*hresgrd[i][2]/pow(Units::kms,2);
	}
//	connect(x,y,ng); setltype(2);
//	connect(x,z,ng); setltype(0);
	double Jrres=Jres[0];
	spline(Jrgrd,dOmgrd,ng,2e30,2e30,y2);
	double om=splint(Jrgrd,dOmgrd,y2,ng,Jrres);
	omp=splintp(Jrgrd,dOmgrd,y2,ng,Jrres);
	printf("omp: %f\n",omp);
	while(fabs(om)>1e-5*Omres[0]){//Refine Jrres by Newton-Raphson
		Jrres-=om/omp;
//		Jrres=0.036473;
		om=splint(Jrgrd,dOmgrd,y2,ng,Jrres);
		omp=splintp(Jrgrd,dOmgrd,y2,ng,Jrres);
//		break;
	}
	Jres=changeJr(Jres,Omres,Jrres-Jres[0]);
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][0];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	double omp2=splintp2(Jrgrd,y,y2,ng,Jrres);
	printf("Jrres, G: %f %f %f\n",Jrres,omp,omp2); //omp=omp2;
	spline(Jrgrd,Omgrd,ng,2e30,2e30,y2);
	Omres[0]=splint(Jrgrd,Omgrd,y2,ng,Jrres);
	Omres[1]=Omres[0];
//value, gradient and 2nd deriv of h2
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][1];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	h[0]=splint(Jrgrd,y,y2,ng,Jrres);
	h[1]=splintp(Jrgrd,y,y2,ng,Jrres);
	h[2]=splintp2(Jrgrd,y,y2,ng,Jrres);
	double dJ=Jrgrd[2]-Jrgrd[0],x0=Jrres-dJ,y0=h[0]-dJ*h[1],x1=Jrres+dJ,y1=h[0]+dJ*h[1];
//	relocate(x0/Units::kms,y0/pow(Units::kms,2));
//	draw(x1/Units::kms,y1/pow(Units::kms,2));
//value, gradient and 2nd deriv of h4
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][2];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	h[3]=splint(Jrgrd,y,y2,ng,Jrres);
	h[4]=splintp(Jrgrd,y,y2,ng,Jrres);
	h[5]=splintp2(Jrgrd,y,y2,ng,Jrres);
	printf("h: %f %f %f %f %f %f\n",h[0]/pow(Units::kms,2),h[1]/Units::kms,h[2],
	       h[3]/pow(Units::kms,2),h[4]/Units::kms,h[5]);
	//h[3]=h[4]=h[5]=0;
	rT.reset(Jres,omp); rT.setI(rT.Imin);
	rT.getJs(J1,J2);
	printf("Refined boundaries: %f %f\n",J1[0]/Units::kms,J2[0]/Units::kms);
	delmatrix(hresgrd,ng);
}
bool define_region(Torus *Tgrid,double *Jrgrid,int ng0,Actions &Jres,Frequencies &Omres,
		   Potential *Phi,double tolJ,int nf,double &omp,double *h){
	resTorus rT(Tgrid,Jrgrid,ng0,Jres,omp,h); rT.setI(rT.Imin);
	Actions J,J1,J2; rT.getJs(J1,J2);
	if(fabs(J1[0]-J2[0])<.000001) {
		printf("Js %f %f\n",J1[0],J2[0]); return false;
	}
	Jrgrid[0]=J1[0]; Jrgrid[1]=J2[0];
	Tgrid[0].AutoFit(J1,Phi,tolJ); Tgrid[1].AutoFit(J2,Phi,tolJ);
	int ng=9;
	Torus T; Frequencies Om;
	J=J1; Om=Tgrid[0].omega();
	double Jrgrd[ng],Omgrd[ng],dOmgrd[ng],y[ng],y2[ng];
	double **hresgrd; hresgrd=dmatrix(ng,3);
	for(int i=0;i<ng;i++){
		double f=.3+.4*i/(float)(ng-1);//from 0.2 to 0.8
		Jrgrd[i]=(1-f)*Jrgrid[0]+f*Jrgrid[1];
		T=InterpTorus(Tgrid,Jrgrid,2,Jrgrd[i]);
		Om=T.omega(); Omgrd[i]=Om(0); dOmgrd[i]=Om(0)-Om(1);
		int ifp=1;// if(i==ng/2) ifp=1; else ifp=0;
		find_hres(&T,Phi,nf,hresgrd[i],ifp);
	}
	double x[ng],z[ng];
	for(int i=0;i<ng;i++){
		x[i]=Jrgrd[i]/Units::kms; y[i]=hresgrd[i][1]/pow(Units::kms,2);
		z[i]=hresgrd[i][2]/pow(Units::kms,2);
	}
//	connect(x,y,ng); setltype(2);
//	connect(x,z,ng); setltype(0);
	double Jrres=Jres[0];
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][0];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	double om=splintp(Jrgrd,y,y2,ng,Jrres);
	omp=splintp2(Jrgrd,y,y2,ng,Jrres);
	while(fabs(om)>1e-5*Omres[0]){//Refine Jrres by Newton-Raphson
		Jrres-=om/omp;
		om=splintp(Jrgrd,y,y2,ng,Jrres);
		omp=splintp2(Jrgrd,y,y2,ng,Jrres);
	}
	Jres=changeJr(Jres,Omres,Jrres-Jres[0]);
	printf("Jrres, G: %f %f\n",Jrres,omp);
	spline(Jrgrd,Omgrd,ng,2e30,2e30,y2);
	Omres[0]=splint(Jrgrd,Omgrd,y2,ng,Jrres);
	Omres[1]=Omres[0];
//value, gradient and 2nd deriv of h2
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][1];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	h[0]=splint(Jrgrd,y,y2,ng,Jrres);
	h[1]=splintp(Jrgrd,y,y2,ng,Jrres);
	h[2]=splintp2(Jrgrd,y,y2,ng,Jrres);
	double dJ=Jrgrd[2]-Jrgrd[0],x0=Jrres-dJ,y0=h[0]-dJ*h[1],x1=Jrres+dJ,y1=h[0]+dJ*h[1];
//	relocate(x0/Units::kms,y0/pow(Units::kms,2));
//	draw(x1/Units::kms,y1/pow(Units::kms,2));
//value, gradient and 2nd deriv of h4
	for(int i=0;i<ng;i++) y[i]=hresgrd[i][2];
	spline(Jrgrd,y,ng,2e30,2e30,y2);
	h[3]=splint(Jrgrd,y,y2,ng,Jrres);
	h[4]=splintp(Jrgrd,y,y2,ng,Jrres);
	h[5]=splintp2(Jrgrd,y,y2,ng,Jrres);
	printf("h: %f %f %f %f %f %f\n",h[0]/pow(Units::kms,2),h[1]/Units::kms,h[2],
	       h[3]/pow(Units::kms,2),h[4]/Units::kms,h[5]);
	rT.reset(Jres,omp); rT.setI(rT.Imin);
	rT.getJs(J1,J2);
	printf("Refined boundaries: %f %f\n",J1[0]/Units::kms,J2[0]/Units::kms);
	delmatrix(hresgrd,ng);
	return true;
}
bool setResTor(double E0,double Lz,Potential *Phi,double &Jzmax,double tolJ,
	      Torus *Tgrid,double *Jrgrid,int ng,Actions &Jres,double &omp,double *h){
/* Given E,Lz we locate the 1:1 resonant orbit Jres and determine the
 * numbers omres', h, h' and h" needed to do p-theory around it.
 * These #s are first determined from tori TM fits in the trapping
 * region with tolJ and then a second time from tori interpolated from tori TM fits
 * with tolJ/2 eithere side of trapping region. These tori are returned
 * in Tgrid. Returns false if no resonance at this (E,Lz), true otherwise */
	//1. Find J at  given E
	//printf("computing Jzmax..");
	double Omz,R=get_closed(6,E0,Jzmax,Omz);
	printf("R Jzmax: %f %f ",R,Jzmax);
	Actions J; J[0]=0.001;  J[1]=Jzmax; J[2]=Lz;
	Torus T; T.AutoFit(J,Phi,tolJ);
	Frequencies Om=T.omega(); double dJr,dJz,dOm=Om(0)-Om(1);
	printf("Om: %f %f %f\n",Om[1],Omz,dOm);
	if(dOm<0) return false;
	//2. Move on E=const to where Omr=Omz
	double Jh=0,Jl=0; int repmax=50, rep=0;
	while(Jh==0 || Jl==0){
		if(dOm>0){//need to decrease Jz
			Jh=J[1]; dJz=-.1*J[1]; J=changeJz(J,Om,dJz);
			if(J[0]<=0){
				printf("Jr gone to zero\n"); return false;
			}
		}else{
			Jl=J[1]; dJz=.1*J[1]; J=changeJz(J,Om,dJz);
			if(J[1]<1e-4){
				printf("Jz going to zero\n"); return false;
			}
		}
		T.AutoFit(J,Phi,tolJ);
		Om=T.omega(); dOm=Om(0)-Om(1);
		if(dOm>0) Jh=J[1]; else Jl=J[1];
		rep++; if(rep>repmax){
			printf("Can't locate resonance\n");
			break;
		}
	}
	rep=0;
	while(fabs(dOm)>2.e-4){
		double Jtry=.5*(Jh+Jl); dJz=Jtry-J[1]; J=changeJz(J,Om,dJz);
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); dOm=Om(0)-Om(1);
		if(dOm>0) Jh=J[1]; else Jl=J[1];
		if(rep>repmax){
			printf("Can't close Omega gap\n");
			break;
		}
		rep++;
		if(rep>repmax-20){
			printf("%f %f %f %f %f %f\n",Jl,Jh,Jtry,Om(0),Om(1),dOm);
		}
	}// Now we have located the resonant torus
	Jres=J; Frequencies Omres=Om;
	int nf=pow(2,7); find_hres(&T,Phi,nf,h,0);
	dJr=.05*J[0]; J=changeJr(Jres,Om,dJr); Jrgrid[0]=J[0];
	Tgrid[0].AutoFit(J,Phi,tolJ);
	Om=Tgrid[0].omega(); dOm=Om(0)-Om(1);
	J=changeJr(Jres,Om,-dJr); Jrgrid[1]=J[0];
	Tgrid[1].AutoFit(J,Phi,tolJ); Om=Tgrid[1].omega();
	omp=.5*(dOm-(Om(0)-Om(1)))/dJr;
	double Delta=2*sqrt(2*h[1]/fabs(omp));
	printf("Jres Omp: %f %f %f\n",Jres[0],omp,h[1]);
	double Jrb=(Jres[0]-Delta)/Units::kms,Jrt=(Jres[0]+Delta)/Units::kms;
	//printf("J,omp: %f %f %f %f\n",Jres[0],Jres[1],omp,Delta);
//Now we have estimate of width of trapped region get good tori on each
//side and from interpolated tori get dependence of h etc on Jr
	h[0]=h[1]; for(int i=1;i<6;i++) h[i]=0;
	if(define_region(Tgrid,Jrgrid,ng,Jres,Omres,Phi,tolJ,nf,omp,h)) return true;
	else return false;
}

int setResTor1(double E0,double Lz,Potential *Phi,double tolJ,
	      Torus *Tgrid,double *Jrgrid,int ng,Actions &Jres,double &omp,double *h){
/* Given E,Lz we locate the 1:1 resonant orbit Jres and determine the
 * numbers omres', h, h' and h" needed to do p-theory around it.
 * These #s are first determined from tori TM fits in the trapping
 * region with tolJ and then a second time from tori interpolated from tori TM fits
 * with tolJ/2 eithere side of trapping region. These tori are returned
 * in Tgrid. Returns -1 if no resonance at this (E,Lz), 0 otherwise */
	//1. Find J at  given E
	Actions J; J[0]=.05; J[1]=.05; J[2]=Lz;
	Torus T; T.AutoFit(J,Phi,tolJ);
	double dE,Jl=0,Jh=0;
	while(Jh==0 || Jl==0){//find bracketing Jrs
		dE=T.energy()-E0;
		if(dE>0){
			Jh=J[0]; if(Jl==0) J[0]*=.8;
		}
		else{
			Jl=J[0]; if(Jh==0) J[0]*=1.2;
		}
		T.AutoFit(J,Phi,tolJ);
	}
	while(fabs(Jh-Jl)>.001){//Now E0 between Jl & Jh
		J[0]=.5*(Jh+Jl);
		T.AutoFit(J,Phi,tolJ);
		if(T.energy()>E0) Jh=J[0];
		else Jl=J[0];
	}
	//2. Move on E=const to where Omr=Omz
	Frequencies Om=T.omega(); double dJr,dOm=Om(0)-Om(1);
	Jh=Jl=0;
	while(Jh==0 || Jl==0){
		if(dOm>0){//need to increase Jr
			Jl=J[0]; dJr=.1*J[0]; J=changeJr(J,Om,dJr);
			if(J[1]<=0){
				printf("Jz gone to zero\n"); return-1;
			}
		}else{
			Jh=J[0]; dJr=-.1*J[0]; J=changeJr(J,Om,dJr);
			if(J[0]<1e-4){
				printf("Jr going to zero\n"); return -1;
			}
		}
		T.AutoFit(J,Phi,tolJ);
		Om=T.omega(); dOm=Om(0)-Om(1);
		if(dOm>0) Jl=J[0]; else Jh=J[0];
	}
	while(fabs(dOm)>1.e-4){
		double Jtry=.5*(Jh+Jl); dJr=Jtry-J[0]; J=changeJr(J,Om,dJr);
		T.AutoFit(J,Phi,tolJ); Om=T.omega(); dOm=Om(0)-Om(1);
		if(dOm>0) Jl=J[0]; else Jh=J[0];
	}// Now we have located the resonant torus
	printf("TM determines Jres=(%f %f %f)\n",
	       J[0]/Units::kms,J[1]/Units::kms,J[2]/Units::kms);
	Jres=J; Frequencies Omres=Om;
	int nf=pow(2,7); find_hres(&T,Phi,nf,h,1);
	dJr=.05*J[0]; J=changeJr(Jres,Om,dJr);
	T.AutoFit(J,Phi,tolJ);
	J=changeJr(Jres,Om,-dJr);
	Om=T.omega(); dOm=Om(0)-Om(1);
	T.AutoFit(J,Phi,tolJ); Om=T.omega();
	omp=.5*(dOm-(Om(0)-Om(1)))/dJr;
	double Delta=2*sqrt(2*h[1]/fabs(omp));
	double Jrb=(Jres[0]-Delta)/Units::kms,Jrt=(Jres[0]+Delta)/Units::kms;
	printf("First estimate of trapped zone (%f %f) kpc km/s\n",Jrb,Jrt);
//Now we have estimate of width of trapped region get good tori on each
//side and from interpolated tori get dependence of h etc on Jr
//	plots(1,Jrgrid,Jrgrid,30,42,0,2.5,"Jr/kpc km s\\u-\\u1",17,
//	      "h/(km s\\u-\\u1)\\u2",17,-.9,10);
	h[0]=h[1]; for(int i=1;i<6;i++) h[i]=0;
	if(ng==2) define_region(Tgrid,Jrgrid,ng,Jres,Omres,Phi,tolJ,nf,omp,h);
	else define_region1(Tgrid,Jrgrid,ng,Jres,Omres,Phi,tolJ,nf,omp,h);
//	setcolour("red");
//Now determine the boundaries J1,J2 of the resonant zone more precisely
	if(ng==2) define_region(Tgrid,Jrgrid,ng,Jres,Omres,Phi,tolJ/3,nf,omp,h);
	else define_region1(Tgrid,Jrgrid,ng,Jres,Omres,Phi,tolJ,nf,omp,h);
//	setcolour("blue");
//	grend(1);
	printf("E E: %f %f\n",Tgrid[0].energy(),Tgrid[1].energy());
	return 0;
}
void plresSoS(ostream& sfile,Torus *Tgrid,double *Jgrid,int ng,Actions Jb,double omp,double *hres){
/* Given in Tgrid[ng] tori on both boundaries of trapped zone and the
 * parameters for p-theory, we plot trapped orbits in the SoS */
	Actions J,Jp; Jp[0]=Jb[0]; Jp[1]=Jb[1]+Jb[0]; Jp[2]=Jb[2];
	sfile << "201" << '\n'; Tgrid[0].SOS(sfile,200);
	sfile << "201" << '\n'; Tgrid[ng-1].SOS(sfile,200);
	sfile << "201" << '\n';
	Torus T; T=InterpTorus(Tgrid,Jgrid,ng,Jb[0]); T.SOS(sfile,200);
	resTorus rT(Tgrid,Jgrid,ng,Jb,omp,hres);
	FILE *ofile; ofile=fopen("res_it.dat","w");
	for(int k=0;k<6;k++){
		double I=.99*(rT.Imin+k*(rT.Imax-rT.Imin)/5.);
		rT.setI(I);
		sfile << 2*NP << '\n';
		rT.SOS(sfile,2*NP);
		printf("libJ: %f\n",rT.librationAction()/Units::kms);
		Angles A; A[0]=PI/2; A[1]=0; A[2]=0; GCY gcy=rT.FullMap(A);
		fprintf(ofile,"%f %f %f %f %f\n",gcy[0],gcy[1],gcy[3],gcy[4],gcy[5]);
		double rOmega=rT.librationOmega();
		Frequencies Om=T.omega();
		if(k==1){
			printf("Now plotting orbit\n");
			FILE *ofile2; ofile2=fopen("res_orb.dat","w");
			FILE *ofile3; ofile3=fopen("res_SOS3.dat","w");
			fprintf(ofile2,"%d\n",k);
			Angles A; float R[NS],z[NS];
			//plots(1,&I,&I,3,8,-3,3,"R/kpc",5,"z/kpc",5,-.9,10);
//			plots(1,&I,&I,6,11,-3,3,"R/kpc",5,"z/kpc",5,-.9,10);
			double pRlast;
			for(int i=0;i<NS;i++){
				double t=2*PI/rOmega*i/(double)NS;
				A[0]=rOmega*t; A[1]=Om[1]*t; A[2]=Om[2]*t;
				GCY gcy=rT.FullMap(A);
				R[i]=gcy[0]; z[i]=gcy[1];
				if(i>0 && z[i-1]*z[i]<0 && gcy[4]>0){
					double fac=z[i]/(z[i]-z[i-1]);
					double Rs=R[i]+fac*(R[i-1]-R[i]);
					double pRs=gcy[3]+fac*(pRlast-gcy[3]);
					fprintf(ofile3,"%f %f\n",Rs,pRs);
				}
				pRlast=gcy[3];
				if(i%4==0) fprintf(ofile2,"%f %f %f\n",t,R[i],z[i]);
			}
//			setcolour("red"); connect(R,z,NS); grend(1);
			fclose(ofile2); fclose(ofile3);
			Position Rzphi; Rzphi[2]=0;
			for(int iR=0;iR<3;iR++){
				Rzphi[0]=4+(7-4)*iR/2.;
				for(int iz=0;iz<3;iz++){
					Rzphi[1]=iz/2.*1.8;
					Angles *Am=new Angles[5]; double tm[5];
					Velocity *Vm=new Velocity[5];
					//int n=rT.containsPoint_Ang(Rzphi,Am,tm);
					int n=rT.containsPoint(Rzphi,Vm,Am,tm);
					if(!n) printf("(%f,%f) not found\n",Rzphi[0],Rzphi[1]);
					else{
						gcy; printf("containsPoint: (%f %f)\n",Rzphi[0],Rzphi[1]);
						for(int i=0;i<n;i++){
							gcy=rT.FullMap(Am[i]);
							printf("    (%f %f)(%f %f %f)\n",gcy[0],gcy[1],Vm[i][0],Vm[i][1],Vm[i][2]);
						}
					}
				}

			}
		}
	}
}
int main(void){
	double R=8,dPR,dPz;
//	double R=4.5,dPR,dPz;
	ifstream Pfile;
	Pfile.open("pot/PJM11_best.Tpot");
	if(!Pfile.is_open()) printf("can't open Tpot file\n");
	Phi = new GalaxyPotential(Pfile);
	Pfile.close();
	double P=(*Phi)(R,0,dPR,dPz);
	double Vc=sqrt(R*dPR);
	double E=.5*Vc*Vc+P, Lz=R*Vc, pih=.5*PI, Jzmax=0;
	Phi->set_Lz(Lz);
	double dv=.6*Vc;//.395*Vc;//.7*Vc;//.32*Vc;
	printf("Vc, dv: %f %f\n",Vc/Units::kms,dv/Units::kms);
	E=.5*dv*dv+(*Phi).eff(R,0); printf("E=%f\n",E);
	//R=8.35;//for untrapped .32Vc
	//R=8.56;//for untrapped .4Vc
	//R=9.46;//for untrapped .6Vc
	//R=4.7;
	int ng=2;
	Actions Jb; double omp,h[6],Jrgrid[ng];
	Torus *Tgrid = new Torus[ng];
	if(setResTor(E,Lz,Phi,Jzmax,.003,Tgrid,Jrgrid,ng,Jb,omp,h)){
	printf("Jrres, G, h: %10.6f %10.6f %10.6f %10.6f %10.6f\n",
	       Jb[0]/Units::kms,omp,h[0]/pow(Units::kms,2),h[1]/Units::kms,h[2]);
	ofstream sfile; sfile.open("resSOS1.dat");
	plresSoS(sfile,Tgrid,Jrgrid,ng,Jb,omp,h);
	sfile.close();
	} else printf("resonance not found\n");
}
