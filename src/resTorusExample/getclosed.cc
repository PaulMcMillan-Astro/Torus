#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "/u/c/spress/hacked/press.h"
#include "falPot.h"
#include "pstuff.h"

extern Potential *Phi;

void derivs(double t,double *y,double *dydt){
	for(int i=0;i<2;i++){
		dydt[i]=y[i+2];
	}
	double dPR,dPz;
	(*Phi).eff(y[0],y[1],dPR,dPz);
	dydt[2]=-dPR; dydt[3]=-dPz;
}
double Ez(double *y){
	return .5*(pow(y[2],2)+pow(y[3],2))+(*Phi).eff(y[0],y[1]);
}
double vsqx(double E,double x){
	return 2*(E-(*Phi).eff(x,0));
}
double goround(double *par,double x,double &Jz,double &Om){//throws up and determines Rdown
	double pi=acos(-1), vsq=vsqx(par[0],x);
	if(vsq<0){
		printf("%f %f\n",x,vsq); exit(0);
	}
	double x2[4]={x,0,0,sqrt(vsq)};
	double dxdt[4],xscal[4],t=0,hdid,hnext,htry=2.e-2,eps=1.e-10;
	Jz=0;
	for(int i=0; i<2; i++){
		xscal[i]=x; xscal[i+2]=fabs(x2[3]);
	}
	double tlast,Rlast,zlast,pzlast=x2[3];
	while(x2[1]>=0){
		tlast=t; Rlast=x2[0]; zlast=x2[1]; pzlast=x2[3];
		derivs(t,x2,dxdt);
		rkqs(x2,dxdt,4,&t,htry,eps,xscal,&hdid,&hnext,&derivs); htry=hnext;
		Jz+=.5*(x2[3]+pzlast)*(x2[1]-zlast);
	}
	double vbar=.5*(x2[3]+pzlast);
	Jz-=vbar*x2[1]; Jz/=pi;
	Om=pi/(t-x2[1]/vbar);
	return x-(x2[0]-x2[1]/(zlast-x2[1])*(Rlast-x2[0]));//difference between up and down
}
double goround(double *par,double x){//throws up and determines Rdown
	double Jz,Om;
	return goround(par,x,Jz,Om);
}
double sorted(double E,double x0){
	int k=0;
	while(vsqx(E,x0)<0){
		double dPR,dPz; (*Phi).eff(x0,0.,dPR,dPz);
		if(dPR<0) x0*=1.1;
		else x0*=.9;
		k++; if(k>50){ printf("vsq: %f %f\n",x0,E); exit(0);}
	}
	return x0;
}
double get_closed(double x0,double E1,double &Jz,double &Om){
	x0=sorted(E1,x0);// ensure vsq>=0
	double x1,dx0,dx1;
	dx0=goround(&E1,x0);
	double dx=.8; x1=sorted(E1,x0/dx); dx1=goround(&E1,x1);
	int k=0,bo=0;
	while(dx0*dx1>0){//we haven't bracketed the root
		if(fabs(dx1)<fabs(dx0)){//more of same
			x0=x1; dx0=dx1;
		}else{// back off
			if(bo==0){
				dx=1/dx; bo=1;
			}else dx=pow(dx,-0.9);
		}
		x1=sorted(E1,x0*dx); dx1=goround(&E1,x1);
		k++;
		if(k==21){
			printf("get_closed: %f %f %g %g %g\n",x0,x1,dx0,dx1,dx);
			if(fabs(dx0)<1.e-7){
				return -x0;
			}else{
				printf("problem in get_closed(): %f %f %f %f\n",x0,x1,dx0,dx1);
				return -1;
			}
		}
	}
	double R0=zbrent(&E1,&goround,x0,x1,dx0,dx1,1.e-6,25);
	goround(&E1,R0,Jz,Om);
	return R0;
}
