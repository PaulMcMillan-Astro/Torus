#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "pstuff.h"

void four1(float *data, unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
void four1(double *data, unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
void fourn(float *data, unsigned long *nn, int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

#undef SWAP
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;
	u = new double[n];

//	u=(double*)calloc(n,sizeof(double));
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	delete[] u;
}

double splint(double *xa, double *ya, double *y2a, int np, double x)
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splint: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splint\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double splintp(double *xa, double *ya, double *y2a, int np, double x)
// returns first deriv
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splintp: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splintp\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return (ya[khi]-ya[klo])/h+((3.*b*b-1.)*y2a[khi]-(3.*a*a-1.)*y2a[klo])*h/6.0;
}

double splintp2(double *xa, double *ya, double *y2a, int np, double x)
// returns second deriv
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splintp2: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splintp\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return y2a[klo]*a + y2a[khi]*b;
}

//To use these RK routines you call rkqs

void rk4(double y[], double dydx[],int n, double x, double h,
	 void (*derivs)(double, double [], double []))
{
	int i;
	double xh,hh,h6,dym[6],dyt[6],yt[6];
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
		y[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

void rkck(double y[], double dydx[],int N, double x, double h, double yout[],
	  double yerr[], void (*derivs)(double, double [], double []))
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
	b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
	b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
	b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
	b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
	c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
	dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double ak2[6],ak3[6],ak4[6],ak5[6],ak6[6],ytemp[6];

	for (i=0;i<N;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);
	for (i=0;i<N;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<N;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double *y,double *dydx,int N,double *x,double htry,double eps,
	  double *yscal, double *hdid, double *hnext,
	  void (*derivs)(double, double*, double*))
{
	int i;
	double errmax,h,htemp,xnew,yerr[6],ytemp[6];
	h=htry;
	for (;;) {
		rkck(y,dydx,N,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=0;i<N;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x)
			printf("At x= %g stepsize underflow in rkqs\n",xnew);

	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=0;i<N;i++) y[i]=ytemp[i];
}

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

#define EPS 3.0e-8
double zbrent(double *pars,double (*func)(double*,double), double x1, double x2,double fa,double fb,double tol,int itmax){
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		printf("Root must be bracketed in zbrent:\nx1,x2,f1,f2: %g %g %g %g",x1,x2,fa,fb);
		exit(0);}
	fc=fb;
	for (iter=1;iter<=itmax;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a; fc=fa; e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s; q=1.0-s;
			} else {
				q=fa/fc; r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d; d=p/q;
			} else {
				d=xm; e=d;
			}
		} else {
			d=xm; e=d;
		}
		a=b; fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(pars,b);
	}
	if(fabs(xm)>100*tol1){
		printf("Too many iters in zbrent: %g %g ",b,fb);// exit(0);
	}
	return b;
}

#undef EPS
