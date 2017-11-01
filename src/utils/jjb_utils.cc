//The indices of data and speq have been decremented by 1 so they start
//from 0 not 1 JJB 10/1/2017
#include <math.h>
#include <stdio.h>
#include "jjb_utils.h"

double *dmatrix(int n){
	double *m1 = new double[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
double **dmatrix(int n,int m){
	double **m1 = new double*[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m);
	return m1;
}
double ***dmatrix(int n,int m,int l){
	double ***m1 = new double**[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l);
	return m1;
}
double ****dmatrix(int n,int m,int l,int k){
	double ****m1 = new double***[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l,k);
	return m1;
}
double ***d3array(int n,int m,int l){//as dmatrix(n,m,l) but guaranteed contiguity of data 
	double *m0 = new double[n*m*l];//allocate all required storage
	double ***m1 = new double**[n];
	for(int i=0;i<n;i++){
		m1[i] = new double*[m];
		for(int j=0; j<m;j++){
			m1[i][j]=m0+l*(i*m+j);
		}
	}
	return m1;
}
void free_d3array(double ***m1){
	delete[] m1[0][0];
	delete[] m1[0];
	delete[] m1;
}
void delmatrix(double **m1,int n){
	for(int i=0;i<n;i++) delete[] m1[i];
	delete [] m1;
}
void delmatrix(double ***m1,int n,int m){
	for(int i=0;i<n;i++) delmatrix(m1[i],m);
	delete [] m1;
}
void delmatrix(double ****m1,int n,int m,int l){
	for(int i=0;i<n;i++) delmatrix(m1[i],m,l);
	delete [] m1;
}
void tidy_angles(Angles &theta){
	for(int i=0;i<3;i++){
		if(fabs(theta[i])>100) theta[i] -= TPi*int(theta[i]/TPi);
		while (theta[i]< 0. ) theta[i] += TPi;
		while (theta[i]> TPi) theta[i] -= TPi;
	}
}
int quadratic(double a,double b,double c,double *roots){
	//smallest root 1st, diff of roots in roots[2]
	double discr=b*b-4*a*c;
	if(discr<0){
		printf("discr<0 in quadratic: %g %g %g %g\n",a,b,c,discr);
		return 1;
	}
	discr=sqrt(MAX(0,discr));
	double q=-.5*(b+SGN(b)*discr),r1=q/a,r2=c/q;
	if(fabs(a)<1e-30) r1=1e30;//effectively a linear equation
	if(fabs(q)<1e-30) r2=-r1;//effectively a perfect square=0
	roots[0]=MIN(r1,r2); roots[1]=MAX(r1,r2); roots[2]=roots[1]-roots[0];
	return 0;
}
double lntp(double *x,double *y,int nin,double xp){
	int bot=0,top=nin-1;
	if((xp-x[0])*(x[nin-1]-xp)<0){
		printf("lntp: %f %f %f \n",x[0],x[nin-1],xp);
		//exit(0);
	}
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((x[top]-xp)*(xp-x[n])>=0) bot=n;
		else top=n;
	}
	double f=(x[top]-xp)/(x[top]-x[bot]);
	return f*y[bot]+(1-f)*y[top];
}
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void fourn(double *data, unsigned long *nn, int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
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
						tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
						tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
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
void rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,
	   unsigned long nn3, int isign)
{
	//void fourn(float data[], unsigned long nn[], int ndim, int isign);
	//void nrerror(char error_text[]);
	unsigned long i1,i2,i3,j1,j2,j3,nn[3],ii3;
	double theta,wi,wpi,wpr,wr,wtemp;
	double c1,c2,h1r,h1i,h2r,h2i;

	if (1+&data[nn1-1][nn2-1][nn3-1]-&data[0][0][0] != nn1*nn2*nn3)
		printf("rlft3: problem with dimensions or contiguity of data array\n");
	c1=0.5;
	c2 = -0.5*isign;
	theta=isign*(6.28318530717959/nn3);
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	nn[0]=nn1;
	nn[1]=nn2;
	nn[2]=nn3 >> 1;
	if (isign == 1) {
		fourn(&data[0][0][0]-1,nn-1,3,isign);
		for (i1=1;i1<=nn1;i1++)
			for (i2=1,j2=0;i2<=nn2;i2++) {
				speq[i1-1][++j2-1]=data[i1-1][i2-1][0];
				speq[i1-1][++j2-1]=data[i1-1][i2-1][1];
			}
	}
	for (i1=1;i1<=nn1;i1++) {
		j1=(i1 != 1 ? nn1-i1+2 : 1);
		wr=1.0;
		wi=0.0;
		for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
			for (i2=1;i2<=nn2;i2++) {
				if (i3 == 1) {
					j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
					h1r=c1*(data[i1-1][i2-1][0]+speq[j1-1][j2-1]);
					h1i=c1*(data[i1-1][i2-1][1]-speq[j1-1][j2]);
					h2i=c2*(data[i1-1][i2-1][0]-speq[j1-1][j2-1]);
					h2r= -c2*(data[i1-1][i2-1][1]+speq[j1-1][j2]);
					data[i1-1][i2-1][0]=h1r+h2r;
					data[i1-1][i2-1][1]=h1i+h2i;
					speq[j1-1][j2-1]=h1r-h2r;
					speq[j1-1][j2]=h2i-h1i;
				} else {
					j2=(i2 != 1 ? nn2-i2+2 : 1);
					j3=nn3+3-(i3<<1);
					h1r=c1*(data[i1-1][i2-1][ii3-1]+data[j1-1][j2-1][j3-1]);
					h1i=c1*(data[i1-1][i2-1][ii3]-data[j1-1][j2-1][j3]);
					h2i=c2*(data[i1-1][i2-1][ii3-1]-data[j1-1][j2-1][j3-1]);
					h2r= -c2*(data[i1-1][i2-1][ii3]+data[j1-1][j2-1][j3]);
					data[i1-1][i2-1][ii3-1]=h1r+wr*h2r-wi*h2i;
					data[i1-1][i2-1][ii3]=h1i+wr*h2i+wi*h2r;
					data[j1-1][j2-1][j3-1]=h1r-wr*h2r+wi*h2i;
					data[j1-1][j2-1][j3]= -h1i+wr*h2i+wi*h2r;
				}
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
	}
	if (isign == -1)
		fourn(&data[0][0][0]-1,nn-1,3,isign);
}
