/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
#include "eTorus.h"

eTorus::eTorus(Torus &Tin,Potential *Phi,bar_pot *bar,const double Ompin){
	tmax=10; T=Tin; Omp=Ompin;
	int nf=pow(2,6);
	getH0Hn(Phi,bar,nf,0);
}
eTorus::eTorus(Actions J,Potential *Phi,bar_pot *bar,const double Ompin,const double tolJ){
	tmax=10;
	AutoFit(J,Phi,bar,Ompin,tolJ);
}
int eTorus::AutoFit(Actions J,Potential *Phi,bar_pot *bar,const double Ompin,const double tolJ){
	tmax=10; Omp=Ompin;
	int ok=T.AutoFit(J,Phi,tolJ);
	int nf=pow(2,6);
	getH0Hn(Phi,bar,nf,0);
	return ok;
}
void eTorus::reset(Torus &Tin,Potential *Phi,bar_pot *bar){
	tmax=10; T=Tin;
	int nf=pow(2,6);
	getH0Hn(Phi,bar,nf,0);
}
double eTorus::H0(Potential *Phi,bar_pot *bar,Angles &thetas){//returns full Hamiltonian
	GCY w=T.FullMap(thetas);
	return .5*(pow(w[3],2)+pow(w[4],2)+pow(w[5],2))+(*Phi)(w[0],w[1])-Omp*T.action(2)
			+bar->Phi(w[0],w[1])*cos(2*w[2]);
}
void eTorus::getH0Hn(Potential *Phi,bar_pot *bar,int nf,int ifp){//Fourier analyses H0
	int nfr=nf,nfz=nf,nfp=nf/4;
	double ***h; h=d3array(nfr,nfz,nfp);
	Angles thetas;
	double dtr=2*Pi/(double)nfr, dtz=2*Pi/(double)nfz, dtp=2*Pi/(double)nfp;
	for(int i=0;i<nfr;i++){
		thetas[0]=i*dtr;
		for(int j=0;j<nfz/2;j++){
			thetas[1]=j*dtz;
			for(int k=0;k<nfp/2;k++){
				thetas[2]=k*dtp;
				double a=H0(Phi,bar,thetas);
				h[i][j][k]=a; h[i][nfz/2+j][k]=a;//N-S symmetry
				h[i][j][nfp/2+k]=a; h[i][nfz/2+j][nfp/2+k]=a;//bi-symmetry
			}
		}
	}
	double **speq; speq=dmatrix(nfr,2*nfz);
	rlft3(h,speq,nfr,nfz,nfp,1);
	double top[tmax];
	int ntop=0;
	double hmax=0;
	for(int i=0;i<nfr;i++){//find largest perturbing terms
		for(int j=0;j<nfz;j++){
			for(int k=0;k<nfp/2;k++){
				double s=sqrt(pow(h[i][j][2*k],2)+pow(h[i][j][2*k+1],2));
				if(ntop<tmax || s>hmax){
					hmax=insert(top,ntop,s,i,j,k);
				}
			}
		}
	}
	double n=(nfr*nfz*nfp);
	for(int i=0;i<ntop;i++){
		hN[i]=h[iN[i]][jN[i]][2*kN[i]]/n;
		if(ifp) printf("(%d %d %d) (%g %g)\n",fq(nfr,iN[i]),fq(nfz,jN[i]),fq(nfp,kN[i]),
			       100*hN[i],h[iN[i]][jN[i]][2*kN[i]+1]/n);
	}
	for(int i=0;i<ntop;i++){
		iN[i]=fq(nfr,iN[i]); jN[i]=fq(nfz,jN[i]); kN[i]=fq(nfp,kN[i]);
	}
	delmatrix(speq,nfr); free_d3array(h);
}
int eTorus::fq(int nf,int i){
	if(i<nf/2) return i;
	else return i-nf;
}
double eTorus::insert(double *top,int &ntop,double s,int i,int j,int k){
	if(ntop==0){
		top[0]=s; iN[0]=i; jN[0]=j; kN[0]=k;
		ntop++; return top[0];
	}
	int l=0;
	while(l<ntop && s<=top[l]) l++;
	if(l==ntop){//line isn't stronger than any previous line
		if(ntop==tmax){
			return top[ntop-1];//no room for this line
		} else {//add line to end of list
			top[l]=s; iN[l]=i; jN[l]=j; kN[l]=k;
			ntop++;
		}
	} else {//we should insert current line
		if(ntop<tmax) ntop++;//we are adding rather than replacing a line
		for(int m=ntop-1; m>l; m--){
			top[m]=top[m-1]; iN[m]=iN[m-1]; jN[m]=jN[m-1]; kN[m]=kN[m-1];
		}
		top[l]=s; iN[l]=i; jN[l]=j; kN[l]=k;
	}
	return top[ntop-1];
}
Frequencies eTorus::omega() const{
	Frequencies Om=T.omega(); Om[2]-=Omp;
	return Om;
}
void eTorus::LevCof3D(const Angles &theta,Matrix<double,3,3> &dJdt,const Position &Rzphi,
		      PSPT &QP,double &chi,
		      Vector<double,3> &B,Matrix<double,3,3> &A,Matrix<double,3,3> &dQdt) const{
	T.LevCof3D(theta,dJdt,Rzphi,QP,chi,B,A,dQdt);
}
void eTorus::LevCof3DTrue(const Angles &theta,Matrix<double,3,3> &dJdt,const Position &Rzphi,
		      PSPT &QP,double &chi,
		      Vector<double,3> &B,Matrix<double,3,3> &dQdt) const{
	T.LevCof3DTrue(theta,dJdt,Rzphi,QP,chi,B,dQdt);
}
double eTorus::DistancetoPoint(const Position &Rzphi,Angles &theta) const{
	return T.DistancetoPoint_Ang(Rzphi,theta);
}
bool eTorus::containsPoint(const Position& Rzphi,Velocity &v1,Velocity &v2,double &d1,
		   Angles &theta1,Angles &theta2,Velocity &v3,Velocity &v4,
		   double &d2,Angles &theta3,Angles &theta4) const{
	return T.containsPoint(Rzphi,v1,v2,d1,theta1,theta2,v3,v4,d2,theta3,theta4);
}
eTorus& eTorus::operator=(const eTorus& T2){
	T=T2.T; hN=T2.hN; iN=T2.iN; jN=T2.jN; kN=T2.kN; Omp=T2.Omp;
	return *this;
}
eTorus& eTorus::operator*=(const double& x){
	T *=x; hN *=x;
	return *this;
}
const eTorus eTorus::operator*(const double& x){
	eTorus T2; T2=*this;
	T2*=x;
	return T2;
}
eTorus& eTorus::operator+=(const eTorus& T2){
	int start=0;
	for(int i=0;i<tmax;i++){
		for(int j=start;j<tmax;j++){
			if(iN[i]==T2.iN[j] && jN[i]==T2.jN[j] && kN[i]==T2.kN[j]){
				hN[i]+=T2.hN[j]; if(j==start) start++; break;
			}
		}
	}
	T+=T2.T;
	return *this;
}
const eTorus eTorus::operator+(const eTorus& T2){
	eTorus T1; T1=*this;
	T1+=T2;
	return T1;
}
eTorus& eTorus::operator-=(const eTorus& T2){
	int start=0;
	for(int i=0;i<tmax;i++){
		for(int j=start;j<tmax;j++){
			if(iN[i]==T2.iN[j] && jN[i]==T2.jN[j] && kN[i]==T2.kN[j]){
				hN[i]-=T2.hN[j]; if(j==start) start++; break;
			}
		}
	}
	T-=T2.T;
	return *this;
}
const eTorus eTorus::operator-(const eTorus& T2){
	eTorus T1; T1=*this;
	T1-=T2;
	return T1;
}
eTorus Interp_eTorus_mn(eTorus **Tgrid,int m,int n,Actions Jbar,Actions dJ,Actions J){
	//m*n rectangles uniformly spaced in Jp and sqrt(Jr).
	//Jbar centre of bottom left rectangle, dJ vector across its diagonal
	double j0=sqrt(MAX(J[0],.00001));
	int np=(J[2]-(Jbar[2]-.5*dJ[2]))/dJ[2];//# of our rectangle
	int nr=(j0-(Jbar[0]-.5*dJ[0]))/dJ[0];//# of our rectangle
	if(np<0 || np>m-2 || nr<0 || nr>n-2){
		printf("\noff grid in Interp_eTorus_mn: cell # n=(%d,%d) (%g %g)\n",np,nr,J[2],J[0]);
		np=MAX(0,np); np=MIN(m-2,np); nr=MAX(0,nr); nr=MIN(n-2,nr);
		exit(0);
	}
	eTorus T; T=Tgrid[np][nr];
	double fp=(Jbar[2]+(np+.5)*dJ[2]-J[2])/dJ[2];//distance from right = weight of left
	double fr=(Jbar[0]+(nr+.5)*dJ[0]-j0)/dJ[0];//distance from top = weight of bottom
	fr=MIN(1,fr);
	for(int i=0;i<2;i++) {
		if(i==1){
			fp=1-fp; fr=1-fr;
		}
		for(int k=0;k<2;k++) {
			if(k==1) fr=1-fr;
			if(fp>1.000001 || fr>1.000001){
				printf("Interp_eTorus_mn error: %f %f\n",fp,fr);
				printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",
				       Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
				exit(0);}
			if(i+k==0) T*=fp*fr;
			else T+=Tgrid[np+i][nr+k]*(fp*fr);
		}
	}
	return T;
}
bool eTorus::write_ebf(const string filename,const string torusname){
	bool ok=T.write_ebf(filename,torusname);
	if(ok){
		ebf::Write(filename,"/"+torusname+"/Omp",&Omp,"a","",1);
		ebf::Write(filename,"/"+torusname+"/hN",&hN[0],"a","",tmax);
		ebf::Write(filename,"/"+torusname+"/iN",&iN[0],"a","",tmax);
		ebf::Write(filename,"/"+torusname+"/jN",&jN[0],"a","",tmax);
		ebf::Write(filename,"/"+torusname+"/kN",&kN[0],"a","",tmax);
	}
	return ok;
}
bool eTorus::read_ebf(const string filename, const string torusname){
	bool ok=T.read_ebf(filename,torusname);
	if(ok){
		ebf::Read(filename,"/"+torusname+"/Omp",&Omp,1);
		ebf::Read(filename,"/"+torusname+"/hN",&hN[0],tmax);
		ebf::Read(filename,"/"+torusname+"/iN",&iN[0],tmax);
		ebf::Read(filename,"/"+torusname+"/jN",&jN[0],tmax);
		ebf::Read(filename,"/"+torusname+"/kN",&kN[0],tmax);
	}
	return ok;
}
