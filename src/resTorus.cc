/* Functionality for 1:1 resonantly trapped tori pewr Binney (2016)
 * Coded by James Binney March 2016 */
#include "resTorus.h"

static double PI=acos(-1),PIH=PI/2,TPI=2*PI;

static double piIt(double diff){//deal with TPI ambiguities
	if(fabs(diff-TPI)<fabs(diff)) return diff-TPI;
	else if(fabs(diff+TPI)<fabs(diff)) return diff+TPI;
	else return diff;
}
static double swap_theta12(double t,Angles &theta1){//get onto best theta on unperturbed torus 
	Angles theta2; theta2[0]=TPI-theta1[0]; theta2[1]=PI-theta1[1]; theta2[2]=theta1[2];
	double diff1=theta1[0]-theta1[1]-t, diff2=theta2[0]-theta2[1]-t;
	diff1=piIt(diff1); diff2=piIt(diff2);
	if(fabs(diff2)<fabs(diff1)){
		theta1=theta2; return diff2;
	} else return diff1;
}
resTorus::resTorus(Torus *Tg,double *Jg,int no,Actions Jbin,double om,double *hr){
// Framework common to all 1:1 resonant torus of given E. Still have to specify I
	Tgrid=Tg; // pointer to 1d array of tori: usually one each side of trapping zone
	Jrgrid=Jg;// pointer to radial actions of above tori
	norb=no;  // number of tori in array
	Jb=Jbin;  // actions of the perfectly resonant torus
	omp=om;   // dOmega_N/dJ
	hres=hr;  // h_(2,-2) and h_(4,-4) and their 1st 2 derivatives
	if(hres[0]>0){
		O=true; t1c=0;
	} else{
		O=false; t1c=PIH;
	}
	lJ=-1;  // We haven't yet evaluated resonant action
	lOmega=-1;// We haven't yet evaluated resonant frequency
	getImin_max();//set limits to I
	nangle=200;//nangle must be even
	t1grd=new double[nangle]; trgrd=new double[nangle]; Dgrd=new double[nangle];
}
void resTorus::reset(Actions Jbin,double om){
	Jb=Jbin; omp=om;
	lJ=-1;  // We haven't yet evaluated resonant action
	lOmega=-1;// We haven't yet evaluated resonant frequency
	getImin_max();//set limits to I
}
int resTorus::setI(double Iin){
	/* This must be called after resTorus::resTorus() to specify
	 * amplitude of libration */
	I=Iin; lJ=-1;
	if(I<Imin || I>Imax){
		printf("I out of range: I, Imin, max %f %f %f\n",I,Imin,Imax);
		return 0;
	}
	set_Dtheta1();
	store_theta(); return 1;
}
void resTorus::getImin_max(){
	double I0=-2*(hres[0]-hres[3])-pow(hres[1]-hres[4],2)/(.5*omp-hres[2]+hres[5]);
	double I1= 2*(hres[0]+hres[3])-pow(hres[1]+hres[4],2)/(.5*omp+hres[2]+hres[5]);
	Imax=.9999*MAX(I0,I1); Imin=MIN(I0,I1);
	if(Imin<0) Imin*=.9999; else Imin*=1.0001;
//	printf("Imin Imax: %f %f\n",Imin,Imax);
}
void resTorus::getDelta(double t1p,double *Delta){//returns J1p-Jb[0]
	double cos2t1p=cos(2*t1p),cos4t1p=cos(4*t1p);
	double a=.5*omp+hres[2]*cos2t1p+hres[5]*cos4t1p;
	double b=2*(hres[1]*cos2t1p+hres[4]*cos4t1p);
	double c=2*(hres[0]*cos2t1p+hres[3]*cos4t1p)-I;
	quadratic(a,b,c,Delta);
}
#define TINY 1e-6
void resTorus::getJs(Actions &J1,Actions &J2){//returns extreme unperturbed actions
	double Delta[3]; getDelta(t1c,Delta);
	J1[2]=J2[2]=Jb[2];
	J1[0]=MAX(TINY,Jb[0]+Delta[0]); J1[1]=MAX(TINY,Jb[1]-Delta[0]);
	J2[0]=MAX(TINY,Jb[0]+Delta[1]); J2[1]=MAX(TINY,Jb[1]-Delta[1]);
}
void resTorus::Dtheta1_fn(double t1p,double &f,double &fp) const{
	double cos2t1p=cos(2*t1p),cos4t1p=cos(4*t1p);
	double a=.5*omp+hres[2]*cos2t1p+hres[5]*cos4t1p;
	double b=2*(hres[1]*cos2t1p+hres[4]*cos4t1p);
	double c=2*(hres[0]*cos2t1p+hres[3]*cos4t1p)-I;
	f=b*b-4*a*c;
	double sin2t1p=sin(2*t1p),sin4t1p=sin(4*t1p);
	double ap=-hres[2]*2*sin2t1p-hres[5]*4*sin4t1p;
	double bp=-2*(hres[1]*2*sin2t1p+hres[4]*4*sin4t1p);
	double cp=-2*(hres[0]*2*sin2t1p+hres[3]*4*sin4t1p);
	fp=2*b*bp-4*(ap*c+a*cp);
	//printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",1e5*hres[0],1e5*hres[1],1e5*hres[2],1e5*hres[3],1e5*hres[4],1e5*hres[5],t1p/PI,1e5*f,1e5*fp);
}
void resTorus::set_Dtheta1(){//find amplitude of excursions in theta1
	//first find limit to excursions
	double tol=1.e-10,f,f1,fp,t1;
	Dtheta1_fn(t1c,f,fp);
	if(f<=0){
		Dt1=0; return;
	}
	if(O){
		t1=.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0){
			t1+=.05*(PIH-t1); Dtheta1_fn(t1,f1,fp);
	//		printf("O: %g %g\n",t1,f1);
		}
	}else{
		t1=t1c-.05; Dtheta1_fn(t1,f1,fp);
		while(f1>0){
			t1-=.1*t1; Dtheta1_fn(t1,f1,fp);
	//		printf("!O: %g %g\n",t1,f1);
		}
	}			
	if(f1>=0){
		Dt1=PIH; return;
	}
	if(fabs(f1)>tol) t1=rtsafe(this,&resTorus::Dtheta1_fn,t1c,t1,tol);
	if(O) Dt1=.999999*t1; else Dt1=.999999*fabs(PIH-t1);
}
double resTorus::librationAction(){//returns libration J
	if(lJ!=-1) return lJ;
	int np=80;
	double dt=Dt1/(float)(np);
	double Delta[3],t1p=t1c; getDelta(t1p,Delta);
	double J1p=0,A=Delta[2];
	for(int i=0;i<np;i++){
		t1p+=dt;
		getDelta(t1p,Delta);
		double B=Delta[2];
		J1p+=.5*(A+B)*dt;
		A=B;
	}
	J1p/=PI; return J1p;
}
double resTorus::theta_dot(double t1p,double D){//with D given
	return omp*D+2*(hres[1]+hres[2]*D)*cos(2*t1p)+2*(hres[4]+hres[5]*D)*cos(4*t1p);
}
double resTorus::rk_step(double &t1p,double h,double &dtheta,int i){
	double th,k[4],Delta[3];
	getDelta(t1p,Delta);
	k[0]=h*theta_dot(t1p,Delta[i]);
	th=t1p+.5*k[0];
	getDelta(th,Delta);
	k[1]=h*theta_dot(th,Delta[i]);
	th=t1p+.5*k[1];
	getDelta(th,Delta);
	k[2]=h*theta_dot(th,Delta[i]);
	th=t1p+k[2];
	getDelta(th,Delta);
	k[3]=h*theta_dot(th,Delta[i]);
	dtheta=(k[0]+2*(k[1]+k[2])+k[3])/6.;
	if(dtheta>0) dtheta*=-1;
	t1p+=dtheta;
	getDelta(t1p,Delta);
	return Delta[i];
}	
void resTorus::store_theta(){//Store what's needed by FullMap
	double D,Delta[3],t1p; 
	double s=0,x=0,dx=PI/(double)(nangle-1);
	t1p=t1c; getDelta(t1p,Delta); D=Delta[0];
	double tdot=theta_dot(t1p,D);
	double A=cos(x)/tdot,B;
	t1grd[0]=t1p; trgrd[0]=s; Dgrd[0]=D;
	for(int i=1;i<nangle;i++){
		x+=dx;
		t1p=t1c+Dt1*sin(x); getDelta(t1p,Delta);
		if(i<nangle/2) D=Delta[0]; else D=Delta[1];
		tdot=theta_dot(t1p,D);
		if(tdot!=0){
			double B1=cos(x)/tdot;
			if(fabs(B-A)<2*A) B=B1; else B=A;
		}
		s+=.5*(A+B); A=B;
		t1grd[i]=t1p; trgrd[i]=s; Dgrd[i]=D;
	}
	lOmega=PI/trgrd[nangle-1];
	for(int i=0;i<nangle;i++) trgrd[i]*=lOmega;
	lOmega/=(Dt1*dx);
}
double resTorus::librationAngle(double t1p,double Delta){//returns libration theta
	double psi,td=theta_dot(t1p,Delta);
	int nh=nangle/2;
	if(td>0)
		psi=lntp(t1grd,trgrd,nh,t1c+fabs(t1p-t1c));
	else
		psi=lntp(t1grd+nh,trgrd+nh,nh,t1c+fabs(t1p-t1c));
	if(t1p<t1c) psi=2*PI-psi;
	return psi;
}
double resTorus::from_librationAngle(double tr,double &D){
	/* Return t1=t_r-t_z and action offset Jr-Jb[0] given the true
	 * libration angle tr */
	while(tr>2*PI) tr-=2*PI;
	while(tr<0) tr+=2*PI;
	bool sn=false;
	if(tr>PI){
		tr=2*PI-tr; sn=true;
	}
	int bot=0,top=nangle-1;
	while(fabs(top-bot)>1){
		int n=(top+bot)/2;
		if((trgrd[top]-tr)*(tr-trgrd[n])>=0) bot=n;
		else top=n;
	}
	double f;
	if(trgrd[bot]>=trgrd[top]) f=0;
	else f=(trgrd[top]-tr)/(trgrd[top]-trgrd[bot]);
	double t1p=f*t1grd[bot]+(1-f)*t1grd[top];
	if(sn) t1p=2*t1c-t1p;
	D=f*Dgrd[bot]+(1-f)*Dgrd[top];
	return t1p;
}
GCY resTorus::FullMap(const Angles& A){
	/* First angle is libration angle, then theta_z, theta_phi.
	 * The main job is to get theta_r from lib angle.  */
	double D,t1p=from_librationAngle(A[0],D),Jr=Jb[0]+D;
	Angles theta; theta[0]=t1p+A[1]; theta[1]=A[1]; theta[2]=A[2];
	Torus T=InterpTorus(Tgrid,Jrgrid,norb,Jr);
	GCY gcy=T.FullMap(theta);
	return gcy;
}
Actions resTorus::actions(const Angles& A){
	double D,t1p=from_librationAngle(A[0],D),Jr=Jb[0]+D;
//	Angles theta; theta[0]=t1p+A[1]; theta[1]=A[1]; theta[2]=A[2];
	Torus T=InterpTorus(Tgrid,Jrgrid,norb,Jr);
	Actions J=T.actions();
	return J;
}	
double resTorus::librationOmega(void){
	return lOmega;
}
void resTorus::SOS(ostream& sfile,int np){//work round with true angle
	double Delta[3],tr;
	for(int i=0;i<np;i++){
		if(O) tr=i/(double)(np-1)*(PI-1.e-6);
		else tr=i/(double)(np-1)*(TPI-1.e-6);
		double D,t1p=from_librationAngle(tr,D),J1p=Jb[0]+D;
		Torus T=InterpTorus(Tgrid,Jrgrid,norb,J1p);
		T.resSOS(sfile,t1p);
	}
}
void resTorus::SOSr(ostream& sfile,int np){//work round with true angle
	double Delta[3],tr;
	for(int i=0;i<np;i++){
		if(O) tr=i/(double)(np-1)*(PI-1.e-6);
		else tr=i/(double)(np-1)*(TPI-1.e-6);
		double D,t1p=PI+from_librationAngle(tr,D),J1p=Jb[0]+D;
		Torus T=InterpTorus(Tgrid,Jrgrid,norb,J1p);
		T.resSOS(sfile,t1p);
	}
}
bool resTorus::search_eval(Position &Rzphi,double tp,int iD,double &diff1, double &diff2,
			   Angles &theta1,Angles &theta2){
	/* Computes differences diff1 and diff2 between tp and t_r-t_z when
	 * unperturbed torus specified by tp passes through Rzphi. There are 4
	 * possible values of tr-tz but TORUS::containsPoint returns only 2.
	 * swap_theta12 picks out the required pair. */
	double Delta[3]; getDelta(tp,Delta);
	Torus T=InterpTorus(Tgrid,Jrgrid,norb,Jb[0]+Delta[iD]);
	bool ok=T.containsPoint_Ang(Rzphi,theta1,theta2);
	if(!ok){
		theta1[0]=0; theta1[1]=0; theta2[0]=0; theta2[1]=0;
	}else{
		diff1=swap_theta12(tp,theta1); diff2=swap_theta12(tp,theta2);
	}
	return ok;
}
#define SMALL 1e-7
Angles resTorus::lockit(Position &Rzphi,double tpa,double tpb,double &diffa,double &diffb,
			int iD,int k,double *tm){
	double t,diff1=diffa,diff2=diffb; Angles theta1,theta2;
	while(fabs(tpb-tpa)>SMALL){
		t=.5*(tpb+tpa);
		search_eval(Rzphi,t,iD,diff1,diff2,theta1,theta2);
		if(k==1){//our quarry lies near theta1
			if(diff1*diffb<=0){
				tpa=t; diffa=diff1;
			}else{
				tpb=t; diffb=diff1;
			}
		}else{//our quarry lies near theta2
			if(diff2*diffb<=0){
				tpa=t; diffa=diff2;
			}else{
				tpb=t; diffb=diff2;
			}
		}
	}
	(*tm)=.5*(tpa+tpb);
	double Delta[3]; getDelta((*tm),Delta);
	if(k==1){
		theta1[0]=librationAngle((*tm),Delta[iD]); return theta1;
	}else{
		theta2[0]=librationAngle((*tm),Delta[iD]); return theta2;
	}
}
#define N 50  //
int resTorus::containsPoint_Ang(Position &Rzphi,Angles *A,double *tm){
	/* Returns number of angles at which trapped torus visits Rzphi, with
	 * angles in A and in tm corresponding values of theta1=theta_r-theta_z.
	 * Returned value can be 0, 2 (normal) or 4*/
	Angles theta1,theta2;
	double tpa,tpb,diff1b,diff2b,diff1a,diff2a;
	int n=0; bool oka=false,okb;
	for(int iD=0;iD<2;iD++){
		for(int i=0;i<N;i++){
			tpb=t1c+Dt1*(-1+2*i/(float)(N-1));
			okb=search_eval(Rzphi,tpb,iD,diff1b,diff2b,theta1,theta2);
			if(i>0 && oka && okb){//found for both values tp
				if(diff1a*diff1b<0){//a root with theta1
					//printf("locking 1: %f %f %f %f\n",tpa,tpb,diff1a,diff1b);
					A[n]=lockit(Rzphi,tpa,tpb,diff1a,diff1b,iD,1,tm+n);
					if(fabs(diff1a-diff1b)<.01) n++;
				}
				if(diff2a*diff2b<0){//a root with theta2
					//printf("locking 2: %f %f %f %f\n",tpa,tpb,diff2a,diff2b);
					A[n]=lockit(Rzphi,tpa,tpb,diff2a,diff2b,iD,2,tm+n);
					if(fabs(diff2a-diff2b)<.01) n++;
				}
				if(n==5){
					printf("too many crossings in resTorus::containsPoint_Ang\n");
					return n;
				}
			}
			tpa=tpb; oka=okb; diff1a=diff1b; diff2a=diff2b;
		}
	}
	return n;
}
int resTorus::containsPoint(Position &Rzphi,Velocity *V,Angles *A,double *tm){
	int n=containsPoint_Ang(Rzphi,A,tm);
	if(n!=0){
		for(int i=0;i<n;i++){
			GCY gcy=FullMap(A[i]);
			V[i][0]=gcy[3]; V[i][1]=gcy[4]; V[i][2]=gcy[5];
		}
	}
	return n;
}
