/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
#include "iTorus.h"
#include <string.h>
#define SMALL 1e-6
/*
void plot_orb(FILE *ofile,iTorus &iT,double R0,double phi0){
	int NS=5000; float x[NS],y[NS],Xmax=10,PI=acos(-1);
	Frequencies Om=iT.omega(); Actions J=iT.actions();
	printf("Acts: %f %f\n",J[0],J[2]);
	double dt=100*PI/fabs(Om(0))/(double)NS;
	fprintf(ofile,"%f %f %f %f %f\n",(NS-1)*dt,J[0],J[2],R0,phi0);
	for(int i=0;i<NS;i++){
		double t=i*dt;
		Angles A;
		A[0]=Om[0]*t; A[1]=Om[1]*t; A[2]=Om[2]*t;
		GCY gcy=iT.FullMap(A);
		x[i]=gcy[0]*cos(gcy[2]); y[i]=gcy[0]*sin(gcy[2]);
		for(int j=0;j<6;j++){
			fprintf(ofile,"%f ",gcy[j]);
		}
		fprintf(ofile,"\n");
	}
	//plots(NS,x,y,-Xmax,Xmax,-Xmax,Xmax,"x/kpc",5,"y/kpc",5,-.9,10);
	plots(NS,x,y,-9,-6,2,5,"x/kpc",5,"y/kpc",5,-.9,10);
}
*/
double fix_Jbar(double *Jrs,int nJr,eTorus **Tgrid,int np,int nr,Actions Jgrid,
		Actions dJ,Position Rzphi,double Jbot,double Jtop){
	//Finds Jp for the closed orbit through Rzphi. Jtop/Jbot are
	//largest/smallest Jp values to consider
	double Jbar=.5*(Jtop+Jbot),Japo,Jperi;
	Actions J;
	int dj=nJr/4,j=nJr-1-2*dj;
	for(int i=0;i<3;i++){//work down from the most eccentric torus finding Jperi & Japo
		J[0]=Jrs[j]; J[1]=.0025; J[2]=Jbar;
		iTorus T(J,Tgrid,np,nr,Jgrid,dJ);
		T.get_crit_Jp(Rzphi,Japo,Jperi,Jbot,Jtop);
		Jbar=.5*(Japo+Jperi); j+=dj;
	}
	printf("Jr Jbar %f %f\n",J[0],Jbar);
	return Jbar;
}
void tabulate_Jr_Jp(GalaxyPotential *Phi,const double R,double *Jrs,double *Jps,int n){
	//For input Jphi values returns Jr values that put peri/apo at R
	double eps=1e-10;
	Torus T;
	for(int i=0;i<n;i++){
		double Jp=Jps[i],E=.5*pow(Jp/R,2)+(*Phi)(R,0);
		double Rc=Phi->RfromLc(Jp);
		double dPdR,dPdz,Rp=2*Rc-R;
		double f=E-.5*pow(Jp/Rp,2)-(*Phi)(Rp,0,dPdR,dPdz);
		int k=0,nint=40;
		while(fabs(f)>eps && k<10){
			double df=pow(Jp/Rp,2)/Rp-dPdR;
			double dR=-f/df;
			while(fabs(dR/Rp)>.4) dR*=.9;
			Rp+=dR; f=E-.5*pow(Jp/Rp,2)-(*Phi)(Rp,0,dPdR,dPdz); k++;
		}
		double Jr=0,A=0,B,theta=0,dtheta=Pi/(double)nint;
		double DR=.5*fabs(R-Rp),Rbar=.5*(Rp+R);
		for(int j=0;j<nint;j++){//trapezium rule
			theta+=dtheta; double r=Rbar-DR*cos(theta);
			if(j<nint-1) B=sin(theta)*sqrt(2*(E-.5*pow(Jp/r,2)-(*Phi)(r,0)));
			else B=0;
			Jr+=(A+B); A=B;
		}
		Jrs[i]=Jr*.5*DR*dtheta/Pi;
	}
}
int main(int narg,char **argv){
	if(narg==1){
		printf("You must specify Omega_p\nAlso: ensure that a subdirectory with that value as a name exists\n"); return 0;
	}
	int omp; sscanf(argv[1],"%d",&omp);
	char fname[40],dirname[5]; sprintf(dirname,"%d/",omp);
	double Omp=0.001*omp;
	int resN[3]={0,0,2}; double PI=acos(-1),TPI=2*PI,torad=PI/180;
	ifstream Pfile; Pfile.open("../pot/PJM11_best.Tpot");
	if(!Pfile.is_open()) printf("can't open Tpot file\n");
	GalaxyPotential *Phi = new GalaxyPotential(Pfile);
	Pfile.close();
	//printf("unit %f\n",1/Units::kms); if(PI>0) return 0;
	DiskPar p0=Phi->DiskParameter(0),p1=Phi->DiskParameter(1);
	double Rd=(p0[0]*p0[1]+p1[0]*p1[1])/(p0[0]+p1[0]);//weighted mean of disc scales
	double Rb=0.7*Rd;
	double Ab=3*exp(2)/20.*pow(220./239.,2)*pow(1.5/Rb,3)*.4;
	bar_pot bar(.9,Rb,239./978.,Ab);
	double tolJ=0.003,Jrmin=.0002,Jrmax=.19;
	double Jpmin=1.,Jpmax=3;
	Actions J,Jgrid,dJ; J[1]=.0025;
	double R0=8; J[0]=.07; J[2]=.5;
	Position Rzphi; Rzphi[0]=R0; Rzphi[1]=.01; Rzphi[2]=155*torad;
	double dPhiR,dPhiz; (*Phi)(R0,0,dPhiR,dPhiz);
	double Vc=sqrt(R0*dPhiR),JpSun=R0*Vc,Usun=11,Vsun=12+Vc/Units::kms; printf("U,Vsun: %f %f\n",Usun,Vsun);
	int nr=30,np=20;
	dJ[0]=(sqrt(Jrmax)-sqrt(Jrmin))/(double)(nr-1); dJ[1]=0; dJ[2]=(Jpmax-Jpmin)/(double)(np-1);
	Jgrid[0]=sqrt(Jrmin)+.5*dJ[0];  Jgrid[1]=J[1]; Jgrid[2]=Jpmin+.5*dJ[2];
	eTorus **Tgrid=PJM::matrix<eTorus>(np,nr);
	strcpy(fname,dirname); strcat(fname,"nores.ebf");
	bool writeit=true;
	if(writeit){
		ebf::Write(fname,"/dJ",&dJ[0],"w","",1);
		printf("Computing eTorus grid:\n");
		for(int i=0;i<np;i++){
			printf("%d ",i);
			J[2]=Jpmin+i*dJ[2];
			for(int j=0;j<nr;j++){
				J[0]=pow(sqrt(Jrmin)+j*dJ[0],2);
				if(i==0) printf("%f ",J[0]);
				Tgrid[i][j].AutoFit(J,Phi,&bar,Omp,tolJ);
				char lab[7]; sprintf(lab,"eT%d-%d",i,j); const string tname(lab);
				Tgrid[i][j].write_ebf(fname,tname);
			}
		}
	}else{
		ebf::Read(fname,"/dJ",&dJ[0],1);
		printf("Reading eTorus grid:\n");
		for(int i=0;i<np;i++){
			printf("%d ",i);
			for(int j=0;j<nr;j++){
				char lab[7];sprintf(lab,"eT%d-%d",i,j); string tname(lab);
				Tgrid[i][j].read_ebf(fname,tname);
			}
		}
	}
	printf("\n");
	Velocity v[9]; Angles thetas[9]; double dens[9];
//Load boundaries of OLR trapping
	strcpy(fname,dirname); strcat(fname,"lindblad.acts");
	FILE *ifile; ifile=fopen(fname,"r");
	int i=0,npO=20; double JOinr[npO],JOinp[npO],x1,x2,x3,x4,x5;
	while(!feof(ifile) && i<npO){
		if(6!=fscanf(ifile,"%f %f %f %f %f %f",&x1,JOinr+i,&x2,&x3,&x4,&x5)) break;
		if(6!=fscanf(ifile,"%f %f %f %f %f %f",&x1,JOinp+i,&x2,&x3,&x4,&x5)) break;
		i++;
	}
	fclose(ifile);
//add extra point to inner edge of trappping zone
	double gr=(JOinp[i-1]-JOinp[i-2])/(JOinr[i-1]-JOinr[i-2]);
	JOinr[i]=0.164; JOinp[i]=JOinp[i-1]+(JOinr[i]-JOinr[i-1])*gr;
	npO=i+1;
	int nmax=50;
	double Jrs[nmax],JCinr[nmax],JCinp[nmax],JCoutr[nmax],JCoutp[nmax];
//load bounaries of CR trapping
	strcpy(fname,dirname); strcat(fname,"corot.acts");
	ifile=fopen(fname,"r"); int ni=0;
	while(!feof(ifile)){
		if(7!=fscanf(ifile,"%f %f %f %f %f %f %f",JCinr+ni,&x1,JCinp+ni,&x2,JCoutp+ni,&x3,&x4)) break;
		JCoutr[ni]=JCinr[ni];
		ni++;
	}
	printf("%d and %d points in Jr read\n",npO,ni); fclose(ifile);
	int nJr=50;
	Jrmax=.12;
	for(int i=0;i<nJr;i++)//Jr values drop from largest value to Jrmin
		Jrs[nJr-i-1]=pow(sqrt(Jrmin)+i*sqrt(Jrmax)/(double)(nJr-1),2);
	int Napo=30;
	double JrPeri[Napo],JpPeri[Napo],JrApo[Napo],JpApo[Napo];
	for(int i=0;i<Napo;i++){
		JpApo[i]=(.1+i*.9/(double)Napo)*JpSun;
		JpPeri[i]=(1.1+i/(double)Napo)*JpSun;
	}
	tabulate_Jr_Jp(Phi,Rzphi[0],JrApo,JpApo,Napo);
	tabulate_Jr_Jp(Phi,Rzphi[0],JrPeri,JpPeri,Napo);
	strcpy(fname,dirname); strcat(fname,"nores8.web");
	FILE *web_file=fopen(fname,"w");
	strcpy(fname,dirname); strcat(fname,"nores.diag");
	FILE *sfile=fopen(fname,"w");
	fprintf(sfile,"       Jr        Jp      Japo    Jperi      JC       JO    psi_min  psi_max\n");
	double J1,J2,Jperi=0,Japo=0;
	double Jbar=fix_Jbar(Jrs,nJr,Tgrid,np,nr,Jgrid,dJ,Rzphi,JCoutp[ni/2],JOinp[npO/2]);
	printf("Jbar = %f\n",Jbar);
	J[1]=.0025;
	for(int i=0;i<nJr;i++){
		J[0]=Jrs[i]; J[2]=Jbar;
		double dJr=2*sqrt(Jrs[i]*Jrmax)/(double)(nJr-1);
		if(J[0]>JCinr[0]){
			J1=lntp(JCinr,JCinp,ni,J[0]); J2=lntp(JCoutr,JCoutp,ni,J[0]);//edges of trapping region
		}else{
			J1=JCinp[0]; J2=JCoutp[0];
		}//Now J2 is Jp(Jr) at top of CR region; J1 isn't used
		double JO1;
		if(J[0]>JOinr[npO-1]) JO1=JOinp[npO-1];
		else if(J[0]<JOinr[0]) JO1=JOinp[0];
		else JO1=lntp(JOinr,JOinp,npO,J[0]);
		//Now JO1 is Jp(Jr) at inside of OLR region
		iTorus T1(J,Tgrid,np,nr,Jgrid,dJ);
		T1.get_crit_Jp(Rzphi,Japo,Jperi,J2,JO1);
		double Japo0=Japo;
		if(Jperi==JO1) Jperi=lntp(JrPeri,JpPeri,Napo,J[0]);
		if(Japo==J2) Japo=lntp(JrApo,JpApo,Napo,J[0]);
		Jbar=.5*(Japo+Jperi);
		//Rzphi can be reached for Japo<Jp<Jperi
		//printf("Japo, Jperi %f %f\n",Japo,Jperi);
		double Delta=MAX(Jbar-Japo,Jperi-Jbar);
		double psi_min,psi_max;
		if(Jbar+Delta>Jperi){
			psi_min=0; psi_max=acos((Jbar-Jperi)/Delta);
		} else {
			psi_min=acos((Jbar-Japo)/Delta); psi_max=PI;
		}
		if(Jbar-Delta*cos(psi_min)<J2) psi_min=acos((Jbar-J2)/Delta);
		if(Jbar-Delta*cos(psi_max)>JO1) psi_max=acos((Jbar-JO1)/Delta);
		int napo,nperi;
		if(i>16 && Japo>J2){
			J[2]=Japo; T1.setJ(J);  napo=T1.InOrbit(Rzphi);
		} else napo=-1;
		if(i>16 && Jperi<JO1){
			J[2]=Jperi; T1.setJ(J); nperi=T1.InOrbit(Rzphi);
		} else nperi=-1;
		fprintf(sfile,"%d %d %f %f %f %f %f %f %f %f\n",
			  napo,nperi,J[0],Japo0,Japo,Jperi,J2,JO1,psi_min,psi_max);
		int nJp=4*nJr-3*i;
		double dpsi=(psi_max-psi_min)/(double)(nJp);
		int jj=0;
		printf("(%6.4f %6.4f %6.4f %6.4f %6.4f)\n",Jrs[i],psi_max,psi_min,Japo,Jperi);
		for(int j=0;j<nJp;j++){
			double psi=psi_min+((double)j+.5)*dpsi;
			J[2]=Jbar-Delta*cos(psi);
			double dJp=Delta*sin(psi)*dpsi;
			iTorus T(J,Tgrid,np,nr,Jgrid,dJ);
			int nv=T.containsPoint(Rzphi,v,thetas,dens,9);
			int k=0;
			if(nv){
				double x[4],y[4],Jac[4];
				for(int m=0; m<nv; m++){
					if(v[m][1]<0) continue;
					x[k]=-Usun-v[m][0]/Units::kms, y[k]=v[m][2]/Units::kms-Vsun;
					Jac[k]=dens[m];
					k++;
				}
				fprintf(web_file,"%d %d %d\n",i,jj,k);
				fprintf(web_file,"%f %f %f %g\n",J[0],J[2],J[2],dJr*dJp);
				for(int m=0;m<k;m++) fprintf(web_file,"%f %f %g\n",x[m],y[m],Jac[m]);
				jj++;
			} else if(j<6)
				fprintf(sfile,"nv=0 %d %d %f\n",i,j,J[2]);
		}
	}
}
