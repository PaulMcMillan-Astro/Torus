/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
//Orbits trapped at a Lindblad resonance
#include <string.h>
#include "resTorus_L.h"
#define SMALL 1e-6

double R0=8,phi0=155,Usun,Vsun;
char dirname[5];

int do_web(Torus **Tgrid,int nr,Actions J,GalaxyPotential *Phi,bar_pot *bar,double Omp,
	   int *resN,double tolJ){
	double x[5],y[5],Jac[5],Jstart=0;
	Actions J0,J1,J2;
	Position Rzphi;
	Rzphi[0]=R0; Rzphi[1]=.01; Rzphi[2]=phi0;
	Velocity Vm[9]; Angles Am[9]; double dm[9]; int nv=0;
	char fname[40];
	strcpy(fname,dirname); strcat(fname,"lindblad.acts");
/*//Uncoment to compute boundaries of trapped region
	FILE *jfile=fopen(fname,"w");
	J[0]=.02;
	while(J[0]<0.19){
		J[0]+=.01;
		Actions Jin2,Jin,Jout;
		resTorus_L rT(Tgrid,nr,J,Phi,bar,Omp,resN,tolJ); J0=rT.get_resJp();
		rT.getJs(J1,J2,1);
		rT.setI(1.01*rT.Imin,-1);
		double Jrin=rT.librationAction();
		Jin[0]=Jrin; Jin[1]=J0[1]; Jin[2]=J0[2]; Jin=rT.prime_it(Jin,1);
		rT.setI(1.01*rT.Imin,1);
		double Jrout=rT.librationAction();
		Jout[0]=Jrout; Jout[1]=J0[1]; Jout[2]=J0[2]; Jout=rT.prime_it(Jout,1);
		rT.setI(1.1*rT.Imin,-1);
		Jrin=rT.librationAction();
		Jin2[0]=Jrin; Jin2[1]=J0[1]; Jin2[2]=J0[2]; Jin2=rT.prime_it(Jin2,1);
		J0=rT.prime_it(J0,1);
		fprintf(jfile,"%f %f %f %f %f %f\n",J0[0],Jin2[0],Jin[0],Jout[0],J1[0],J2[0]);
		fprintf(jfile,"%f %f %f %f %f %f\n",J0[2],Jin2[2],Jin[2],Jout[2],J1[2],J2[2]);
		nv=rT.containsPoint(Rzphi,Vm,Am,dm);
		if(Jstart==0 && nv) Jstart=J[0];
		printf("nv: %d %f\n",nv,J[0]);
	}
	fclose(jfile);
	if(J[0]>0) return 1;*/
//This section computes velocities of visits to Sun
	int nJr=50; double Jrmin=.0002,Jrmax=.12,Jrs[nJr];
	for(int i=0;i<nJr;i++) Jrs[nJr-i-1]=pow(sqrt(Jrmin)+i*sqrt(Jrmax)/(float)(nJr-1),2);
	strcpy(fname,dirname); strcat(fname,"lindblad8.web");
	FILE *web_file=fopen(fname,"w");
	bool done=false;
	for(int j=0;j<nJr;j++){
		if(done) break;
		float dJr=2*sqrt(Jrs[j]*Jrmax)/(float)(nJr-1);
		J[0]=Jrs[j];
		resTorus_L rT(Tgrid,nr,J,Phi,bar,Omp,resN,tolJ);
		Actions resJp=rT.get_resJp();
		int nI=40;
		for(int i=0;i<nI;i++){
			double Ires=.98*rT.Imin+i/(double)(nI-1)*.98*(rT.Imax-rT.Imin);
			rT.setI(Ires,-1); double Jl=rT.librationAction();
			nv=rT.containsPoint(Rzphi,Vm,Am,dm);
			printf("nv: %d %f\n",nv,Jl);
			if(i==0 && nv==0) done=true; if(!nv) break;
			int k=0;
			for(int m=0;m<nv;m++){
				if(Vm[m][1]<0) continue;
				x[k]=-Usun-Vm[m][0]/Units::kms; y[k]=Vm[m][2]/Units::kms-Vsun;
				Jac[k]=dm[m];
				double z=Vm[m][1]/Units::kms;
				printf("(%f %f %f)\n",x[k],y[k],z); k++;
			}
			fprintf(web_file,"%d %d %d\n",j,i,k);
			fprintf(web_file,"%f %f %f %g\n",resJp[0],Jl,resJp[2],dJr);
			for(int m=0;m<k;m++) fprintf(web_file,"%f %f %g\n",x[m],y[m],Jac[m]);
		}
	}
	fclose(web_file);
	return 1;
}
int main(int narg,char **argv){
	if(narg==1){
		printf("You must specify Omega_p\nAlso: ensure that a subdirectory with that value as a name exists\n"); return 0;
	}
	printf("Initialising... Please ensure that you run this from the src directory of the Torus code,\nand that a directory named %s exists\n",argv[1]);


	int omp; sscanf(argv[1],"%d",&omp);
	char fname[40]; sprintf(dirname,"%d/",omp);
	double Omp=0.001*omp;
	int sc[3]={1,0,2}; int3 resN(sc); double PI=acos(-1),torad=PI/180;;
	ifstream Pfile; Pfile.open("../pot/PJM11_best.Tpot");
	if(!Pfile.is_open()) printf("can't open Tpot file\n");
	GalaxyPotential *Phi = new GalaxyPotential(Pfile);
	Pfile.close();
	DiskPar p0=Phi->DiskParameter(0),p1=Phi->DiskParameter(1);
	double Rd=(p0[0]*p0[1]+p1[0]*p1[1])/(p0[0]+p1[0]);//weighted mean of disc scales
	double Rb=0.7*Rd;
	double Ab=3*exp(2)/20.*pow(220./239.,2)*pow(1.5/Rb,3)*.4;
	bar_pot bar(.9,Rb,239./978.,Ab);
	phi0*=torad;
	double R=R0,Jr=.035,Jz=.0025,Jphi=R*240*Units::kms,tolJ=0.003;
	double dPhiR,dPhiz; (*Phi)(R0,0,dPhiR,dPhiz);
	Usun=11; Vsun=12+sqrt(R0*dPhiR)/Units::kms; printf("U,Vsun: %f %f\n",Usun,Vsun);
	Actions J; J[0]=Jr; J[1]=Jz; J[2]=Jphi;
	strcpy(fname,dirname); strcat(fname,"Lindblad.out");
	FILE *ofile; ofile=fopen(fname,"w");
	int nr=5;
	Torus **Tgrid; Tgrid = PJM::matrix<Torus>(nr,2);
//	if(do_web(Tgrid,nr,J,Phi,&bar,Omp,resN,tolJ)) return 0;
	J[0]=.1;
	resTorus_L rT(Tgrid,nr,J,Phi,&bar,Omp,resN,tolJ);
	double Ires=.6*rT.Imin;
	rT.setI(Ires,-1); //rT.check_angle();
	Angles A;
	int NS=20000;
	double x[NS],y[NS],Xmax=15;
	Frequencies Omres=rT.omega();
	double Jl=rT.librationAction();
	printf("Jr Jp, Jl, lOmega: %f %f %f %g\n",J[0],J[2],Jl,Omres[0]);
	double dtr=200*PI/fabs(Omres[2])/(double)NS,dtm=3*PI/fabs(Omres[0]),dt=.5*MIN(dtr,dtm);
	fprintf(ofile,"%f %f %f %f %f %f\n",(NS-1)*dt,J[0],Jl,R0,phi0,Omp);
	for(int i=0;i<NS;i++){
		double t=i*dt;
		A[0]=Omres[0]*t; A[1]=Omres[1]*t; A[2]=Omres[2]*t;
		GCY gcy=rT.FullMap(A);
		x[i]=gcy[0]*cos(gcy[2]); y[i]=gcy[0]*sin(gcy[2]);
		for(int j=0;j<6;j++){
			fprintf(ofile,"%f ",gcy[j]);
		}
		fprintf(ofile,"\n");
	}
	printf("\n");
	ofstream ostrm; ostrm.open("Lindblad.sos");
	bool fit=false;
	if(fit){
		int nI=5;
		for(int i=0;i<nI;i++){
			Ires=1.05*rT.Imin+i/(double)(nI-1)*.99*(rT.Imax-rT.Imin);
			printf("Ires: %f\n",Ires/rT.Imin);
			rT.setI(Ires);
			ostrm << "200" <<'\n';
			rT.SOS(ostrm,200);
		}
		for(int i=0;i<nr;i++){
			Torus T=Tgrid[i][0]*.5;
			for(int j=0;j<1;j++){
				for(int k=0;k<2;k++){
					if(j+k>0) T+=Tgrid[i+j][k]*.5;
				}
			}
			ostrm << "100" <<'\n';
			T.RpSOS(ostrm,100);
		}
	}else{
		Actions J0,J1,J2;
		Torus T;
		J0=rT.get_resJp();
		J0=rT.prime_it(J0,-1);
		T.AutoFit(J0,Phi,tolJ);
		ostrm << "100" <<'\n';
		T.RpSOS(ostrm,100);
		rT.getJs(J1,J2,1);
		T.AutoFit(J1,Phi,tolJ);
		ostrm << "100" <<'\n';
		T.RpSOS(ostrm,100);
		T.AutoFit(J2,Phi,tolJ);
		ostrm << "100" <<'\n';
		T.RpSOS(ostrm,100);
	}
}
