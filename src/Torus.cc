/*******************************************************************************
*                                                                              *
* Torus.cc                                                                     *
*                                                                              *
* C++ code written by Walter Dehnen, 1995,                                     *
*                     Paul McMillan, 2007-                                     *
*                     James Binney, 2015-
* * e-mail: paul@astro.lu.se                                                     *
* github: https://github.com/PaulMcMillan-Astro/Torus                          *
*                                                                              *
*******************************************************************************/

#include <iomanip>
#include "Torus.h"
#include "Point_None.h"
#include "Point_ClosedOrbitCheby.h"
#include "ebf.hpp"
#include "Numerics.h"
#include <cmath>
#define MIN(A,B) ((A)<(B)?(A):(B))

typedef Vector<double,2>   DB2;
typedef Matrix<double,2,2> DB22;
typedef Vector<double,4>   DB4;
typedef Matrix<double,4,4> DB44;

#define MAX(A,B) ((A)>(B)?(A):(B))

////////////////////////////////////////////////////////////////////////////////
// class Torus ************************************************************** //
////////////////////////////////////////////////////////////////////////////////

Torus& Torus::operator= (const Torus& T)
{
	PT=0;
	TM=0;
	J =T.J; E =T.E; Fs=T.Fs; Om=T.Om; dc=T.dc; scale=T.scale;
	if((T.PT)->NumberofParameters()){
		double tmp[(T.PT)->NumberofParameters()];
		(T.PT)->parameters(tmp);
		SetMaps(tmp,
			  (T.TM)->parameters(),
			  (T.GF).parameters(),
			  (T.AM).parameters());
	}  else
		SetMaps((T.TM)->parameters(),
			  (T.GF).parameters(),
			  (T.AM).parameters()); // May cause trouble as PT/TM != 0 initially
	return *this;
}
Torus& Torus::operator*=(const double& x){
	J *=x; E *=x; Fs*=x; Om*=x; dc*=x; Rmin*=x; Rmax*=x; scale*=x;
	vec4 TPtest=TP()*x;
	GenPar SNtest=SN()*x;
	AngPar APtest=AP()*x;
	int k0=PT->NumberofParameters();
	if(k0){
		double tmp0[k0];
		PT->parameters(tmp0);
		int ncx=int(tmp0[2]), ncy=int(tmp0[3+ncx]), ncz=int(tmp0[4+ncx+ncy]);
		//for(int i=0;i<2;i++) tmp0[i]*=x;
		tmp0[1]*=x;//scale r0 but not thmax
		for(int i=0; i<ncx; i++) tmp0[3+i]*=x;
		for(int i=0; i<ncy; i++) tmp0[4+ncx+i]*=x;
		for(int i=0; i<ncz; i++) tmp0[5+ncx+ncy+i]*=x;
		SetMaps(tmp0,TPtest,SNtest,APtest);
	}  else
		SetMaps(TPtest,SNtest,APtest);
	return *this;
}
const Torus Torus::operator* (const double &x){
	Torus T2; T2=*this;
	T2*=x;
	return T2;
}
Torus& Torus::operator+=(const Torus& T){
	J +=T.J; E +=T.E; Fs+=T.Fs; Om+=T.Om; dc+=T.dc; Rmin+=T.Rmin; Rmax+=T.Rmax; scale+=T.scale;
	vec4 TPtest=TP()+T.TP();
	GenPar SNtest=SN()+T.SN();
	AngPar APtest=AP()+T.AP();
	int k=(T.PT)->NumberofParameters(),k0=PT->NumberofParameters(),km=MAX(k,k0);
	if(km){
		double tmp[MAX(k,9)]; vec4 IP;
		if(k) (T.PT)->parameters(tmp);
		else {//an identity map
			IP=T.TM->parameters();
			tmp[0]=Pih;//thmax
			tmp[1]=IP[3]*T.scale;//r0
			tmp[2]=1; tmp[3]=T.scale;//xa const fn
			tmp[4]=2; tmp[5]=0; tmp[6]=T.scale;//ya linear fn
			tmp[7]=1; tmp[8]=T.scale;//za const fn
		}
		int ncx=int(tmp[2]), ncy=int(tmp[3+ncx]), ncz=int(tmp[4+ncx+ncy]);
		double tmp0[MAX(k0,9)];
		if(k0) PT->parameters(tmp0);
		else {//an identity map
			IP=TM->parameters();
			tmp0[0]=Pih;//thmax
			tmp0[1]=IP[3]*scale;
			tmp0[2]=1; tmp0[3]=scale;//xa const fn
			tmp0[4]=2; tmp0[5]=0; tmp0[6]=scale;//ya linear fn
			tmp0[7]=1; tmp0[8]=scale;//za const fn
		}
		int ncx0=int(tmp0[2]), ncy0=int(tmp0[3+ncx0]), ncz0=int(tmp0[4+ncx0+ncy0]);
		int ncxma=MAX(ncx,ncx0),ncyma=MAX(ncy,ncy0),nczma=MAX(ncz,ncz0);
		int kn=7+ncxma+ncyma+nczma;
		double tmp1[kn];
		tmp1[2]=ncxma; tmp1[3+ncxma]=ncyma; tmp1[4+ncxma+ncyma]=nczma;
		if(k!=0 && k0!=0){
			tmp1[0]=MIN(tmp0[0],tmp[0]);//thmax
		}
		else if(k!=0) tmp1[0]=tmp[0];
		else if(k0!=0) tmp1[0]=tmp0[0];//thmax that of only real PT
		tmp1[1]=tmp0[1]+tmp[1];
		for(int i=0; i<ncx0; i++) tmp1[3+i]=tmp0[3+i];
		for(int i=0; i<ncx; i++)  tmp1[3+i]+=tmp[3+i];
		for(int i=0; i<ncy0; i++) tmp1[4+ncxma+i]=tmp0[4+ncx0+i];
		for(int i=0; i<ncy; i++)  tmp1[4+ncxma+i]+=tmp[4+ncx+i];
		for(int i=0; i<ncz0; i++) tmp1[5+ncxma+ncyma+i]=tmp0[5+ncx0+ncy0+i];
		for(int i=0; i<ncz; i++)  tmp1[5+ncxma+ncyma+i]+=tmp[5+ncx+ncy+i];
		SetMaps(tmp1,TPtest,SNtest,APtest);
	}  else
		SetMaps(TPtest,SNtest,APtest);
	return *this;
}
const Torus Torus::operator+ (const Torus& T){
	Torus T2; T2=*this;
	T2+=T;
	return T2;
}
Torus& Torus::operator-=(const Torus& T){
	J -=T.J; E -=T.E; Fs-=T.Fs; Om-=T.Om; dc-=T.dc; Rmin-=T.Rmin; Rmax-=T.Rmax; scale-=T.scale;
	vec4 TPtest=TP()-T.TP();
	GenPar SNtest=SN()-T.SN();
	AngPar APtest=AP()-T.AP();
	int k=(T.PT)->NumberofParameters(),k0=PT->NumberofParameters(),km=MAX(k,k0);
	if(km){
		double tmp[MAX(k,11)];
		if(k) (T.PT)->parameters(tmp);
		else {
			tmp[4]=1; tmp[5]=1;//xa const fn
			tmp[6]=2; tmp[7]=0; tmp[8]=1;//ya linear fn
			tmp[9]=1; tmp[10]=1;//za const fn
		}
		int ncx=int(tmp[4]), ncy=int(tmp[5+ncx]), ncz=int(tmp[6+ncx+ncy]);
		double tmp0[MAX(k0,1)];
		if(k0) PT->parameters(tmp0);
		else {
			tmp0[4]=1; tmp0[5]=1;//xa const fn
			tmp0[6]=2; tmp0[7]=0; tmp0[8]=1;//ya linear fn
			tmp0[9]=1; tmp0[10]=1;//za const fn
		}
		int ncx0=int(tmp0[4]), ncy0=int(tmp0[5+ncx0]), ncz0=int(tmp0[6+ncx0+ncy0]);
		int ncxma=MAX(ncx,ncx0),ncyma=MAX(ncy,ncy0),nczma=MAX(ncz,ncz0);
		int kn=7+ncxma+ncyma+nczma;
		double tmp1[kn];
		tmp1[4]=ncxma; tmp1[5+ncxma]=ncyma; tmp1[6+ncxma+ncyma]=nczma;
		for(int i=0;i<4;i++) tmp1[i]=tmp0[i]-tmp[i];
		for(int i=0; i<ncx0; i++) tmp1[5+i]=tmp0[5+i];
		for(int i=0; i<ncx; i++) tmp1[5+i]-=tmp[5+i];
		for(int i=0; i<ncy0; i++) tmp1[6+ncxma+i]=tmp0[6+ncx0+i];
		for(int i=0; i<ncy; i++) tmp1[6+ncxma+i]-=tmp[6+ncx+i];
		for(int i=0; i<ncz0; i++) tmp1[7+ncxma+ncyma+i]=tmp0[7+ncx0+ncy0+i];
		for(int i=0; i<ncz; i++) tmp1[7+ncxma+ncyma+i]-=tmp[7+ncx+ncy+i];
		SetMaps(tmp1,TPtest,SNtest,APtest);
	}  else
		SetMaps(TPtest,SNtest,APtest);
	return *this;
}
const Torus Torus::operator- (const Torus& T){
	Torus T2; T2=*this;
	T2-=T;
	return T2;
}
void Torus::SetMaps(const double* pp,
	            const vec4 &tp,
	            const GenPar &sn,
	            const AngPar &ap)
{
  delete PT;
  PT = new PoiClosedOrbit(pp);
  if(!TM) TM = new ToyIsochrone;
  TM->set_parameters(tp);
  GF.set_parameters(sn);
  AM.set_parameters(ap);
}
void Torus::SetMaps(const vec4 &tp,
	            const GenPar &sn,
	            const AngPar &ap)
{
  delete PT;
  PT = new PoiNone;
  if(!TM) TM = new ToyIsochrone;
  TM->set_parameters(tp);
  GF.set_parameters(sn);
  AM.set_parameters(ap);
}

void Torus::SetPP(Potential * Phi, const Actions JPT)
{
  delete PT;
  PT = new PoiClosedOrbit(Phi,JPT,TM);
  printf("SetPP: NumberofParameters %d\n",PT->NumberofParameters());
}

void Torus::SetPP(double * param)
{
  delete PT;
  PT = new PoiClosedOrbit(param);
}

void Torus::SetPP(Cheby cb1, Cheby cb2, Cheby cb3, double thmax)
{
  delete PT;
  PT = new PoiClosedOrbit(cb1,cb2,cb3,thmax);
}

void Torus::SetPP()
{
  delete PT;
  PT = new PoiNone;
}

void Torus::SetTP(const vec4& tp)
{
    if(!TM) TM = new ToyIsochrone;
    TM->set_parameters(tp);
}


void Torus::DelMaps()
{
    delete PT;
    delete TM;
}

void Torus::show(ostream& out) const
{
    out <<" Actions                = "<<J<<'\n'
	<<" E, dE                  = "<<E<<' ';
    if(J(0) && J(1))
      out << (  hypot(Om(0),Om(1)) * sqrt(J(0)*J(1)) * dc(0) ) <<'\n';
    else
      out << (  hypot(Om(0),Om(1)) * (J(0) + J(1)) * dc(0) ) <<'\n';
    out <<" Frequencies            = "<<Om<<'\n'
	<<" dJ, chi_rms            = "<<dc<<'\n'
      //<<" parameters of PoiTra   = "<<PP()<<'\n'
	<<" parameters of ToyMap   = "<<TP()<<'\n'
	<<" number of Sn, log|Sn|  : ";
	SN().write_log(out);
    out <<"\n log|dSn/dJr|           : ";
	AM.dSdJ1().write_log(out);
    out <<"\n log|dSn/dJl|           : ";
	AM.dSdJ2().write_log(out);
    out <<"\n log|dSn/dJp|           : ";
	AM.dSdJ3().write_log(out);
}

bool Torus::write_ebf(const string filename, const string torusname, const string mode)
{
	if(mode=="w") {
    ebf::WriteString(filename,"/log","Torus list","w");
  } else if(mode !="a") {
    cerr << "Unknown input mode to write_ebf: " << mode << "\n";
    return false;
  }
  if(ebf::ContainsKey(filename,"/"+torusname)) { return false; }

  int NPTp = PT->NumberofParameters();
  double *PTp;
  if(NPTp) {
    PTp = new double[NPTp];
    PT->parameters(PTp);
  } else {
    PTp = new double[2]; // just to give myself something to delete
  }
  //vec6 PTp = PT->parameters();
  vec4 TMp = TM->parameters();
  int nSn = GF.NumberofParameters();
  int n1_tmp[nSn], n2_tmp[nSn];
  double Sn_tmp[nSn], dSn1[nSn], dSn2[nSn], dSn3[nSn];
  for(int i=0;i!=nSn;i++) {
    n1_tmp[i] = GF.n1(i);
    n2_tmp[i] = GF.n2(i);
    Sn_tmp[i] = GF.coeff(i);
    dSn1[i]   = AM.dSdJ1(i);
    dSn2[i]   = AM.dSdJ2(i);
    dSn3[i]   = AM.dSdJ3(i);
  }

  ebf::Write(filename,"/"+torusname+"/actions",&J[0],"a","",J.NumberofTerms());
  ebf::Write(filename,"/"+torusname+"/df",&Fs,"a","",1);
  ebf::Write(filename,"/"+torusname+"/E",&E,"a","",1);
  ebf::Write(filename,"/"+torusname+"/Rmin",&Rmin,"a","",1);
  ebf::Write(filename,"/"+torusname+"/Rmax",&Rmax,"a","",1);
  ebf::Write(filename,"/"+torusname+"/zmax",&zmax,"a","",1);
  ebf::Write(filename,"/"+torusname+"/frequencies",&Om[0],"a","",Om.NumberofTerms());
  ebf::Write(filename,"/"+torusname+"/errors",&dc[0],"a","",dc.NumberofTerms());

  // Where there is the possibility of alternatives, subdirectories
  if(NPTp)
    ebf::Write(filename,"/"+torusname+"/pointtransform/shellorbit",
	       &PTp[0],"a","",NPTp);
  ebf::Write(filename,"/"+torusname+"/toymap/isochrone",
	     &TMp[0],"a","",TMp.NumberofTerms());
  ebf::Write(filename,"/"+torusname+"/generatingfunction/N1",
	     &n1_tmp[0],"a","",nSn);
  ebf::Write(filename,"/"+torusname+"/generatingfunction/N2",
	     &n2_tmp[0],"a","",nSn);
  ebf::Write(filename,"/"+torusname+"/generatingfunction/Sn",
	     &Sn_tmp[0],"a","",nSn);
  ebf::Write(filename,"/"+torusname+"/anglemap/dSndJ1",&dSn1[0],"a","",nSn);
  ebf::Write(filename,"/"+torusname+"/anglemap/dSndJ2",&dSn2[0],"a","",nSn);
  ebf::Write(filename,"/"+torusname+"/anglemap/dSndJ3",&dSn3[0],"a","",nSn);
  delete[] PTp;
  return true;
}

bool Torus::read_ebf(const string filename, const string torusname) {
  ebf::EbfDataInfo dinfo;
  double *PP;
  vec4 TMp;
  bool got_PT=false;
  // Read data
  if(ebf::ContainsKey(filename,"/"+torusname+"/actions")) {
    // one can allocate memory as y1=(double *)malloc(dinfo.elements*8)
    ebf::Read(filename,"/"+torusname+"/actions", &J[0],J.NumberofTerms());
    ebf::Read(filename,"/"+torusname+"/df",&Fs,1);
    ebf::Read(filename,"/"+torusname+"/E",&E,1);
    ebf::Read(filename,"/"+torusname+"/Rmin",&Rmin,1);
    ebf::Read(filename,"/"+torusname+"/Rmax",&Rmax,1);
    ebf::Read(filename,"/"+torusname+"/zmax",&zmax,1);
    ebf::Read(filename,"/"+torusname+"/frequencies",&Om[0],Om.NumberofTerms());
    ebf::Read(filename,"/"+torusname+"/errors",&dc[0],dc.NumberofTerms());
    if(ebf::ContainsKey(filename,"/"+torusname+"/pointtransform/shellorbit",
			dinfo)) {
      got_PT = true;
      PP = new double[dinfo.elements];
      ebf::Read(filename,"/"+torusname+"/pointtransform/shellorbit",
		&PP[0],dinfo.elements);
    } else PP = new double[2];
    //GenPar::read(int *N1in, int *N2in, double *Snin, short newtot)
    // only one choice of toymap for now, so leave that be
    ebf::Read(filename,"/"+torusname+"/toymap/isochrone",
	      &TMp[0],TMp.NumberofTerms());
    ebf::ContainsKey(filename,"/"+torusname+"/generatingfunction/N1",dinfo);
    int nSn = dinfo.elements;
    int n1_tmp[nSn], n2_tmp[nSn];
    double Sn_tmp[nSn], dSn1[nSn], dSn2[nSn], dSn3[nSn];
    ebf::Read(filename,"/"+torusname+"/generatingfunction/N1",&n1_tmp[0],nSn);
    ebf::Read(filename,"/"+torusname+"/generatingfunction/N2",&n2_tmp[0],nSn);
    ebf::Read(filename,"/"+torusname+"/generatingfunction/Sn",&Sn_tmp[0],nSn);
    ebf::Read(filename,"/"+torusname+"/anglemap/dSndJ1",&dSn1[0],nSn);
    ebf::Read(filename,"/"+torusname+"/anglemap/dSndJ2",&dSn2[0],nSn);
    ebf::Read(filename,"/"+torusname+"/anglemap/dSndJ3",&dSn3[0],nSn);
    GenPar SN,S1,S2,S3;
    SN.read(n1_tmp,n2_tmp,Sn_tmp,nSn);
    S1.read(n1_tmp,n2_tmp,dSn1,nSn);
    S2.read(n1_tmp,n2_tmp,dSn2,nSn);
    S3.read(n1_tmp,n2_tmp,dSn3,nSn);
    AngPar AP(S1,S2,S3);
    if(got_PT)
      SetMaps(PP,TMp,SN,AP);
    else
      SetMaps(TMp,SN,AP);

    delete[] PP;
    return true;
  } else return false;

}


////////////////////////////////////////////////////////////////////////////////
// Deletes all |Sn| < a*|max_Sn|, creates new terms around all |Sn| > b*|max_Sn|
// up to maximum Nmax, sets all |Sn| < off*|max_Sn| to zero
void Torus::TailorAndCutSN(const double ta, const double tb, const double off,
		           const int Nmax)
{
    GenPar SN=GF.parameters();
    SN.tailor(ta,tb,Nmax);
    SN.cut(off);
    GF.set_parameters(SN);
    SN=0.;
    AM.set_parameters(AngPar(SN,SN,SN));
}

////////////////////////////////////////////////////////////////////////////////
void Torus::get_derivs(const Angles &Theta,//input true angles
			 Matrix<double,3,3> &dJdt,//derivs from p-theory
			 Matrix<double,3,3>  &dQdt) const
{
// computes derivatives of coords w.r.t. theta_true
	PSPT QP3,J3(J[0],J[1],J[2],Theta[0],Theta[1],Theta[2]);
	PSPT Jt=AM.Forward3D(J3);//convert to toy
	double dQPdqp[4][4], dQdT3[3][3], dQdJ3[3][3], djdt[2][2];
	PSPT JT3=GF.Forward3DWithDerivs(Jt,djdt); //compute toy actions
	QP3 = PT->Forward3D(TM->Forward3DwithDerivs(JT3,dQdT3));//find indicated point
	PT->Derivatives(dQPdqp);
	TM->Derivatives3D(dQdJ3);
	double dq[3][3],M[3][3]; AM.getMatrix(Jt,M);//M_ij=d\theta_j/d\thetaT_i
	double BB[3],**iM; iM=dmatrix(3,3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) iM[i][j]=M[j][i];
	GaussJordan(iM,3,BB);//now iM=M^-1 i.e. iM_ij = d\thetaT_i/d\theta_j
	Matrix<double,3,3> M1,M2;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M1[i][j]=0;
			for(int k=0;k<3;k++) M1[i][j]+=M[i][k]*dJdt[k][j];
		}//Now M1=M*dJdt
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M2[i][j]=0;
			for(int k=0;k<3;k++) M2[i][j]+=M1[i][k]*M[j][k];
		}//now M2=M*dJdt*M^T
//	for(int i=0;i<3;i++)
//		printf("%f %f %f\n",dJdt[i][0],dJdt[i][1],dJdt[i][2]);
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++) M2[i][j]+=djdt[i][j];
	}//M2 now (d JT/d thetaT)_J'
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) {
			dq[i][j] = dQdT3[i][j];
			for(int k=0;k<3;k++) dq[i][j]+=dQdJ3[i][k]*M2[k][j];
		}
	for(int i=0; i<3; i++){//now deal with point transf
		for(int j=0; j<3; j++){
			if(i<2 && j<2){
				dQdt[i][j]=dQPdqp[i][0]*dq[0][j]+dQPdqp[i][1]*dq[1][j];
			}else{
				dQdt[i][j]=dq[i][j];
			}
		}
	}
	for(int i=0;i<3;i++){//transform derivs so wrt to true angles
		for(int j=0;j<3;j++){
			M1[i][j]=0;
			for(int k=0;k<3;k++) M1[i][j]+=dQdt[i][k]*iM[k][j];
		}
	}
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) dQdt[i][j]=M1[i][j];
	delmatrix(iM,3);
}
void Torus::LevCof3DTrue(const Angles     &Theta,//input true angles
			 Matrix<double,3,3> &dJdt,//derivs from p-theory
			 const Position &Q, //input position
			 PSPT	    &QP3,//output position
			 double         &chsq,//chi^2
			 Vector<double,3>    &B,//2*dchisq/dtheta_i
		     //Matrix<double,3,3>  &A,//approx matrix of 2nd derivs
			 Matrix<double,3,3>  &dQdt) const
{
// computes
// chi^2 = (R0 - R[th])^2 + (z0 - z[th])^2 + R0^2(phi0-phi)^2
// and its derivatives w.r.t. theta_true
	PSPT J3(J[0],J[1],J[2],Theta[0],Theta[1],Theta[2]);
	PSPT Jt=AM.Forward3D(J3);//convert to toy
	double dQPdqp[4][4], dQdT3[3][3], dQdJ3[3][3], djdt[2][2];
	PSPT JT3=GF.Forward3DWithDerivs(Jt,djdt); //compute toy actions
	QP3 = PT->Forward3D(TM->Forward3DwithDerivs(JT3,dQdT3));//find indicated point
	PT->Derivatives(dQPdqp);
	TM->Derivatives3D(dQdJ3);
	double dq[3][3],M[3][3]; AM.getMatrix(Jt,M);//M_ij=d\theta_j/d\thetaT_i
	double BB[3],**iM; iM=dmatrix(3,3);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) iM[i][j]=M[j][i];
	GaussJordan(iM,3,BB);//now iM=M^-1 i.e. iM_ij = d\thetaT_i/d\theta_j
	Matrix<double,3,3> M1,M2;
//Test that iM is the inverse
/*	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			M2[i][j]=0;
			for(int k=0;k<3;k++) M2[i][j]+=iM[i][k]*M[j][k];
			printf("%f ",M2[i][j]);
		}
		printf("\n");
	}*/
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M1[i][j]=0;
			for(int k=0;k<3;k++) M1[i][j]+=M[i][k]*dJdt[k][j];
		}//Now M1=M*dJdt
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M2[i][j]=0;
			for(int k=0;k<3;k++) M2[i][j]+=M1[i][k]*M[j][k];
		}//now M2=M*dJdt*M^T
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++) M2[i][j]+=djdt[i][j];
	}//M2 now (d JT/d thetaT)_J'
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) {
			dq[i][j] = dQdT3[i][j];
			for(int k=0;k<3;k++) dq[i][j]+=dQdJ3[i][k]*M2[k][j];
			//if(j<2) dq[i][j]+= dQdJ3[i][0]*djdt[0][j] + dQdJ3[i][1]*M2[1][j];
		}
	for(int i=0; i<3; i++){//now deal with point transf
		for(int j=0; j<3; j++){
			if(i<2 && j<2){
				dQdt[i][j]=dQPdqp[i][0]*dq[0][j]+dQPdqp[i][1]*dq[1][j];
			}else{
				dQdt[i][j]=dq[i][j];
			}
		}
	}
	for(int i=0;i<3;i++){//transform derivs so wrt to true angles
		for(int j=0;j<3;j++){
			M1[i][j]=0;
			for(int k=0;k<3;k++) M1[i][j]+=dQdt[i][k]*iM[k][j];
		}
	}
	double dR      = (QP3(0)==0.)? 1e99 : (Q(0)-QP3(0));
	double dz      = (Q(1)-QP3(1));
	double dp      = Q(0)*(Q(2)-QP3(2));
	chsq      = (dR*dR+dz*dz+dp*dp);
	B[0]    = dR*M1[0][0] + dz*M1[1][0] + dp*Q(0)*M1[2][0];
	if(!QP3(1) && !QP3(4)) chsq = 1e99;//planar orbit!
	B[1]    = dR*M1[0][1] + dz*M1[1][1] + dp*Q(0)*M1[2][1];
	B[2]    = dR*M1[0][2] + dz*M1[1][2] + dp*Q(0)*M1[2][2];
/*	A[0][0] = pow(M1[0][0],2) + pow(M1[1][0],2)+ pow(Q(0)*M1[2][0],2);
	A[0][1] = M1[0][0]*M1[0][1] + M1[1][0]*M1[1][1] + pow(Q(0),2)*M1[2][0]*M1[2][1];
	A[1][0] = A[0][1];
	A[0][2] = M1[0][0]*M1[0][2] + M1[1][0]*M1[1][2] + pow(Q(0),2)*M1[2][0]*M1[2][2];
	A[2][0] = A[0][2];
	A[1][1] = pow(M1[0][1],2) + pow(M1[1][1],2) + pow(Q(0)*M1[2][1],2);
	A[1][2] = M1[0][1]*M1[0][2] + M1[1][1]*M1[1][2] + pow(Q(0),2)*M1[2][1]*M1[2][2];
	A[2][1] = A[1][2];
	A[2][2] = pow(M1[0][2],2) + pow(M1[1][2],2) + pow(Q(0)*M1[2][2],2);*/
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) dQdt[i][j]=M1[i][j];
	delmatrix(iM,3);
}
void Torus::LevCof3D(const Angles     &theta,//input toy angles
		     Matrix<double,3,3> &dJdt,//derivs from p-theory
		     const Position &Q, //input position
		     PSPT	    &QP3,//output position
		     double         &ch,//chi
		     Vector<double,3>    &B,//2*dchisq/dtheta_i
		     Matrix<double,3,3>  &A,//approx matrix of 2nd derivs
		     Matrix<double,3,3>  &dQdt) const
{
// computes chi, where
// chi^2 = (R0 - R[th])^2 + (z0 - z[th])^2 + R0^2(phi0-phi)^2
// and its derivatives w.r.t. th
	PSPT J3(J[0],J[1],J[2],theta[0],theta[1],theta[2]);
	double dQPdqp[4][4], dQdT3[3][3], dQdJ3[3][3], djdt[2][2];
	PSPT JT3=GF.Forward3DWithDerivs(J3,djdt); //compute toy actions
	//djdt=(dJT/dthetaT)J and only 2x2 because = d^2S/dthetT dthetaT
	QP3 = PT->Forward3D(TM->Forward3DwithDerivs(JT3,dQdT3));//find indicated point
	PT->Derivatives(dQPdqp);//only 2x2
	TM->Derivatives3D(dQdJ3);
	double dq[3][3],M[3][3]; AM.getMatrix(J3,M);
	Matrix<double,3,3> M1,M2;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M1[i][j]=0;
			for(int k=0;k<3;k++) M1[i][j]+=M[i][k]*dJdt[k][j];
		}//Now M1=M*dJdt
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			M2[i][j]=0;
			for(int k=0;k<3;k++) M2[i][j]+=M1[i][k]*M[j][k];
		}//now M2=M*dJdt*M^T
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++) M2[i][j]+=djdt[i][j];
	}//M2 now (d JT/d thetaT)_J'
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) {
			dq[i][j] = dQdT3[i][j];
			for(int k=0;k<3;k++) dq[i][j]+=dQdJ3[i][k]*M2[k][j];
			//if(j<2) dq[i][j]+= dQdJ3[i][0]*djdt[0][j] + dQdJ3[i][1]*M2[1][j];
			//if(j<2) dq[i][j]+= dQdJ3[i][0]*djdt[0][j] + dQdJ3[i][1]*djdt[1][j];
		}
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			if(i<2 && j<2){
				dQdt[i][j]=dQPdqp[i][0]*dq[0][j]+dQPdqp[i][1]*dq[1][j];
			}else{
				dQdt[i][j]=dq[i][j];
			}
		}
	}
	double dR      = (QP3(0)==0.)? 1e99 : (Q(0)-QP3(0));
	double dz      = (Q(1)-QP3(1));
	double dp      = Q(0)*(Q(2)-QP3(2));
	ch      = sqrt(dR*dR+dz*dz+dp*dp);
	B[0]    = dR*dQdt[0][0] + dz*dQdt[1][0] +dp*Q(0)*dQdt[2][0];
	if(!QP3(1) && !QP3(4)) ch = 1e99;//planar orbit!
	B[1]    = dR*dQdt[0][1] + dz*dQdt[1][1] + dp*Q(0)*dQdt[2][1];
	B[2]    = dR*dQdt[0][2] + dz*dQdt[1][2] + dp*Q(0)*dQdt[2][2];
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
void Torus::LevCof(const PSPD     &Jt,
		   const Position &Q,
		   const double   &a,
		   const double   &b,
		   PSPD	          &QP,
		   double         &ch,
		   DB2            &B,
		   DB22           &A,
		   DB22           &dQdt) const
{
// May need a rewrite, to include possibility of effectively a>>1
// computes chi, where
//		chi^2 = a^2 * (R0 - R[th])^2 + b^2 * (z0 - z[th])^2
// and its derivatives w.r.t. th
	register int    i,j,k;
	double dQPdqp[4][4], dqdt[2][2], dqpdj[4][2], djdt[2][2];
	DB22   dq;
	register double dR, dz, aq=a*a,bq=b*b;
	QP = TM->ForwardWithDerivs(GF.ForwardWithDerivs(Jt,djdt),dqdt) >> (*PT);
	PT->Derivatives(dQPdqp);
	TM->Derivatives(dqpdj);
	for(i=0; i<2; i++)
		for(j=0; j<2; j++) {
			dq[i][j] = dqdt[i][j] + dqpdj[i][0]*djdt[0][j] + dqpdj[i][1]*djdt[1][j];
		}
	for(i=0; i<2; i++)
		for(j=0; j<2; j++)
			for(k=0,dQdt[i][j]=0.; k<2; k++)
				dQdt[i][j]+= dQPdqp[i][k] * dq[k][j];

	dR      = (QP(0)==0.)? 1e99 : a*(Q(0)-QP(0));
	dz      = b*(Q(1)-QP(1));
	ch      = hypot(dR,dz);
	B[0]    = a*dR*dQdt[0][0] + b*dz*dQdt[1][0];

	if(!QP(1) && !QP(3)) ch = 1e99;
	B[1]    = a*dR*dQdt[0][1] + b*dz*dQdt[1][1];
	A[0][0] = aq*pow(dQdt[0][0],2) + bq*pow(dQdt[1][0],2);
	A[0][1] = aq*dQdt[0][0]*dQdt[0][1] + bq*dQdt[1][0]*dQdt[1][1];
	A[1][0] = A[0][1];
	A[1][1] = aq*pow(dQdt[0][1],2) + bq*pow(dQdt[1][1],2);
}


////////////////////////////////////////////////////////////////////////////////
void Torus::LevCof(const PSPD         &Jt,
    	           const PSPD         &QP_aim,
		   const double       &a,
		   const double       &b,
		   const double       &vscale,
		   PSPD	              &QP,
		   double             &ch,
 		   DB2                &B,
 		   DB22               &A,
 		   Matrix<double,4,2> &dQPdt) const
{
// May need a rewrite, to include possibility of effectively a>>1
// computes chi, where
//		chi^2 = a^2 * (R0 - R[th])^2 + b^2 * (z0 - z[th])^2
//                      + vscale^2 * ( (vR0 - vR[th])^2 + (vz0 - vz[th])^2 )
// and its derivatives w.r.t. th
    register int    i,j,k;
    double dQPdqp[4][4], dqdt[2][2], dpdt[2][2], dqpdj[4][2], djdt[2][2];
    Matrix<double,4,2>   dqp;
    register double dR, dz, dvR, dvz, aq=a*a, bq=b*b, vscale2=vscale*vscale;


    QP = TM->ForwardWithDerivs(GF.ForwardWithDerivs(Jt,djdt),dqdt,dpdt)>> (*PT);
    PT->Derivatives(dQPdqp);
    TM->Derivatives(dqpdj);
    for(i=0; i<2; i++)
      for(j=0; j<2; j++) {
	dqp[i][j]= dqdt[i][j] + dqpdj[i][0]*djdt[0][j] + dqpdj[i][1]*djdt[1][j];
      }
    for(i=2; i<4; i++)
      for(j=0; j<2; j++) {
	dqp[i][j]= dpdt[i-2][j]+dqpdj[i][0]*djdt[0][j] + dqpdj[i][1]*djdt[1][j];
      }

    for(i=0; i<4; i++)
      for(j=0; j<2; j++)
	for(k=0,dQPdt[i][j]=0.; k<4; k++)
	  dQPdt[i][j]+= dQPdqp[i][k] * dqp[k][j];

    dR      = (QP(0)==0.)? 1e99 : a*(QP_aim(0)-QP(0));
    dz      = b*(QP_aim(1)-QP(1));
    dvR     = vscale*(QP_aim(2)-QP(2));
    dvz     = vscale*(QP_aim(3)-QP(3));
    ch      = sqrt(dR*dR+dz*dz+dvR*dvR+dvz*dvz);
    if(!QP(1) && !QP(3)) ch = 1e99;

    B[0]    = a*dR*dQPdt[0][0] + b*dz*dQPdt[1][0]
      + vscale*dvR*dQPdt[2][0] + vscale*dvz*dQPdt[3][0];
    B[1]    = a*dR*dQPdt[0][1] + b*dz*dQPdt[1][1]
      + vscale*dvR*dQPdt[2][1] + vscale*dvz*dQPdt[3][1];

    A[0][0] = aq*pow(dQPdt[0][0],2) + bq*pow(dQPdt[1][0],2)
      + vscale2*(pow(dQPdt[2][0],2) + pow(dQPdt[3][0],2));
    A[0][1] = aq*dQPdt[0][0]*dQPdt[0][1] + bq*dQPdt[1][0]*dQPdt[1][1]
      + vscale2*(dQPdt[2][0]*dQPdt[2][1] + dQPdt[3][0]*dQPdt[3][1]);
    A[1][0] = A[0][1];
    A[1][1] = aq*pow(dQPdt[0][1],2) + bq*pow(dQPdt[1][1],2)
      + vscale2*(pow(dQPdt[2][1],2) + pow(dQPdt[3][1],2));
}

void Torus::LevCof(const PSPD         &Jt,
		   const PSPD         &QP_aim,
		   const Vector<double,4> &sc,
		   PSPD	              &QP,
		   double             &ch,
		   DB2                &B,
		   DB22               &A,
		   Matrix<double,4,2> &dQPdt) const
{
// May need a rewrite, to include possibility of effectively a>>1
// computes chi, where
//		chi^2 = a^2 * (R0 - R[th])^2 + b^2 * (z0 - z[th])^2
//                      + vscale^2 * ( (vR0 - vR[th])^2 + (vz0 - vz[th])^2 )
// and its derivatives w.r.t. th
	register int    i,j,k;
	double dQPdqp[4][4], dqdt[2][2], dpdt[2][2], dqpdj[4][2], djdt[2][2];
	Matrix<double,4,2>   dqp;
	register double dR, dz, dvR, dvz;//, aq=a*a, bq=b*b, vscale2=vscale*vscale;
	register Vector<double,4> scq;
	for(int i=0;i!=4;i++) scq[i] = sc[i]*sc[i];

	QP = TM->ForwardWithDerivs(GF.ForwardWithDerivs(Jt,djdt),dqdt,dpdt)>> (*PT);
	PT->Derivatives(dQPdqp);
	TM->Derivatives(dqpdj);
	for(i=0; i<2; i++)
		for(j=0; j<2; j++) {
			dqp[i][j]= dqdt[i][j] + dqpdj[i][0]*djdt[0][j] + dqpdj[i][1]*djdt[1][j];
		}
	for(i=2; i<4; i++)
		for(j=0; j<2; j++) {
			dqp[i][j]= dpdt[i-2][j]+dqpdj[i][0]*djdt[0][j] + dqpdj[i][1]*djdt[1][j];
		}

	for(i=0; i<4; i++)
		for(j=0; j<2; j++)
			for(k=0,dQPdt[i][j]=0.; k<4; k++)
				dQPdt[i][j]+= dQPdqp[i][k] * dqp[k][j];

	dR      = (QP(0)==0.)? 1e99 : sc[0]*(QP_aim(0)-QP(0));
	dz      = sc[1]*(QP_aim(1)-QP(1));
	dvR     = sc[2]*(QP_aim(2)-QP(2));
	dvz     = sc[3]*(QP_aim(3)-QP(3));
	ch      = sqrt(dR*dR+dz*dz+dvR*dvR+dvz*dvz);
	if(!QP(1) && !QP(3)) ch = 1e99;

	B[0]    = sc[0]*dR*dQPdt[0][0] + sc[1]*dz*dQPdt[1][0]
		  + sc[2]*dvR*dQPdt[2][0] + sc[3]*dvz*dQPdt[3][0];
	B[1]    = sc[0]*dR*dQPdt[0][1] + sc[1]*dz*dQPdt[1][1]
		  + sc[2]*dvR*dQPdt[2][1] + sc[3]*dvz*dQPdt[3][1];

	A[0][0] = scq[0]*powf(dQPdt[0][0],2) + scq[1]*powf(dQPdt[1][0],2)
		  + scq[2]*powf(dQPdt[2][0],2) + scq[3]*powf(dQPdt[3][0],2);
	A[0][1] = scq[0]*dQPdt[0][0]*dQPdt[0][1] + scq[1]*dQPdt[1][0]*dQPdt[1][1]
		  + scq[2]*dQPdt[2][0]*dQPdt[2][1] + scq[3]*dQPdt[3][0]*dQPdt[3][1];
	A[1][0] = A[0][1];
	A[1][1] = scq[0]*powf(dQPdt[0][1],2) + scq[1]*powf(dQPdt[1][1],2)
		  + scq[2]*powf(dQPdt[2][1],2) + scq[3]*powf(dQPdt[3][1],2);
}





////////////////////////////////////////////////////////////////////////////////
inline bool velocities_are_dependent(		// return: v1, v2 dependent?
				     const double    norm_det,		// input:  z/R*det(D1)/r^2
				     const Velocity& v1,			// input:  v1
				     const Velocity& v2,			// input:  v2
				     const double    tolerance,              // input: tolerance
				     double&         stat)                   // output: discriminant
{
	if( sign(v1(0)) * sign(v1(1)) * sign(v2(0)) * sign(v2(1)) < 0) return false;
	register double x   = hypot(fabs(v1(0))-fabs(v2(0)),
				    fabs(v1(1))-fabs(v2(1))),
	y   = hypot(v1(0),v1(1)),
	eps = tolerance * sqrt(norm_det);
    //if(tolerance<0.05) cerr << x << " " << y << " " <<sqrt(norm_det) << "\n";
	if ( x > eps * y ) return false;
	stat = x/(eps*y);
	return true;
}

bool Torus::containsPoint(
    const Position &Q,       // input:      (R,z,phi)
          Velocity &v1,      // output:	    (vR,vz,vphi)_1
	  DB22     &D1,	     // output:     {d(R,z)/d(T1,T2)}_1
          Angles   &A1,      // output:     T
          Velocity &v2,      // output:	    (vR,vz,vphi)_2
          DB22     &D2,	     // output:     {d(R,z)/d(T1,T2)}_2
          Angles   &A2,      // output:     T
          bool     needA,    // input:      angles out?
          bool     toy,      // input:      toy angles?
          bool     useA,     // input:      use input angles?
          double   delr)     // input:      tolerance in position
    const
// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the
// torus passes through the point given, this happens four times, in each case
// with a different velocity. However, only two of these are independent, since
// changing the sign of both vR and vz simultaneously gives the same orbit. For
// each of these two both possible velocities and the determinant
// | d(x,y,z)/d(Tr,Tl,phi) | is returned. The latter vanishes on the edge of the
// orbit, such that its inverse, the density of the orbit, diverges there
// (that's the reason why the density itself is not returned).
//
// We'll use Levenberg-Marquardt to minimize [R-R(th)]^2 + [z-z(th)]^2
// If the minimum is zero the angles found yield our points (R,z),
// if otherwise the minimum is non-zero (R,z) is never reached by the orbit
{
  // Zeroth: Avoid bugs causing crashes/loops
  if(J(0)<0. || J(1) < 0.) {
    cerr << "Warning: negative actions in containsPoint\n";
    return false;
  }
  // First: Special case of Jl=0. Doing it normally creates infinite loop.
  if(J(1)==0.){
    if(Q(1)!=0.) return false;
    Angles Ang = 0.;
    PSPD tmp = Map(Ang), QP =0., JT, Jt, Jtry;
    register PSPT   QP3D, Jt3D, JT3D;
    double chi,chio,rmin,rmax, dTdt[2][2];;
    DB22   A, Atry, dQdt, dQdtry;
    DB2    B, Btry, dt;
    rmin=tmp(0);
    Ang[0]=Pi;   tmp = Map(Ang);
    rmax=tmp(0);
    if(Q(0)<rmin || Q(0)>rmax) return false;
    // Check that it's in range. If not...
    Jt = PSPD(J(0),J(1),1.,0.);
    const int maxit1=100;
    int it=0;
    const double rtiny=Q(0)*1.e-4;
    while(fabs(QP(0)-Q(0))>rtiny){
      LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
      Jt[2] -= (QP(0)-Q(0))/sqrt(A(0,0));
      while (Jt(2)< 0. ) Jt[2] += TPi;
      while (Jt(2)> TPi) Jt[2] -= TPi;
      if(Jt(2)>Pi) Jt[2] = TPi - Jt(2);
      it++;
      if((rmin-QP(0))*(QP(0)-rmax) < 0.) {
	cerr << "out of range in Newton-Raphson within containsPoint\n"
	     << "rmin=" << rmin << " R="<<QP(0) << "  rmax="<< rmax << '\n';
	return false;
      }
      if(it == maxit1) {
	cerr << "too many iterations in Newton-Raphson within containsPoint\n";
	return false;
      }
    }
    v1[0]    = QP(2);
    v1[1]    = QP(3);
    v1[2]    = J(2)/QP(0);
    if(needA) {
      QP3D.Take_PSPD(QP);
      QP3D[2] = Q(2); QP3D[5] = v1(2);
      Jt3D = QP3D << (*PT) << (*TM); // find theta_phi
      Jt3D.Take_PSPD(Jt);
      Jt3D[2] = J(2); // just in case.

      JT3D       = AM.Backward3DWithDerivs(Jt3D,dTdt);
      JT = JT3D.Give_PSPD();
      if(toy) {A1[0] = Jt3D(3); A1[1] = Jt3D(4); A1[2] = Jt3D(5);}
      else    {A1[0] = JT3D(3); A1[1] = JT3D(4); A1[2] = JT3D(5);}
      A2 = A1;
    } else
      JT = AM.BackwardWithDerivs(Jt,dTdt);

    D1       = 0.;
    D1[0][0] = dQdt(0,0)/dTdt[0][0];
    v2 = v1;
    D2 = 0.;
    D2[0][0] = D1[0][0];
    return true;
  }
//------------------------------------------------------------------------------
// If it isn't  the special case. Do this the hard way.
  const    int    maxit1=100,maxit2=32;
  const    double tiny=1.e-8,
  hit[64]={1/4.,1.,1/8.,15/8.,3/4.,5/4.,0.,1/2.,3/2.,7/4.,
  3/8.,5/8.,7/8.,9/8.,11/8.,13/8.,
  1/16.,3/16.,5/16.,7/16.,9/16.,11/16.,13/16.,15/16.,
  17/16.,19/16.,21/16.,23/16.,25/16.,27/16.,29/16.,
  31/16.,
  1/32.,3/32.,5/32.,7/32.,9/32.,11/32.,13/32.,15/32.,
  17/32.,19/32.,21/32.,23/32.,25/32.,27/32.,29/32.,
  31/32.,33/32.,35/32.,37/32.,39/32.,41/32.,43/32.,
  45/32.,47/32.,49/32.,51/32.,53/32.,55/32.,57/32.,
  59/32.,61/32.,63/32.}; // for finding 2nd theta
  register int    it=0, tried=0;
  DB22   A, Atry, dQdt, dQdtry;                   // for Lev-Mar
  DB2    B, Btry, dt;                             // for Lev-Mar
  double chi, chio, dTdt[2][2];                   // ditto
  PSPD   QP;
  register PSPD   JT, Jt, Jtry;
  register PSPT   QP3D, Jt3D, JT3D;
  register double lam=0.5, lam1, det, JT3_0, det1,
  rq   = Q(0)*Q(0) + Q(1)*Q(1),
  rtin = (delr)? delr : sqrt(rq)*tiny;          // tolerance in position



  Jt = PSPD(J(0),J(1),Pih,0.);                   // guess
  if(useA) { //cerr << A1 << "\n";
	  Jt[2] = A1[0]; Jt[3] = A1[1]; }     // if guess is given
  LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);           // find detc/ detc
  if(std::isnan(B(0))) return false;
  while(chio>rtin && maxit1>it++ && lam < 1.e20 ) {  // Lev Mar iteration
      //cerr << it << " ";
	  lam1  = 1.+lam;
	  det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
	  dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
	  dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
	//cerr << A << B << "\n" << dt << " ";
	  Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
	  if(std::isnan(Jtry(2)) || std::isinf(Jtry(2)) || fabs(Jtry(2))>INT_MAX)
		  Jtry[2] = 0.; // give up
	  if(std::isnan(Jtry(3)) || std::isinf(Jtry(3)) || fabs(Jtry(3))>INT_MAX)
		  Jtry[3] = 0.; // give up
	  if(fabs(Jtry(2))>100.) Jtry[2] -= TPi*int(Jtry(2)*iTPi);
	  if(fabs(Jtry(3))>100.) Jtry[3] -= TPi*int(Jtry(3)*iTPi);
	//cerr << Jtry << "\n";
	  while (Jtry(2)< 0. ) Jtry[2] += TPi;
	  while (Jtry(3)< 0. ) Jtry[3] += TPi;
	  while (Jtry(2)> TPi) Jtry[2] -= TPi;
	  while (Jtry(3)> TPi) Jtry[3] -= TPi;
	//cerr << "here\n";
	  LevCof(Jtry,Q,1.,1.,QP,chi,Btry,Atry,dQdtry);
	  if(chi<chio  && !std::isnan(Btry(0))) {
		  lam *= 0.125;
		  chio = chi;
		  Jt   = Jtry;
		  A    = Atry;
		  B    = Btry;
		  dQdt = dQdtry;
	  } else
		  lam *= 8.;
  }
    //cerr << "\n";
  if(chio > rtin) return false;

  v1[0]    = QP(2);
  v1[1]    = QP(3);
  v1[2]    = J(2)/QP(0);
  if(needA) {
	  QP3D.Take_PSPD(QP);
	  QP3D[2] = Q(2); QP3D[5] = v1(2);
	  Jt3D = QP3D << (*PT) << (*TM);                     // find theta_phi
	  Jt3D.Take_PSPD(Jt); Jt3D[2] = J(2); // just in case
	  JT3D       = AM.Backward3DWithDerivs(Jt3D,dTdt);   // always needed
	  JT = JT3D.Give_PSPD();
	  if(toy) {A1[0] = Jt3D(3); A1[1] = Jt3D(4); A1[2] = Jt3D(5);}
	  else    {A1[0] = JT3D(3); A1[1] = JT3D(4); A1[2] = JT3D(5);}
  } else
	  JT = AM.BackwardWithDerivs(Jt,dTdt);
//Now multiply dQdt by dtdT = dTdt^{-1}
  D1[0][0] = dQdt(0,0)*dTdt[1][1] - dQdt(0,1)*dTdt[1][0];
  D1[0][1] =-dQdt(0,0)*dTdt[0][1] + dQdt(0,1)*dTdt[0][0];
  D1[1][0] = dQdt(1,0)*dTdt[1][1] - dQdt(1,1)*dTdt[1][0];
  D1[1][1] =-dQdt(1,0)*dTdt[0][1] + dQdt(1,1)*dTdt[0][0];
  D1      /= dTdt[0][0]*dTdt[1][1]-dTdt[0][1]*dTdt[1][0];
  det1     = fabs(Q(1)/Q(0)) * fabs(D1(0,0)*D1(1,1)-D1(0,1)*D1(1,0)) / rq;

  JT[2] = TPi-JT(2);
  double Jt2_0 = Jt[2],
  Jt3_0 = Jt[3];
// Try to find other independent velocity.
// It must not fulfill the criterion in the do-while loop

  if(Q(1) == 0.) {	// in symmetry plane, second pair of Vs is dependent.
	  v2 = v1; v2[0] *=-1.;
	  D2 = D1; D2[0][0] *=-1.; D2[0][1] *=-1.;
	  A2 = A1; A2[1] = (A2(1)>Pi)? -Pi+A2(1) : Pi+A2(1);
	  return true;
  }
  bool usedA=false;                  // used supplied guess
  bool notdone=true;
  double depend_tol = 0.1;           // tolerance for whether v are dependent
  double stat, beststat=0;
  PSPD bestJt=0.;
  do {
	  it    = 0;
	  lam   = 0.5;
	  do {
		  if(useA && !usedA) {        // use supplied guess
	    //cerr << A2 << "\n";
			  Jt[2] = A2[0]; Jt[3] = A2[1]; usedA = true;
		  } else {
			  Jt[2] = Jt2_0;           // Or take a shot based on other theta
			  Jt[3] = Jt3_0 + Pi*hit[tried++];
		  }
		  LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
	  } while (QP(0) == 0. && tried<64); // in case of negative actions
	  if(QP(0) == 0.) it = maxit2;       // abort

	  while(chio>rtin && maxit2>it++ && lam < 1.e20 ) {  // Lev-Mar iteration
		  lam1  = 1.+lam;
		  det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
		  dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
		  dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
		  Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
		  while (Jtry(2)< 0. ) Jtry[2] += TPi;
		  while (Jtry(3)< 0. ) Jtry[3] += TPi;
		  while (Jtry(2)> TPi) Jtry[2] -= TPi;
		  while (Jtry(3)> TPi) Jtry[3] -= TPi;
		  LevCof(Jtry,Q,1.,1.,QP,chi,Btry,Atry,dQdtry);

		  if(chi<chio  && !std::isnan(Btry(0))) { // better
			  lam *= 0.125;
			  chio = chi;
			  Jt   = Jtry;
			  A    = Atry;
			  B    = Btry;
			  dQdt = dQdtry;
		  } else                                  // worse
			  lam *= 8.;
	  }
	  v2[0] = QP(2);
	  v2[1] = QP(3);
	  v2[2] = J(2)/Q(0);
	  while(JT(3)>TPi) JT[3]-=TPi;
	  depend_tol = (tried<5)? 0.1 : (tried<10)? 0.05 : 0.01; // tolerance
	  if(chio<=rtin) {
		  notdone = velocities_are_dependent(det1,v1,v2,depend_tol,stat);
		  if(notdone && stat>beststat) {
			  beststat = stat; bestJt = Jt;
		  }
	  }
  } while ( (chio>rtin || notdone ) && 64>tried);

  if(tried>=64 && notdone) {
    //if( chio>rtin || velocities_are_dependent(det1,v1,v2,depend_tol))
	  cerr<<" containsPoint() failed at (R,z)=("<<Q(0)<<","<<Q(1)<<")\n";
	  Jt = bestJt;
	  LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
  }

  if(needA) {
	  QP3D.Take_PSPD(QP);
	  QP3D[2] = Q(2); QP3D[5] = v2(2);
	  Jt3D = QP3D << (*PT) << (*TM);       // find theta_phi
	  Jt3D.Take_PSPD(Jt);  Jt3D[2] = J(2); // just in case.

	  JT3D  = AM.Backward3DWithDerivs(Jt3D,dTdt);
	  JT    = JT3D.Give_PSPD();
	  if(toy) {A2[0] = Jt3D(3); A2[1] = Jt3D(4); A2[2] = Jt3D(5);}
	  else    {A2[0] = JT3D(3); A2[1] = JT3D(4); A2[2] = JT3D(5);}
  } else
	  JT = AM.BackwardWithDerivs(Jt,dTdt);

  D2[0][0] = dQdt(0,0)*dTdt[1][1] - dQdt(0,1)*dTdt[1][0];
  D2[0][1] =-dQdt(0,0)*dTdt[0][1] + dQdt(0,1)*dTdt[0][0];
  D2[1][0] = dQdt(1,0)*dTdt[1][1] - dQdt(1,1)*dTdt[1][0];
  D2[1][1] =-dQdt(1,0)*dTdt[0][1] + dQdt(1,1)*dTdt[0][0];
  D2      /= (dTdt[0][0]*dTdt[1][1]-dTdt[0][1]*dTdt[1][0]);

    //if(tried>=64) cerr << D2;

  if(chio < rtin) return true; else return false;
}

void Torus::CheckLevCof(PSPD QP_aim, Angles A_in) {
  PSPD Jt = PSPD(J(0),J(1),A_in[0],A_in[1]), Jtry,QP,oQP;
  double small = 1.e-4;
  double chi,chio,sc;
  DB22   A;
  Matrix<double,4,2> dQPdt;
  DB2    B;
  double tmpx2 = powf(Rmax-Rmin,2) + powf(zmax,2), tmpv2;
  Angles tmpA = 0.; tmpA[0] = Pih;
  //Position tmpq;
  //PSPT   tmp3Dqp = MapfromToy3D(A_in);
  PSPD   tmpqp = MapfromToy(tmpA);
  tmpv2 = tmpqp(2)*tmpqp(2) + tmpqp(3)*tmpqp(3);
  sc = sqrt(tmpx2/tmpv2);

  LevCof(Jt,QP_aim,1.,1.,sc,oQP,chio,B,A,dQPdt);
  //cerr << dQPdt;
  Jtry = Jt; Jtry[2] = A_in[0] + small;

  LevCof(Jtry,QP_aim,1.,1.,sc,QP,chi,B,A,dQPdt);
  //for(int i=0;i!=4;i++) cerr << (QP[i]-oQP[i])/small << ' ';
  //cerr << '\n';
  cerr << -2*B(0) << ' ' << (chi*chi-chio*chio)/small << '\n';

  Jtry = Jt; Jtry[3] = A_in[1] + small;
  LevCof(Jtry,QP_aim,1.,1.,sc,QP,chi,B,A,dQPdt);
  //for(int i=0;i!=4;i++) cerr << (QP[i]-oQP[i])/small << ' ';
  //cerr << '\n';

  cerr << -2*B(1) << ' ' << (chi*chi-chio*chio)/small << '\n';

}

void Torus::CheckLevCof(Position Q_aim, Angles A_in) {
  PSPD Jt = PSPD(J(0),J(1),A_in[0],A_in[1]), Jtry,QP,oQP;
  double small = 1.e-5;
  double chi,chio,sc;
  DB22   A ;
  Matrix<double,2,2> dQdt;
  DB2    B;

  LevCof(Jt,Q_aim,1.,1.,oQP,chio,B,A,dQdt);
  //cerr << dQdt(0,0) << ' ' << dQdt(1,0) << ' ';
  Jtry = Jt; Jtry[2] = A_in[0] + small;
  LevCof(Jtry,Q_aim,1.,1.,QP,chi,B,A,dQdt);
  cerr << -2*B(0) << ' '
       << (chi*chi-chio*chio)/small << '\n';
  //<< (QP(0)-oQP(0))/small << ' ' <<  (QP(1)-oQP(1))/small << '\n';
  Jtry = Jt; Jtry[3] = A_in[1] + small;
  LevCof(Jtry,Q_aim,1.,1.,QP,chi,B,A,dQdt);
  cerr << -2*B(1) << ' ' << (chi*chi-chio*chio)/small << '\n';

}

////////////////////////////////////////////////////////////////////////////////
DB2 Torus::DistancetoPSP(const PSPD &QP_aim, double &scale, Angles &Aclosest) const
// We'll use Levenberg-Marquardt to minimize
// [R-R(th)]^2 + [z-z(th)]^2 + scale^2 * ( [vR-vR(th)]^2 + [vz-vz(th)]^2 )
{
  const    int    maxit=100;
  const    double tiny=1.e-8;
  register int    it=0;
           DB22   A, Atry;
	   Matrix<double,4,2> dQPdt, dQPdtry;
           DB2    B, Btry, dt;
           double chi, chio=1.e99;
  register double lam=0.5, lam1, det, r=hypot(QP_aim(0),QP_aim(1)), rtin=r*tiny;
  PSPD   QP, QP_best;
  register PSPD   Jt, Jtry;
  register double sc;
  DB2 out=0.;
  Angles Astart;
  if(scale==0.) {
    //cerr << Rmin << ' ' << Rmax << ' ' << zmax << '\n';
    double tmpx2 = powf(0.5*(Rmax-Rmin),2) + powf(zmax,2), tmpv2;
    Angles tmpA = 0.; tmpA[0] = Pih;
    PSPD   tmpqp = MapfromToy(tmpA);
    //cerr << tmpqp << ' ' <<  Rmin << ' ' << Rmax << ' ' << zmax << '\n';
    tmpv2 = tmpqp(2)*tmpqp(2) + tmpqp(3)*tmpqp(3);
    sc = sqrt(tmpx2/tmpv2);
    scale = sc;
    //cerr << sc << '\n';
  } else sc = scale;


  // Avoid bugs causing crashes/loops
  if(J(0)<0. || J(1) < 0.) {
    cerr << "Warning: negative actions in DistancetoPoint\n";
    out = 0.;
    return out;
  }

  // find a starting point using a course grid
  const int ngridr = 30, ngridz=30;
  Angles Atest;
  for(int i=0;i!=ngridr;i++) {
    Atest[0] = TPi*i/double(ngridr);
    for(int j=0;j!=ngridz;j++) {
      Atest[1] = TPi*j/double(ngridz);
      QP = MapfromToy(Atest);
      chi = pow(QP(0)-QP_aim(0),2) +pow(QP(1)-QP_aim(1),2) +
	sc*sc*(pow(QP(2)-QP_aim(2),2) +pow(QP(3)-QP_aim(3),2));
      if(chi<chio) {
	Astart = Atest;
	chio = chi;
      }
    }
    //std::cout << '\n';
  }

  //Astart = 1.;

  Jt = PSPD(J(0),J(1),Astart[0],Astart[1]);
  LevCof(Jt,QP_aim,1.,1.,sc,QP,chio,B,A,dQPdt);
  if(std::isnan(B(0))) {
    out = 0.; return out;
  }
  QP_best = QP;
  //cerr << Jt << '\n';
  //cerr << QP_aim << ' ' << QP_best << '\n';
  while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
    lam1  = 1.+lam;
    det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
    dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
    dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
    Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
    //cerr << Jtry << ' ' << chio << '\n';
    while (Jtry(2)< 0. ) Jtry[2] += TPi;
    while (Jtry(3)< 0. ) Jtry[3] += TPi;
    while (Jtry(2)> TPi) Jtry[2] -= TPi;
    while (Jtry(3)> TPi) Jtry[3] -= TPi;
    LevCof(Jtry,QP_aim,1.,1.,sc,QP,chi,Btry,Atry,dQPdtry);
    //cerr << chi << ' ' << Jtry[2] << ' ' << Jtry[3] << '\n';
    if(chi<chio  && !std::isnan(Btry(0))) {
      lam *= 0.125;
      chio = chi;
      Jt   = Jtry;
      A    = Atry;
      B    = Btry;
      QP_best = QP;
      dQPdt = dQPdtry;
    } else
      lam *= 8.;
  }
  Aclosest[0] = Jt[2]; Aclosest[1] = Jt[3];
  //cerr << QP_aim << ' ' << QP_best << '\n';
  if(chio < rtin) { out = 0.; return out; }
  else {
    out[0] = hypot(QP_best(0)-QP_aim(0),QP_best(1)-QP_aim(1));
    out[1] = hypot(QP_best(2)-QP_aim(2),QP_best(3)-QP_aim(3));
    return out;
  }
}


////////////////////////////////////////////////////////////////////////////////
Vector<double,4> Torus::DistancetoPSP(const PSPD &QP_aim,
				      Vector<double,4> &scales,
				      Angles &Aclosest) const
// We'll use Levenberg-Marquardt to minimize
// (scale[0]*[R-R(th)])^2   + (scale[1]*[z-z(th)])^2 +
// (scale[2]*[vR-vR(th)])^2 + (scale[3]*[vz-vz(th)])^2
{
  const    int    maxit=100;
  const    double tiny=1.e-8;
  register int    it=0;
           DB22   A, Atry;
	   Matrix<double,4,2> dQPdt, dQPdtry;
           DB2    B, Btry, dt;
           double chi, chio=1.e99;
  register double lam=0.5, lam1, det, r=hypot(QP_aim(0),QP_aim(1)), rtin=r*tiny;
  PSPD   QP, QP_best;
  register PSPD   Jt, Jtry;
  register Vector<double,4> sc;
  Vector<double,4>  out=0.;
  Angles Astart;
  if(scales==0.) {
    //FindLimits();
    //cerr << Rmin << ' ' << Rmax << ' ' << zmax << '\n';
    Angles tmpA = 0.; tmpA[0] = Pih;
    PSPD   tmpqp = MapfromToy(tmpA);
    double tmpR2 = powf(0.5*(Rmax-Rmin),2), tmpz2 = zmax*zmax,
      tmpvR2 = tmpqp(2)*tmpqp(2), tmpvz2 =tmpqp(3)*tmpqp(3);
    //cerr << tmpqp << ' ' <<  Rmin << ' ' << Rmax << ' ' << zmax << '\n';
    sc[0] = 1./sqrt(tmpR2);
    sc[1] = 1./sqrt(tmpz2);
    sc[2] = 1./sqrt(tmpvR2);
    sc[3] = 1./sqrt(tmpvz2);
      //sc = sqrt(tmpR2/tmpv2);
    scales = sc;
    //cerr << sc << '\n';
  } else sc = scales;


  // Avoid bugs causing crashes/loops
  if(J(0)<0. || J(1) < 0.) {
    cerr << "Warning: negative actions in DistancetoPoint\n";
    out = 0.;
    return out;
  }

  // find a starting point using a course grid
  const int ngridr = 30, ngridz=30;
  Angles Atest;
  for(int i=0;i!=ngridr;i++) {
    Atest[0] = TPi*i/double(ngridr);
    for(int j=0;j!=ngridz;j++) {
      Atest[1] = TPi*j/double(ngridz);
      QP = MapfromToy(Atest);
      chi = sc[0]*sc[0]*powf(QP(0)-QP_aim(0),2) +
	sc[1]*sc[1]*powf(QP(1)-QP_aim(1),2) +
	sc[2]*sc[2]*powf(QP(2)-QP_aim(2),2) +
	sc[3]*sc[3]*powf(QP(3)-QP_aim(3),2);
      if(chi<chio) {
	Astart = Atest;
	chio = chi;
      }
      //if(chi<0.1)
      //std::cout << chi << ' ';
	//else cerr << "0.1 ";
    }
    //std::cout << '\n';
  }

  //Astart = 1.;

  Jt = PSPD(J(0),J(1),Astart[0],Astart[1]);
  LevCof(Jt,QP_aim,sc,QP,chio,B,A,dQPdt);
  if(std::isnan(B(0))) {
    out = 0.; return out;
  }
  QP_best = QP;
  //cerr << Jt << '\n';
  //cerr << QP_aim << ' ' << QP_best << '\n';
  while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
    lam1  = 1.+lam;
    det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
    dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
    dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
    Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
    //cerr << Jtry << ' ' << chio << '\n';
    while (Jtry(2)< 0. ) Jtry[2] += TPi;
    while (Jtry(3)< 0. ) Jtry[3] += TPi;
    while (Jtry(2)> TPi) Jtry[2] -= TPi;
    while (Jtry(3)> TPi) Jtry[3] -= TPi;
    LevCof(Jtry,QP_aim,1.,1.,sc[0],QP,chi,Btry,Atry,dQPdtry);
    //cerr << chi << ' ' << Jtry[2] << ' ' << Jtry[3] << '\n';
    if(chi<chio  && !std::isnan(Btry(0))) {
      lam *= 0.125;
      chio = chi;
      Jt   = Jtry;
      A    = Atry;
      B    = Btry;
      QP_best = QP;
      dQPdt = dQPdtry;
    } else
      lam *= 8.;
  }
  Aclosest[0] = Jt[2]; Aclosest[1] = Jt[3];
  //cerr << QP_aim << ' ' << QP_best << '\n';
  if(chio < rtin) { out = 0.; return out; }
  else {
    out[0] = hypot(QP_best(0)-QP_aim(0),QP_best(1)-QP_aim(1));
    out[1] = hypot(QP_best(2)-QP_aim(2),QP_best(3)-QP_aim(3));
    return out;
  }
}







////////////////////////////////////////////////////////////////////////////////
double Torus::DistancetoPoint(const Position &Q) const
// We'll use Levenberg-Marquardt to minimize [R-R(th)]^2 + [z-z(th)]^2
{
    const    int    maxit=100;
    const    double tiny=1.e-8;
    register int    it=0;
	     DB22   A, Atry, dQdt, dQdtry;
	     DB2    B, Btry, dt;
	     double chi, chio;
    register double lam=0.5, lam1, det, r=hypot(Q(0),Q(1)), rtin=r*tiny;
	     PSPD   QP;
    register PSPD   Jt, Jtry;

  // Avoid bugs causing crashes/loops
  if(J(0)<0. || J(1) < 0.) {
    cerr << "Warning: negative actions in DistancetoPoint\n";
    return 0;
  }

    Jt = PSPD(J(0),J(1),Pih,0.);
    LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
    if(std::isnan(B(0))) return 0;

    while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
	lam1  = 1.+lam;
	det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
	dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
	dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
	Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
        LevCof(Jtry,Q,1.,1.,QP,chi,Btry,Atry,dQdtry);
	if(chi<chio  && !std::isnan(Btry(0))) {
	    lam *= 0.125;
	    chio = chi;
	    Jt   = Jtry;
            A    = Atry;
	    B    = Btry;
	    dQdt = dQdtry;
	} else
	    lam *= 8.;
    }

    return (chio > rtin)? chio : 0.;
}


////////////////////////////////////////////////////////////////////////////////
double Torus::DistancetoPoint(const Position &Q, double &thr, double &thz) const
// We'll use Levenberg-Marquardt to minimize [R-R(th)]^2 + [z-z(th)]^2
{
	const    int    maxit=100;
	const    double tiny=1.e-8;
	register int    it=0, itb=0;
	DB22   A, Atry, dQdt, dQdtry;
	DB2    B, Btry, dt;
	double chi, chio, chih;
	register double lam=0.5, lam1, det, r=hypot(Q(0),Q(1)), rtin=r*tiny;
	PSPD   QP;
	register PSPD   Jt, Jtry;

  // Avoid bugs causing crashes/loops
	if(J(0)<0. || J(1) < 0.) {
		cerr << "Warning: negative actions in DistancetoPoint\n";
		return 0;
	}
	if(thr>TPi || thz>TPi || thr<0. || thz<0. ) { thr=Pih; thz=0.; }
	Jt = PSPD(J(0),J(1),Pih,0.);
  //Jt = PSPD(J(0),J(1),thr,thz);
	LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
	if(std::isnan(B(0))) return 0;

	while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
		if(it==30) chih = chio;
		lam1  = 1.+lam;
		det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
		dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
		dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
		Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
		LevCof(Jtry,Q,1.,1.,QP,chi,Btry,Atry,dQdtry);
		if(chi<chio  && !std::isnan(Btry(0))) {
			itb  = it;
			lam *= 0.125;
			chio = chi;
			Jt   = Jtry;
			A    = Atry;
			B    = Btry;
			dQdt = dQdtry;
		} else
			lam *= 8.;
	}
	thr = Jt[2]; thz = Jt[3];
	return (chio > rtin)? chio : 0.;
}
////////////////////////////////////////////////////////////////////////////////
double Torus::DistancetoPoint_Ang(const Position &Q, Angles &theta) const
// Finds true angles theta of a point of approach to Position Q
{
	double thr=Pih, thz=0;
	const    int    maxit=100;
	const    double tiny=1.e-8;
	register int    it=0, itb=0;
	DB22   A, Atry, dQdt, dQdtry;
	DB2    B, Btry, dt;
	double chi, chio, chih;
	register double lam=0.5, lam1, det, r=hypot(Q(0),Q(1)), rtin=r*tiny;
	PSPD   QP;
	register PSPD   Jt, Jtry;

  // Avoid bugs causing crashes/loops
	if(J(0)<0. || J(1) < 0.) {
		cerr << "Warning: negative actions in DistancetoPoint\n";
		return 0;
	}
	Jt = PSPD(J(0),J(1),Pih,0.);
	LevCof(Jt,Q,1.,1.,QP,chio,B,A,dQdt);
	if(std::isnan(B(0))) return 0;

	while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
		if(it==30) chih = chio;
		lam1  = 1.+lam;
		det   = A(0,0)*A(1,1)*pow(lam1,2) - pow(A(0,1),2);
		dt[0] = (B(0)*lam1*A(1,1)-B(1)*A(0,1)) / det;
		dt[1] = (B(1)*lam1*A(0,0)-B(0)*A(0,1)) / det;
		Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
		LevCof(Jtry,Q,1.,1.,QP,chi,Btry,Atry,dQdtry);
		if(chi<chio  && !std::isnan(Btry(0))) {
			itb  = it;
			lam *= 0.125;
			chio = chi;
			Jt   = Jtry;
			A    = Atry;
			B    = Btry;
			dQdt = dQdtry;
		} else
			lam *= 8.;
	}
	PSPT QP3D,Jt3D,JT3D;
	QP3D.Take_PSPD(QP);
	QP3D[2] = Q(2); QP3D[5] = J(2)/QP(0);
	Jt3D = QP3D << (*PT) << (*TM); //toy actions and angles
	Jt3D.Take_PSPD(Jt); Jt3D[2] = J(2); // now true actions and toy angles
	JT3D = AM.Backward3D(Jt3D);   // now true actions and angles
	theta[0] = JT3D(3); theta[1] = JT3D(4); theta[2] = JT3D(5);
	//GCY gcy=FullMap(theta);
	//printf("Torus: %f %f %f %f %g %g\n",gcy[0],QP[0],gcy[1],QP[1],chio,rtin);
	return (chio > rtin)? sqrt(chio) : 0.;
	//return sqrt(chio);
}


////////////////////////////////////////////////////////////////////////////////
double Torus::DistancetoRadius(const double R) const
// We'll use Levenberg-Marquardt to minimize [R-R(th)]^2
{
    const    int      maxit=1000;
    const    double   tiny=1.e-4;
    register int      it=0;
	     DB22     A, Atry, dQdt, dQdtry;
	     DB2      B, Btry, dt;
	     double   chi, chio;
    register double   lam=0.5, lam1, rtin=R*tiny;
	     PSPD     QP;
	     Position Q=0.;
    register PSPD     Jt, Jtry;
             Angles Ang = 0.;
	     double rmin,rmax;

  // Avoid bugs causing crashes/loops
  if(J(0)<0. || J(1) < 0.) {
    cerr << "Warning: negative actions in DistancetoRadius\n";
    return 0;
  }
    Ang[1] = Pih;
    rmin = (Map(Ang))(0);  // Questionable, but certainly
    Ang[0] = Pi; Ang[1] = 0.;
    rmax = (Map(Ang))(0);
    if(R<=rmax && R>=rmin) return 0.; // Added for a quicker, easier thing

    Q[0] = R;
    Jt   = PSPD(J(0),J(1),Pih,0.);
    LevCof(Jt,Q,1.,0.,QP,chio,B,A,dQdt);
    while(chio>rtin && maxit>it++ && lam < 1.e20 ) {
	lam1  = 1.+lam;
	dt[0] = B(0)/A(0,0)/lam1;
	dt[1] = B(1)/A(1,1)/lam1;
	Jtry  = PSPD(J(0),J(1),Jt(2)+dt(0),Jt(3)+dt(1));
	LevCof(Jtry,Q,1.,0.,QP,chi,Btry,Atry,dQdtry);
	if(chi<chio && !(std::isnan(Btry[0])) ) {
	  lam *= 0.125;
	  chio = chi;
	  Jt   = Jtry;
	  A    = Atry;
	  B    = Btry;
	  dQdt = dQdtry;
	} else
	  lam *= 8.;
    }

    return (chio > rtin)? chio : 0.;
}





////////////////////////////////////////////////////////////////////////////////
void Torus::SOSroot(const double t2, double& z, double& dz) const//t2 is toy angle
{
	double          dQPdqp[4][4], dqdt[2][2], dqpdj[4][2], djdt[2][2];
	register PSPD   jt, QP;
	register int    i,j;
	register double dqdt2;
	PSPD Jtroot=PSPD(Jtr[0],Jtr[1],Jtr[2],t2);
//	Jtroot[3] = t2;
	jt        = GF.ForwardWithDerivs(Jtroot, djdt);//get toy actions
	QP        = TM->ForwardWithDerivs(jt, dqdt) >> (*PT);//dqdt derivs wrt toy angles
	PT->Derivatives(dQPdqp);
	TM->Derivatives(dqpdj);//dqpdj derivs wrt toy actions
	z  = QP(1);
	dz = 0.;
	for(i=0; i<2; i++) {
		dqdt2 = dqdt[i][1];//deriv of r or theta wrt to toy theta_z
		for(j=0; j<2; j++)
			dqdt2 += dqpdj[i][j] * djdt[j][1];
		dz += dqdt2 * dQPdqp[1][i];
	}
}
////////////////////////////////////////////////////////////////////////////////
// Files (R,z,pR,pz,theta_z) where torus cuts z=0
////////////////////////////////////////////////////////////////////////////////
void Torus::SOS(ostream& to, const int Nthr)
// runs _toy_ theta_r from 0 - Pi finding _toy_ theta_z that brings
// orbit to z=0
{
//	Jtroot = PSPD(J(0),J(1),0.,0.);
	Jtr[0]=J[0]; Jtr[1]=J[1]; Jtr[2]=Jtr[3]=0;
	PSPD Jtroot=PSPD(J(0),J(1),0.,0.);
	to << (Jtroot>>GF>>(*TM)>>(*PT)) <<"    "<<Jtroot(3)<<'\n';
	const double allowed = -1.e-8;
	double tempz1,tempz2,tempdz;
	for(int ithr=1; ithr<Nthr; ithr++) {
		Jtr[2]=Jtroot[2] = double(ithr) * Pi / double(Nthr);
		SOSroot(-Pih,tempz1,tempdz); SOSroot(Pih,tempz2,tempdz);
		if(tempz1*tempz2<0.){
			Jtroot[3] = rtsafe(this,&Torus::SOSroot, -Pih, -(-Pih), -allowed);
			to << (Jtroot>>GF>>(*TM)>>(*PT)) <<"    "<<Jtroot(3)<<'\n';
		}
	}
	Jtroot[2] = Pi;
	Jtroot[3] = 0.;
	to << (Jtroot>>GF>>(*TM)>>(*PT)) <<"    "<<Jtroot(3)<<'\n';
}

////////////////////////////////////////////////////////////////////////////////
// Called by resSOS to compute z and dz/dtheta_z subject to theta_r-theta_z=t1p
////////////////////////////////////////////////////////////////////////////////
void Torus::resSOSroot(const double t2, double& z, double& dz) const
{
	double          dQPdqp[4][4], dqdt[2][2], dqpdj[4][2], djdt[2][2];
	register PSPD   JT,jt, QP;
	register int    i,j;
	register double t1=t2+t1p,dqdt2[2],dz2[2];
	//t1,t2 are the true thetas
//	Jtroot[2]=t1; Jtroot[3] = t2;
	PSPD Jtroot=PSPD(Jtr[0],Jtr[1],t1,t2);
	JT	  = AM.Forward(Jtroot);//get toy angles
	jt        = GF.ForwardWithDerivs(JT, djdt);//get toy actions & derivs wrt toy angles
	QP        = TM->ForwardWithDerivs(jt, dqdt) >> (*PT);//dqdt derivs q wrt toy angles
	PT->Derivatives(dQPdqp);
	TM->Derivatives(dqpdj);//derivs qp wrt toy actions
	z  = QP(1);
	for(int k=0;k<2;k++){
		dz2[k] = 0.;//deriv of z wrt tr and tz
		for(i=0; i<2; i++) {
			dqdt2[k] = dqdt[i][k];
			for(j=0; j<2; j++)
				dqdt2[k] += dqpdj[i][j] * djdt[j][k];
			dz2[k] += dqdt2[k] * dQPdqp[1][i];
		}
	}
	dz=dz2[1]+dz2[0];//for N=(1,-1,0)
}
////////////////////////////////////////////////////////////////////////////////
//Finds where the torus cuts z=0 subject to the restricton theta_r-theta_z=t1pp
////////////////////////////////////////////////////////////////////////////////
void Torus::resSOS(ostream& to,const double t1pp)
//for given difference of _true_ theta_r - theta_z find _true_ theta_r
//and theta_z that brings orbit to z=0
{
	t1p=t1pp;
	PSPD Jtroot = PSPD(J(0),J(1),0.,0.);
	Jtr[0]=J[0]; Jtr[1]=J[1];
	const double allowed = -1.e-8;
	double tempz1,tempz2,tempdz;
	resSOSroot(-Pih,tempz1,tempdz); resSOSroot(Pih,tempz2,tempdz);
	if(tempz1*tempz2<0.){
		Jtr[3]=Jtroot[3] = rtsafe(this,&Torus::resSOSroot, -Pih, -(-Pih), -allowed);
		Jtr[2]=Jtroot[2]=t1p+Jtroot[3];
		to << std::setw(12) << (Jtroot>>AM>>GF>>(*TM)>>(*PT)) <<"   "<<Jtroot(3)<<'\n';
	}
}
////////////////////////////////////////////////////////////////////////////////
void Torus::SOS_z_root(const double t2, double& RmRS, double& dR) const
{
    double          dQPdqp[4][4], dqdt[2][2], dqpdj[4][2], djdt[2][2];
    register PSPD   jt, QP;
    register int    i,j;
    register double dqdt2;
    PSPD Jtroot=PSPD(Jtr[0],Jtr[1],t2,Jtr[3]);
    jt        = GF.ForwardWithDerivs(Jtroot, djdt);
    QP        = TM->ForwardWithDerivs(jt, dqdt) >> (*PT);
    PT->Derivatives(dQPdqp);
    TM->Derivatives(dqpdj);
    RmRS  = QP(0)-RforSOS;
    dR = 0.;
    for(i=0; i<2; i++) {
	dqdt2 = dqdt[i][0];
	for(j=0; j<2; j++)
	    dqdt2 += dqpdj[i][j] * djdt[j][0];
        dR += dqdt2 * dQPdqp[0][i];
    }
}
////////////////////////////////////////////////////////////////////////////////
void Torus::SOS_z(ostream& to, const double RSOS, const int Nthz)
{
	PSPD Jtroot = PSPD(J(0),J(1),0.,0.);
	Jtr[0]=J[0]; Jtr[1]=J(1);
  RforSOS = RSOS;
  int ntry=100;
  double tmpR1,tmpR2,thmin,thmax;
  const double allowed = -1.e-8;
  // do this slowly but safely-ish
  for(int ithz=0; ithz<Nthz; ithz++) {
    Jtr[3]=Jtroot[3] = double(ithz) * TPi / double(Nthz);
    // OK first, can I do this the easy way?
    Jtroot[2] = thmin = 0.; tmpR1 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    Jtroot[2] = thmax = Pi; tmpR2 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    for(int i=0;i!=ntry && (tmpR1-RSOS)*(tmpR2-RSOS)>0;i++) {
      if(i) { thmin = thmax; tmpR1 = tmpR2; }
      Jtroot[2] = thmax = double(i+1)*TPi/double(ntry);
      tmpR2 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    }
    if((tmpR1-RSOS)*(tmpR2-RSOS)<0) {
      Jtroot[2] = rtsafe(this,&Torus::SOS_z_root, thmin, thmax, -allowed);
      to <<  (Jtroot>>GF>>(*TM)>>(*PT)) <<"    "<<Jtroot(3)<<'\n';
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
int Torus::SOS_z(double* outz, double* outvz,
		 const double RSOS, const int Nthz)
{
	PSPD Jtroot = PSPD(J(0),J(1),0.,0.);
	Jtr[0]=J[0]; Jtr[1]=J[1];
  RforSOS = RSOS;
  int ntry=100,nout=0;
  double tmpR1,tmpR2,thmin,thmax;
  const double allowed = -1.e-8;
  // do this slowly but safely-ish
  for(int ithz=0; ithz<Nthz; ithz++) {
    Jtr[3]=Jtroot[3] = double(ithz) * TPi / double(Nthz);
    // OK first, can I do this the easy way?
    Jtroot[2] = thmin = 0.; tmpR1 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    Jtroot[2] = thmax = Pi; tmpR2 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    for(int i=0;i!=ntry && (tmpR1-RSOS)*(tmpR2-RSOS)>0;i++) {
      if(i) { thmin = thmax; tmpR1 = tmpR2; }
      Jtroot[2] = thmax = double(i+1)*TPi/double(ntry);
      tmpR2 = (Jtroot>>GF>>(*TM)>>(*PT))(0);
    }
    if((tmpR1-RSOS)*(tmpR2-RSOS)<0) {
      Jtroot[2] = rtsafe(this,&Torus::SOS_z_root, thmin, thmax, -allowed);
      PSPD outQP=(Jtroot>>GF>>(*TM)>>(*PT));
      outz[nout]=outQP[1]; outvz[nout]=outQP[3]; nout++;
    }
  }
  return nout;
}
void Torus::RpSOS(ostream& sfile,const int np){//work round with true angle
	double Delta[3],tr;
	for(int i=0;i<np;i++){
		Angles A; A[1]=0; A[2]=0;
		A[0]=i/(double)(np-1)*(TPi-1.e-6);
		GCY gcy=FullMap(A);
		while(gcy[1]<-0.01){
			A[1]+=.001; gcy=FullMap(A);
		}
		while(gcy[1]>0.01){
			A[1]-=.001; gcy=FullMap(A);
		}
		if(gcy[2]>Pi) gcy[2]-=TPi;
		double Ttop,Tbot,Atop,Abot,dA=.1,eps=.001;
		if(fabs(gcy[2])<eps){
			sfile << std::setw(12) << gcy <<'\n'; continue;
		}
		if(gcy[2]<0){
			Abot=A[2]; Tbot=gcy[2];
			while(gcy[2]<0){
				A[2]+=dA;
				gcy=FullMap(A); if(gcy[2]>Pi) gcy[2]-=TPi;
				if(gcy[2]<0){
					Abot=A[2]; Tbot=gcy[2];
				}
			}
			Atop=A[2]; Ttop=gcy[2];
		}else{
			Atop=A[2]; Ttop=gcy[2];
			while(gcy[2]>0){
				A[2]-=dA;
				gcy=FullMap(A); if(gcy[2]>Pi) gcy[2]-=TPi;
				if(gcy[2]>0){
					Atop=A[2]; Ttop=gcy[2];
				}
			}
			Abot=A[2]; Tbot=gcy[2];
		}
		while(1){
			double f=Tbot/(Tbot-Ttop);
			dA=(f*Atop+(1-f)*Abot)-A[2];
			A[2]+=dA;
			gcy=FullMap(A); if(gcy[2]>Pi) gcy[2]-=TPi;
			if(fabs(gcy[2])<eps) break;
			if(gcy[2]<0){
				Abot=A[2]; Tbot=gcy[2];
			}else{
				Atop=A[2]; Ttop=gcy[2];
			}
		}
		sfile << std::setw(12) << gcy <<'\n';
	}
}

/*
////////////////////////////////////////////////////////////////////////////////
// Called by resSOS_p to compute phi and dphi/dtheta_p subject to
// resN_r*theta_r+resN_p*theta_p=t1p
////////////////////////////////////////////////////////////////////////////////
void Torus::resSOSroot_p(const double t2, double& z, double& dz) const
{
	double          dQPdqp[4][4], dqdt[2][2], dqpdj[4][2], djdt[2][2];
	PSPD   JT,jt, QP;
	int    i,j;
	double t1=t2+t1p,dqdt2[2],dz2[2];
	//t1,t2 are the true thetas
//	Jtroot[2]=t1; Jtroot[3] = t2;
	PSPD Jtroot=PSPD(Jtr[0],Jtr[1],t1,t2);
	JT	  = AM.Forward(Jtroot);//get toy angles
	jt        = GF.ForwardWithDerivs(JT, djdt);//get toy actions
	QP        = TM->ForwardWithDerivs(jt, dqdt) >> (*PT);
	PT->Derivatives(dQPdqp);
	TM->Derivatives(dqpdj);
	z  = QP(1);
	for(int k=0;k<2;k++){
		dz2[k] = 0.;
		for(i=0; i<2; i++) {
			dqdt2[k] = dqdt[i][k];
			for(j=0; j<2; j++)
				dqdt2[k] += dqpdj[i][j] * djdt[j][k];
			dz2[k] += dqdt2[k] * dQPdqp[1][i];
		}
	}
	dz=dz2[1]+dz2[0];//for N=(1,-1,0)
}
////////////////////////////////////////////////////////////////////////////////
//Finds where the torus cuts phi=0 subject to the restricton
//resN_r*theta_r+resN_p*theta_p=t1pp
////////////////////////////////////////////////////////////////////////////////
void Torus::resSOS_p(ostream& to,const double t1pp)
//for given _true_ resN_r*theta_r + resN_p*theta_p find _true_ theta_r
//and theta_p that brings orbit to phi=0, z=0
{
	t1p=t1pp;
	PSPD Jtroot = PSPD(J(0),J(1),0.,0.);
	Jtr[0]=J[0]; Jtr[1]=J[1];
	const double allowed = -1.e-8;
	double tempz1,tempz2,tempdz;
	resSOSroot(-Pih,tempz1,tempdz); resSOSroot(Pih,tempz2,tempdz);
	if(tempz1*tempz2<0.){
		Jtr[3]=Jtroot[3] = rtsafe(this,&Torus::resSOSroot, -Pih, -(-Pih), -allowed);
		Jtr[2]=Jtroot[2]=t1p+Jtroot[3];
		to << std::setw(12) << (Jtroot>>AM>>GF>>(*TM)>>(*PT)) <<"   "<<Jtroot(3)<<'\n';
	}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Torus::SetToyPot(Potential * Phi, const Actions Jin, const double R0) {

  IsoPar IP = 0.;
  double Rc = Phi->RfromLc(Jin(1)+fabs(Jin(2)));
  IP[1] = (R0)? sqrt(R0) : sqrt(Rc);
  double dPdR,dPdz, atmp, btmp;
  (*Phi)(Rc,0.,dPdR,dPdz);
  btmp = (R0)? R0 : Rc;
  atmp = (R0)? sqrt(R0*R0+Rc*Rc) : Rc*sqrt(2.);
  IP[0] = sqrt(dPdR*atmp*powf(btmp+atmp,2.)/Rc);
  IP[2] = fabs(Jin(2));
  GenPar SN;
  SN.MakeGeneric();
  GenPar s1(SN);
  s1=0.;
  SetMaps(IP,SN,AngPar(s1,s1,s1));
  //SetActions(Jin);
  J = Jin;
  E = 0.;
  Om = 0.;
  dc = 0.;
}*/
inline double Torus::GM(double r,double F,double bi){//returns GM of isochrone that gives dP/dr=F at r
	double s=1+sqrt(1+pow(r/bi,2));
	return s*pow(s*(s-1)*bi,2)*F/(s-2);
}
void Torus::SetToyPot(Potential * Phi, const Actions Jin, const double R0) {
	IsoPar IP = 0.;
	double Rc = Phi->RfromLc(Jin(1)+fabs(Jin(2)));
	double r1=.75*Rc,r2=1.25*Rc,F1,F2,dPdz;
	(*Phi)(r1,0,F1,dPdz); (*Phi)(r2,0,F2,dPdz);
	double bl=.1,bh=10,eps=1.e-4;
	double dl=GM(r2,F2,bl)-GM(r1,F1,bl),dh=GM(r2,F2,bh)-GM(r1,F1,bh);
	while(dl<0){
		bl*=.5; dl=GM(r2,F2,bl)-GM(r1,F1,bl);
	}
	while(dh>0){
		bh*=1.5; dh=GM(r2,F2,bh)-GM(r1,F1,bh);
	}//now we have bracketed b where masses agree
	while(bh-bl>eps){
		double b=sqrt(bl*bh), d=GM(r2,F2,b)-GM(r1,F1,b);
		if(d<0) bh=b; else bl=b;
	}
	IP[1]=sqrt(.5*(bl+bh));
	IP[0]=sqrt(GM(r2,F2,pow(IP[1],2)));
	IP[2] = fabs(Jin(2));
	GenPar SN;
	SN.MakeGeneric();
	GenPar s1(SN);
	s1=0.;
	SetMaps(IP,SN,AngPar(s1,s1,s1));
	//printf("SetToyPot: %f %f %f %f\n",IP[0],IP[1],IP[2],IP[3]);
  //SetActions(Jin);
	J = Jin;
	E = 0.;
	Om = 0.;
	dc = 0.;
}

void Torus::AutoPTTorus(Potential *Phi, const Actions Jin, const IsoPar IP) {
  GenPar SN;
  SN.MakeGeneric();
  GenPar s1(SN);
  s1=0.;
  J = Jin;
  SetMaps(IP,SN,AngPar(s1,s1,s1));
  SetPP(Phi,Jin);
}

void Torus::FindLimits() { // Approximate.
  Angles Ang = 0.;
  Ang[1] = Pih;
  Rmin = (MapfromToy(Ang))(0);
  Ang[0] = Pi;
  zmax = (MapfromToy(Ang))(1);
  Ang[0] = Pi; Ang[1] =0.;
  Rmax = (MapfromToy(Ang))(0);
}
Torus InterpTorus(Torus *Tgrid,double *Jgrid,int np,double Jr){
	int top,bot;
	if(Jgrid[0]>Jgrid[np-1]){
		top=0; bot=np-1;
	}else{
		top=np-1; bot=0;
	}//now top should point to largest Jr
//	if(Jr>Jgrid[top] || Jr<Jgrid[bot]) printf("off grid in InterpTorus: %f %f %f\n",
//		Jr,Jgrid[bot],Jgrid[bot]);
	Torus T;
	if((Jgrid[top]-Jr)*(Jr-Jgrid[bot])<0.){//Jr lies out of grid
		if(Jr<Jgrid[bot]) T= Tgrid[bot];
		else T=Tgrid[top];
	}else{
		while(abs(top-bot)>1){
			int n=(top+bot)/2;
			if((Jgrid[top]-Jr)*(Jr-Jgrid[n])>=0) bot=n;
			else top=n;
		}
		double f=(Jgrid[top]-Jr)/(Jgrid[top]-Jgrid[bot]);//distance from top
		//printf("InterpTorus: %f ",f);
		if(f){
			T=Tgrid[bot];
			if(1-f) {T*=f; T+=Tgrid[top]*(1-f);
			}
		} else T=Tgrid[top];
	}
	return T;
}
Torus InterpTorus2(Torus *Tgrid,double *Jgrid,int np,double Jr){
	//as above but grid defined only on even-numbered slots
	int top,bot;
	if(np%2==0) printf("InterpTorus2: np must be odd\n");
	if(Jgrid[0]>Jgrid[np-1]){
		top=0; bot=(np-1)/2;
	}else{
		top=(np-1)/2; bot=0;
	}//now top should point to largest Jr
	if(Jr>Jgrid[2*top] || Jr<Jgrid[2*bot]) printf("off grid in InterpTorus: %f %f %f\n",
		Jr,Jgrid[2*bot],Jgrid[2*bot]);
	Torus T;
	if((Jgrid[2*top]-Jr)*(Jr-Jgrid[2*bot])<0.){
		if(Jr<Jgrid[2*bot]) T= Tgrid[2*bot];
		else T=Tgrid[2*top];
	}else{
		while(abs(top-bot)>1){
			int n=(top+bot)/2;
			if((Jgrid[2*top]-Jr)*(Jr-Jgrid[2*n])>=0) bot=n;
			else top=n;
		}
		double f=(Jgrid[2*top]-Jr)/(Jgrid[2*top]-Jgrid[2*bot]);//distance from top
		if(f){
			T=Tgrid[2*bot];
			if(1-f) {T*=f; T+=Tgrid[2*top]*(1-f);
			}
		} else T=Tgrid[2*top];
	}
	return T;
}
Torus InterpTorus(Torus **Tgrid,Actions Jbar,Actions dJ,Actions J){
	//Jbar centre of rectangle, dJ vector across diagonal
	Torus T; T=Tgrid[0][0];
	double fi=(Jbar[0]+.5*dJ[0]-J[0])/dJ[0];//distance from top = weight of bottom
	//printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
	for(int i=0;i<2;i++) {
		if(i==1) fi=1-fi;
		double fk=(Jbar[2]+.5*dJ[2]-J[2])/dJ[2];
		for(int k=0;k<2;k++) {
			if(k==1) fk=1-fk;
			if(fi>1 || fk>1){
				printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",
				       Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
				printf("IT: %f %f\n",fi,fk); exit(0);}
			if(i+k==0) T*=fi*fk;
			else T+=Tgrid[i][k]*(fi*fk);
		}
	}
	return T;
}
Torus InterpTorus_n2(Torus **Tgrid,int nr,Actions Jbar,Actions dJ,Actions J){
	//nr-1 rectangles. Jbar centre of first rectangle, dJ vector across
	//its diagonal
	double j0=sqrt(MAX(.00001,J[0]));
	int n=(j0-(Jbar[0]-.5*dJ[0]))/dJ[0];//# of our rectangle
	if(n<0 || n>nr-2){
		printf("off grid in InterpTorus_n2: cell # n=%d\n",n);
		printf("%f %f\n",J[0],J[2]);
		n=MAX(0,MIN(nr-2,n));
		exit(0);
	}
	Torus T; T=Tgrid[n][0];
	double fi=(Jbar[0]+(n+.5)*dJ[0]-j0)/dJ[0];//distance from top = weight of bottom
	//printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
	for(int i=0;i<2;i++) {
		if(i==1) fi=1-fi;
		if(fi>1.000001) fi=1;
		double fk=(Jbar[2]+.5*dJ[2]-J[2])/dJ[2];//distance from top
		for(int k=0;k<2;k++) {
			if(k==1) fk=1-fk;
			if(fk>1.000001) fk=1;
			if(i+k==0) T*=fi*fk;
			else if(fi*fk>0) T+=Tgrid[n+i][k]*(fi*fk);
		}
	}
	return T;
}
/*Torus InterpTorus_n2(Torus **Tgrid,int nr,Actions Jbar,Actions dJ,Actions J){
	//nr-1 rectangles. Jbar centre of first rectangle, dJ vector across
	//its diagonal
	int n=(J[0]-(Jbar[0]-.5*dJ[0]))/dJ[0];//# of our rectangle
	if(n<0 || n>nr-1){
		printf("off grid in InterpTorus_n2: cell # n=%d\n",n);
		exit(0);
	}
	Torus T; T=Tgrid[n][0];
	double fi=(Jbar[0]+(n+.5)*dJ[0]-J[0])/dJ[0];//distance from top = weight of bottom
	//printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
	for(int i=0;i<2;i++) {
		if(i==1) fi=1-fi;
		double fk=(Jbar[2]+.5*dJ[2]-J[2])/dJ[2];
		for(int k=0;k<2;k++) {
			if(k==1) fk=1-fk;
			if(fi>1.000001 || fk>1.000001){
				printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",
				       Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
				printf("InterpTorus_n2 error: %f %f\n",fi,fk);// exit(0);
			}
			if(i+k==0) T*=fi*fk;
			else T+=Tgrid[n+i][k]*(fi*fk);
		}
	}
	return T;
}*/
Torus InterpTorus_2n(Torus **Tgrid,int np,Actions Jbar,Actions dJ,Actions J){
	//np rectangles. Jbar centre of first rectangle, dJ vector across
	//its diagonal
	int n=(J[2]-(Jbar[2]-.5*dJ[2]))/dJ[2];//# of our rectangle
	if(n<0 || n>np-2){
		printf("off grid in InterpTorus_2n cell # n=%d\n",n);
		printf("%f (%f %f)\n",J[2],Jbar[2]-.5*dJ[2],Jbar[2]+((double)np-1.5)*dJ[2]);
		exit(0);
	}
	Torus T; T=Tgrid[0][n];
	double fi=(Jbar[2]+(n+.5)*dJ[2]-J[2])/dJ[2];//distance from right = weight of left
	//printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
	for(int i=0;i<2;i++) {
		if(i==1) fi=1-fi;
		double fk=(Jbar[0]+.5*dJ[0]-J[0])/dJ[0];
		for(int k=0;k<2;k++) {
			if(k==1) fk=1-fk;
			if(fi>1.001 || fk>1.001){
				printf("IT error: %d %f %f\n",n,fi,fk);
				printf("(%f %f %f) (%f %f %f) (%f %f %f)\n",
				       Jbar[0],Jbar[1],Jbar[2],dJ[0],dJ[1],dJ[2],J[0],J[1],J[2]);
				//exit(0);
			}
			if(i+k==0) T*=fi*fk;
			else T+=Tgrid[k][n+i]*(fi*fk);
		}
	}
	return T;
}
Torus InterpTorus(Torus ***Tgrid,Actions Jbar,Actions dJ,Actions J){
	//Jbar centre of cube, dJ the cube's diagonal vector
	Torus T; T=Tgrid[0][0][0];
	double fi=(Jbar[0]+.5*dJ[0]-J[0])/dJ[0];//distance from top = weight of bottom
	for(int i=0;i<2;i++) {
		if(i==1) fi=1-fi;
		double fj=(Jbar[1]+.5*dJ[1]-J[1])/dJ[1];
		for(int j=0;j<2;j++) {
			if(j==1) fj=1-fj;
			double fk=(Jbar[2]+.5*dJ[2]-J[2])/dJ[2];
			for(int k=0;k<2;k++) {
				if(k==1) fk=1-fk;
				if(i+j+k==0) T*=fi*fj*fk;
				else T+=Tgrid[i][j][k]*(fi*fj*fk);
			}
		}
	}
	return T;
}
