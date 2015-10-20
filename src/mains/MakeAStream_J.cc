
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
using std::cout;

#include "PJM_cline.h"
#include "falPot.h"
#include "Torus.h"
#include "Random.h"
#include "FreeMemory.h"
#include "TNT_JAMA/jama_eig.h"
#include "WDMath.h"


// Torus InterpTorus(Torus ***Tgrid,Actions Jbar,Actions dJ,Actions J){
//     //Jbar centre of cube, dJ side length of cube
//     Torus T; T=Tgrid[0][0][0];
//     for(int i=0;i<2;i++) {
//         double fi=(J[0]-(Jbar[0]-.5*dJ[0]))/dJ[0];
//         for(int j=0;j<2;j++) {
//             double fj=(J[1]-(Jbar[1]-.5*dJ[1]))/dJ[1];
//             for(int k=0;k<2;k++) {
//                 double fk=(J[2]-(Jbar[2]-.5*dJ[2]))/dJ[2];
//                 if(i+j+k==0) T*=fi*fj*fk;
//                 else T+=Tgrid[i][j][k]*(fi*fj*fk);
//             }
//         }
//     }
//     return T;
// }



Vector<double,3> find_eigen(Matrix <double,3,3> in, Matrix <double,3,3> &out) {
  // Find eigenvalues and eigenvectors using JAMA
  Vector<double,3> eigval;
  Array2D <double> covar_eig(3,3), eigvec_jama,  eigval_jama;
  
  for(int j=0;j!=3;j++)
    for(int k=0;k!=3;k++)
      covar_eig[j][k] = in[j][k];
  
  JAMA::Eigenvalue<double> solver(covar_eig);
  solver.getV(eigvec_jama);
  solver.getD(eigval_jama);
  
  for(int j=0;j!=3;j++) {
    eigval[j]    = eigval_jama[j][j];
    for(int k=0;k!=3;k++) {
      out[j][k] = eigvec_jama[j][k];
    }
  }
  
  return eigval;
     
}
Matrix <double,3,3> Find_Hessian(Potential *Phi, Actions J, double dJbyJ) {
  // to find Hessian, I need 6 points on faces of cube around the point J
  // Spacing should be of a reasonable size, without getting too close to 0 
  // at any of the edges
  Matrix <double,3,3> dOmdJ;
  Actions dJ, Om1, dOm, J_t;
  Torus T;

  for(int i=0;i!=3;i++) {
    J_t = J;
    dJ[i] = std::min(0.1,J[i]*0.25);
    J[i] -= dJ[i];
    T.AutoFit(J,Phi,dJbyJ);
    Om1 = T.omega();
    J[i] += 2*dJ[i];
    T.AutoFit(J,Phi,dJbyJ);
    dOm = T.omega()-Om1;
    for(int j=0;j!=3;j++) dOmdJ[j][i] = 0.5*dOm[j]/dJ[i]; 
  }

  {
    // The Hessian should be symmetric. Force it to be.
    double tmp = (dOmdJ[0][1]+dOmdJ[1][0])*0.5;
    dOmdJ[1][0] = dOmdJ[0][1] = tmp;
    tmp = (dOmdJ[0][2]+dOmdJ[2][0])*0.5;
    dOmdJ[2][0] = dOmdJ[0][2] = tmp;
    tmp = (dOmdJ[1][2]+dOmdJ[2][1])*0.5;
    dOmdJ[2][1] = dOmdJ[1][2] = tmp;
  }


  return dOmdJ;
}



template<class T, int N>
bool safe_GaussInvert(Matrix<T,N,N>& A)
{
  register int    i, icol=0, irow=0, j, k, l, ll;
  register T      big, dum, pivinv, One=T(1), Zero=T(0);
  Vector<int,N>   ipiv=0,indxr,indxc;
  for(i=0; i<N; i++) {
    big = Zero;
    for(j=0; j<N; j++)
      if(ipiv(j)!=1)
	for(k=0; k<N; k++) {
	  if(ipiv(k)==0) {
	    if(fabs(A(j,k)) >= big) {
	      big  = fabs(A(j,k));
	      irow = j;
	      icol = k;
	    }
	  } else if(ipiv(k)>1) return false;
	}
    ++(ipiv[icol]);
    if(irow != icol) {
      for(l=0; l<N; l++) std::swap(A[irow][l], A[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(A(icol,icol)==Zero) return false;
    pivinv = One/A(icol,icol);
    A[icol][icol] = One;
    for(l=0; l<N; l++) A[icol][l] *= pivinv;
    for(ll=0; ll<N; ll++)
      if(ll != icol) {
	dum         = A(ll,icol);
	A[ll][icol] = Zero;
	for(l=0; l<N; l++) A[ll][l] -= A(icol,l) * dum;
      }
  }
  for(l=N-1; l>=0; l--)
    if (indxr(l) != indxc(l) )
      for(k=0; k<N; k++) std::swap(A[k][indxr(l)], A[k][indxc(l)]);

  return true;
}

 

int main(int argc,char *argv[])
{
  ifstream from;
  Actions J,Jp;
  Torus T,T2;

  bool show_orbit = true;
  
  int nstars = 10000; // to do - input parameter

  
  if(argc<3) {
    cerr << argv[0] << ": A programme to construct a stream model "
	 << "by interpolation between tori\n"
	 << "Input: potential outputfile_root\n";
    return 0;
  }

  
  Angles A0 = 2.;
  double tmax = 12000.; 
  
  Random3 R3(7),R3b(78);
  Gaussian Gau(&R3,&R3b);

  my_open(from,argv[1]);
  GalaxyPotential Phi(from);
  from.close();

  string Output = string(argv[2]);
  
  // Progenitor J, sig_x, sig_v, mu
  // Future: leave as inputs
  Jp[0] = 0.5;
  Jp[1] = 0.5;
  Jp[2] = 3.5;
  
  double sig_x = 0.05,      // spread in position
    sig_v = 0.5*Units::kms, // spread in velocity
    mu =5.;                 // Used roughly as per Bovy 2014
  

  // Things needed to pick a principle direction in action 
  
  Matrix<double,3,3> Hess = Find_Hessian(&Phi,Jp,0.001), // dOm/dJ
    InvHess, covarience_J =0., HessEigVecs,
    covarience_om=0., eigenvecs;

  Vector<double,3>   eigenvals, HessEigVals, testev,
    testbasis[3], OmVectors[3], MainOmVec, MainJVec, NormMainJVec;

  T.AutoFit(Jp,&Phi,0.001);  

  

  
  // sigJr = sigv(rapo-rperi)/pi
  // sigJz = 2*sigv*Zmax/pi
  // sigJp = sigv*rperi

  covarience_J[0][0] =  sig_v*(T.maxR()-T.minR())/Pi; 
  covarience_J[1][1] =  2.*sig_v*(T.maxz())/Pi; 
  covarience_J[2][2] = sig_v*T.minR();

  Angles delA;
  // Similar 
  for(int i=0;i!=3;i++) delA[i] = sig_x*sig_v/covarience_J[i][i];
  
  // covarience -> sigma^2
  covarience_J = covarience_J*covarience_J;

  // Propogate covariece in J to one in omega with Hessian
  covarience_om = Hess*(covarience_J*!Hess);

  // Find eigenvalues of covarience matrix 
  eigenvals = find_eigen(covarience_om,eigenvecs);


  // Find principle direction in J space
  {
    // eigenvalues are, unhelpfully, not given in order of absolute value
    // Find the largest
    int whichFreqVec = (fabs(eigenvals[0])>fabs(eigenvals[1]))? 0 : 1;
    if( fabs(eigenvals[2])>fabs(eigenvals[0]) &&
        fabs(eigenvals[2])>fabs(eigenvals[1]) ) whichFreqVec = 2;

    // Principle direction in frequency space
    for(int j=0;j!=3;j++) MainOmVec[j] = eigenvecs[j][whichFreqVec];

    // Propagate to action space...

    // First find dJ/dOm
    InvHess = Hess;
    if(!safe_GaussInvert(InvHess)) {
      cerr << "Unable to invert Inverse Hessian. What do I do now?\n";
    }

    // Then convert vector in action space
    MainJVec = InvHess*MainOmVec;
    // normalise
    NormMainJVec = MainJVec/sqrt(MainJVec.norm());
    // cerr << NormMainJVec << '\n';
  }


  // Pick two perpendicular directions, and find related dispersion
  
  Vector<double,3> JVecPerp1=0.,JVecPerp2=0.;
  double VariencePar=0.,VariencePerp=0., DispPar,DispPerp;

  {
    // A (simple) perpendicular vector to the main axis
    JVecPerp1[0] = NormMainJVec[1];
    JVecPerp1[1] = -NormMainJVec[0];
    double tmp = sqrt(JVecPerp1.norm());
    JVecPerp1 = JVecPerp1/tmp;
    
    // A third axis
    JVecPerp2 = NormMainJVec^JVecPerp1;

    
    Vector<double,3> tmpV = covarience_J * NormMainJVec;
    
    // varience in main direction
    for(int i=0;i!=3;i++) VariencePar += tmpV[i]*NormMainJVec[i];
    
    // variance in perpendicular directions
    VariencePerp = sqrt(covarience_J[0][0]*covarience_J[1][1]*
			covarience_J[2][2]/VariencePar);
    
    DispPar  = sqrt(VariencePar);
    DispPerp = sqrt(VariencePerp);
    //cerr << DispPar << ' ' << DispPerp << '\n';
  }



  
  // Now I need to check that my grid will be large enough.
  // What's the largest likely value of JR,Jz,Jp.

  double dJRmax = DispPar*(mu+3.)*fabs(NormMainJVec[0]) +
    DispPerp*5.*( fabs(JVecPerp1[0]) + fabs(JVecPerp2[0]) );
  double dJzmax = DispPar*(mu+3.)*fabs(NormMainJVec[1]) +
    DispPerp*5.*( fabs(JVecPerp1[1]) + fabs(JVecPerp2[1]) );
  double dJpmax = DispPar*(mu+3.)*fabs(NormMainJVec[2]) +
    DispPerp*5.*( fabs(JVecPerp1[2]) + fabs(JVecPerp2[2]) );

  
  // for now, 2 by 2 grid only
  // If it expands to a larger grid, then I *suspect* that the best spacing
  // of points is linear in Jphi, logarithmic in Jr, Jz

  int nJR = 2;
  int nJz = 2;
  int nJp = 2;

  double JRmin = (dJRmax<0.8*Jp[0])?  Jp[0] - dJRmax : 0.2*Jp[0];
  double JRmax =  Jp[0] + dJRmax;
  double JRfac = powf(JRmax/JRmin,1./double(nJR-1));
  double Jzmin = (dJzmax<0.8*Jp[1])? Jp[1] - dJzmax : 0.2*Jp[1];
  double Jzmax =  Jp[1] + dJzmax;
  double Jzfac = powf(Jzmax/Jzmin,1./double(nJz-1));
  double Jpmin = (dJpmax<0.8*Jp[2])?  Jp[2] - dJpmax : 0.2*Jp[2];
  double Jpmax =  Jp[2] + dJpmax;
  double Jpdiff = (Jpmax-Jpmin)/double(nJp-1);
  
  double JRgrid[nJR];
  double Jzgrid[nJz];
  double Jpgrid[nJp]; 

  for(int i=0;i!=nJR;i++) { JRgrid[i] = JRmin * powf(JRfac,i); }
  for(int i=0;i!=nJz;i++) { Jzgrid[i] = Jzmin * powf(Jzfac,i); }
  for(int i=0;i!=nJp;i++) { Jpgrid[i] = Jpmin + i*Jpdiff; }

  // Grids of parameters
  
  vec4 ***TPgrid = PJM::matrix<vec4>(nJR,nJz,nJp);
  GenPar ***SNgrid = PJM::matrix<GenPar>(nJR,nJz,nJp);
  AngPar ***APgrid = PJM::matrix<AngPar>(nJR,nJz,nJp);
  Frequencies ***Freqgrid = PJM::matrix<Frequencies>(nJR,nJz,nJp);

  vec4 TPtest=0.;
  GenPar SNtest=0.;
  AngPar APtest=0.;
  Frequencies Freqtest=0.;
  double fac[3];

  // Fill the grids of parameters  
  for(int i=0;i!=nJR;i++) {
    J[0] = JRgrid[i];
    for(int j=0;j!=nJz;j++) {
      J[1] = Jzgrid[j];
      for(int k=0;k!=nJp;k++) {
  	J[2] = Jpgrid[k];
  	T.SetToyPot(&Phi,J);
  	vec4 TP = T.TP();
  	T.FitWithFixToyPot(J,TP,&Phi,0.001);
  	TPgrid[i][j][k] = T.TP();
  	SNgrid[i][j][k] = T.SN();
  	APgrid[i][j][k] = T.AP();
  	Freqgrid[i][j][k] = T.omega();
      }
    }
  }

    Frequencies Om0;
    for(int ii=0;ii!=2;ii++) { 
      for(int jj=0;jj!=2;jj++) {
	for(int kk=0;kk!=2;kk++) {
	  TPtest += 0.125*TPgrid[ii][jj][kk];
	  SNtest += SNgrid[ii][jj][kk]*0.125;
	  APtest += APgrid[ii][jj][kk]*0.125;
	  Om0 += 0.125*Freqgrid[ii][jj][kk];
	}
      }
    }
    
  //  cerr << Om0 << '\n';
    
  // Output path of the orbit?
  if(show_orbit) {
    
    T.SetActions(Jp);
    T.SetTP(TPtest);
    T.SetSN(SNtest);
    T.SetAP(APtest);
    T.SetFrequencies(Om0);
  
    
    ofstream to_orb;

    
    my_open(to_orb, Output+".orb");
    
    for(int i=-6000;i!=6000;i++) {
      Angles Atmp = A0;
      double dA0 = 0.001*i; 
      Atmp[0] += dA0;
      Atmp[1] += dA0*Om0[1]/Om0[0];
    Atmp[2] += dA0*Om0[2]/Om0[0];
    to_orb << T.Map3D(Atmp) << '\n';
    }

  }

  ofstream to;
  my_open(to, Output+".st");
  
  // Create stars
  for(int i=0;i!=nstars;i++) {
    Angles A;
    Actions dJ;
    J = Jp;
    string type = "bowtie";
    if(type=="nobowtie") {
      // mu + Gaussian in principle direction
      dJ = (R3()>0.5)? (mu+Gau()) * DispPar * NormMainJVec :
	(-mu-Gau()) * DispPar * NormMainJVec;
      // Gaussian in perpendicular directions
      dJ += Gau() * DispPerp * JVecPerp1;
      dJ += Gau() * DispPerp * JVecPerp2;
    } else {
      // bowtie
      double tmp = Gau(), PerpFac = 1.+tmp/3.;
      PerpFac = 1+tanh(tmp/3.);
      dJ = (R3()>0.5)? (mu+tmp) * DispPar * NormMainJVec :
	(-mu-tmp) * DispPar * NormMainJVec;
      dJ += Gau() * DispPerp * JVecPerp1 * PerpFac;
      dJ += Gau() * DispPerp * JVecPerp2 * PerpFac;
      
    }
    J += dJ;
    for(int j=0;j!=3;j++) A[j] = A0[j] + delA[j]*Gau();
    double t = tmax * R3();
    
    double tt,uu,vv;
    int JRlo=0,Jzlo=0,Jplo=0;

    
    TPtest=0.;
    SNtest=0.;
    APtest=0.;
    Freqtest=0.;
    tt = (J[0] - JRgrid[JRlo])/(JRgrid[JRlo+1] - JRgrid[JRlo]); 
    uu = (J[1] - Jzgrid[Jzlo])/(Jzgrid[Jzlo+1] - Jzgrid[Jzlo]); 
    vv = (J[2] - Jpgrid[Jplo])/(Jpgrid[Jplo+1] - Jpgrid[Jplo]); 
    for(int ii=0;ii!=2;ii++) {
      fac[0] = (ii)? tt : 1-tt;
      for(int jj=0;jj!=2;jj++) {
	fac[1] = (jj)? uu : 1-uu;
	for(int kk=0;kk!=2;kk++) {
	  fac[2] = (kk)? vv : 1-vv;
	  TPtest += fac[0]*fac[1]*fac[2]*TPgrid[JRlo+ii][Jzlo+jj][Jplo+kk];
	  SNtest += SNgrid[JRlo+ii][Jzlo+jj][Jplo+kk]*fac[0]*fac[1]*fac[2];
	  APtest += APgrid[JRlo+ii][Jzlo+jj][Jplo+kk]*fac[0]*fac[1]*fac[2];
	  Freqtest += Freqgrid[JRlo+ii][Jzlo+jj][Jplo+kk]*fac[0]*fac[1]*fac[2];
	}
      }
    }
    T.SetActions(J);
    T.SetTP(TPtest);
    T.SetSN(SNtest);
    T.SetAP(APtest);
    T.SetFrequencies(Freqtest);
    
      

      
    //T.AutoFit(J,&Phi);
    Frequencies Om = T.omega();
    
    // Add the drift away in angle
    for(int j=0;j!=3;j++) {
      A[j] += (Om[j]-Om0[j])*t; 
    }
      
    to << J << ' ' << A << ' ' << T.FullMap(A) << '\n'<< std::flush; 
  }
}
