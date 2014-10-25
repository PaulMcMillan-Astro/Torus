/*******************************************************************************
*                                                                              *
* PaperGRAPH.cc                                                                *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2007                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "Orb.h"
#include "Torus.h"
#include "falPot.h"
#include "Random.h"
#include "PJMCoords.h"
#include "PJM_cline.h"

using std::cout;
using std::cerr;
using std::cin;

/*
SuperMongo macro to plot output in Rz plane from this programme:

TorRz 1		 data $1.tororb
      		 read {_RT 3 _zT 4}
		 ctype default
		 lweight 4
		 limits _RT _zT box 
		 lweight 5 connect _RT _zT
		 data $1.intorb
      		 read {_RI 3 _zI 4}
		 ctype red ltype 2
		 lweight 2 connect _RI _zI
		 ltype 0
		 ctype default    

 */


int main(int argc,char *argv[])
{
 
  int NT=1000, Ndyn=15,Nlarge=500000,flag;
  time_t second=time(NULL);
  long cp=second;
  float ti[Nlarge],Ri[Nlarge], zi[Nlarge],phii[Nlarge];
  PSPD StPo, oJT, JT;
  Actions J;
  Torus    *T, T0;
  string basename;
  ifstream from;
  ofstream to;
  Potential *Phi;
  Random3 r3(cp);

  T = new Torus;
 
//---------------------------------------------------------------------------
  if( argc !=6 && argc !=7 ){
    cerr << "Outputs orbit found from torus and from orbit integration\n"
	 <<"Input: potential JR Jz Jphi output_file_root (dJ/J)\n";
    exit(1);
  }
//---------------------------------------------------------------------------
//===========================================================================
// 1. Get target potential
  my_open(from, argv[1]);
  Phi = new GalaxyPotential(from);
  from.close();

// input actions
  J[0] = atof(argv[2]); J[1] = atof(argv[3]); J[2] = atof(argv[4]);
  double dJ = (argc==7)? atof(argv[6]) : 0.003;
  bool err = false;
// Use the Torus member functions AutoTorus to take a first (weak) guess -------
  
  flag = T->AutoFit(J,Phi,dJ);
  // std::cout << J << ", flag = " << flag<< "\n" << std::flush;

  Phi->set_Lz(J(2));


// The Torus part -------------------------------------------------------------
    Angles A0;
    PSPT Q;
    Frequencies om=T->omega(); double omratz=1.,omratp=1.,t=0; 
    if(om(0) && om(1)) omratz = om(1)/om(0);
    if(om(0) && om(2)) omratp = om(2)/om(0);
    Angles A =0.;
    A[0] = r3.RandomDouble()*TPi; A[1] = r3.RandomDouble()*TPi;
    //A[0] = Pi/3.; A[1] = 0.;

    A0 = A;
    double dOr=TPi*double(Ndyn)/double(NT),dOl=omratz*dOr,dOp=omratp*dOr;
  
    JT[0] = T->action(0); JT[1] = T->action(1);

    double Ezmin = 1e99, Ezmax = -1e99;

    // Output points from the orbit starting at A0
    for(int i=0; i<NT;) {
      A[0] += dOr; A[1] += dOl; A[2] += dOp;
      for(int j=0;j!=3;j++) if(A(j) > TPi) A[j] -=TPi;
      Q  = T->Map3D(A);
      t += dOr/om(0);

      ti[i] = t;  Ri[i] =  Q(0); zi[i] =  Q(1); phii[i++] =  Q(2);
      // double Ez = 0.5*Q(4)*Q(4) + (*Phi)(Q(0),Q(1)) - (*Phi)(Q(0),0.);
      // if(Ez<Ezmin) Ezmin = Ez; 
      // if(Ez>Ezmax) Ezmax = Ez; 
      // cerr << Q(0) << ' ' << Q(1) << ' ' << Ez << '\n';
    }
    float tT[NT],RT[NT],zT[NT],phiT[NT];
    for(int i=0; i!=NT; i++) {
      tT[i] = ti[i]; RT[i] = Ri[i]; zT[i] = zi[i]; phiT[i] = phii[i];
    }
    //cerr << "Ezmin/max: " << Ezmin << ' ' << Ezmax << '\n'; 
 
// The orbit integration part --------------------------------------------------
    int NO = 0;
    t=0.;
    double dt = 1.e-5;
    Q  = T->Map3D(A0);

    Record3D X(Q,Phi);
    X.set_tolerance(1.e-12);
    X.set_maxstep(0.1);
    // Output points from same orbit found by integration.
    for(;t<tT[NT-1];) {
      X.stepRK_by(dt);
      t += dt;
      if(fabs(X.QP3D(2)-Q(2)) > 0.01*Pi && t<tT[NT-1]) {
	Q = X.QP3D();
	ti[NO] = t;  Ri[NO] =  Q(0); zi[NO] =  Q(1); phii[NO++] =  Q(2);

      } 
    }

    float tO[NO], RO[NO],zO[NO],phiO[NO];
    for(int i=0; i!=NO; i++) {
      tO[i] = ti[i]; RO[i] = Ri[i]; zO[i] = zi[i]; phiO[i] = phii[i];
    } 

    // correct so they end at the same time
    // tO[NO-1] < tT[NT-1] always
    {
      int nstop=NT-1;
      while(tT[nstop-1]>tO[NO-1] && nstop>0) nstop--;
      // interpolate R,z,phi. Don't worry too much about phi
      double frac = (tO[NO-1]-tT[nstop-1])/(tT[nstop]-tT[nstop-1]);
      RT[nstop] = RT[nstop-1] + frac*(RT[nstop]-RT[nstop-1]);
      zT[nstop] = zT[nstop-1] + frac*(zT[nstop]-zT[nstop-1]);
      phiT[nstop] = phiT[nstop-1] + frac*(phiT[nstop]-phiT[nstop-1]);
      tT[nstop] = tO[NO-1];
      NT = nstop+1;
    }

    string TorOut = string(argv[5]) + ".tororb",
      OrbOut = string(argv[5]) + ".intorb";

    my_open(to,TorOut);
    for(int i=0;i!=NT;i++) {
      if(i)
	if(phiT[i]< phiT[i-1]) {
	  // interpolate
	  double diff = 1. - (phiT[i]/(phiT[i]-(phiT[i-1]-TPi)));
	  to << i-1+diff << ' ' <<  tT[i-1]+diff*(tT[i]-tT[i-1]) << ' ' 
	     << RT[i-1]+diff*(RT[i]-RT[i-1]) << ' ' 
	     << zT[i-1]+diff*(zT[i]-zT[i-1]) << " 6.283\n"
	     << i-1+diff << ' ' <<  tT[i-1]+diff*(tT[i]-tT[i-1]) << ' ' 
	     << RT[i-1]+diff*(RT[i]-RT[i-1]) << ' ' 
	     << zT[i-1]+diff*(zT[i]-zT[i-1]) << " 0\n";
	}
      to << i  << ' ' <<  tT[i] << ' ' <<  RT[i]  << ' ' << zT[i]  << ' ' << phiT[i]<<'\n';

    }
    to.close();

    my_open(to,OrbOut);
    for(int i=0;i!=NO;i++) {
      if(i)
	if(phiO[i]< phiO[i-1]) {
	  // interpolate
	  double diff = 1. - (phiO[i]/(phiO[i]-(phiO[i-1]-TPi)));
	  to << tO[i-1]/(dOr/om(0))+diff << ' ' 
	     <<  tO[i-1]+diff*(tO[i]-tO[i-1]) << ' ' 
	     << RO[i-1]+diff*(RO[i]-RO[i-1]) << ' ' 
	     << zO[i-1]+diff*(zO[i]-zO[i-1]) << " 6.283\n"
	     << tO[i-1]/(dOr/om(0))+diff << ' ' 
	     <<  tO[i-1]+diff*(tO[i]-tO[i-1]) << ' ' 
	     << RO[i-1]+diff*(RO[i]-RO[i-1]) << ' ' 
	     << zO[i-1]+diff*(zO[i]-zO[i-1]) << " 0\n";
	    }
      to << tO[i]/(dOr/om(0)) << ' ' <<  tO[i] << ' ' 
	 <<  RO[i]  << ' ' << zO[i]  << ' ' << phiO[i]<<'\n';
    }
    to.close();

}
