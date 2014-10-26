//
// C++ code by Paul McMillan
//
// Using the Torus machinery to populate a very simple stream model
//

/**
\file SmoothStreamMaker.cc
\brief Make a stream, moderately simple prescription

Creates a stream. Parameters input are potential, progenitor actions & angles,
initial velocity dispersion & size of satellite,
time since first stripping (Myr) and number of stars wanted.

We assume that stars have been stripped at a constant rate since the
first one was stripped, and that the actions and initial angles are
independant of the time of stripping.

Following Eyre & Binney 2012 we take the spread in J to be 

sig_JR = sig_v*(r_ap-r_peri)/pi
sig_Jz = 2*sig_v*Zmax/pi
sig_J  = sig_v*r_peri

And the initial spread in theta is (in terms of an effective radius)

sig_thR   = R_e * pi/(r_ap-r_peri)
sig_thz   = R_e * pi/(2*Zmax)
sig_thphi = R_e / r_peri

Which rather pleasingly means that we have 
del^3 J del^3 th = del^3 x del^3 v 
without any fine tuning.

Output is in code units and galactocentric cylindrical polar coordinates.

 */
#include <iostream>
#include <iomanip>
#include <fstream>



#include "Random.h"
#include "LogPot.h"
#include "falPot.h"
#include "Torus.h"
#include "PJM_utils.h"

using std::cout;
using std::cerr;




int main(int argc,char *argv[])
{
  Random3 R3(6), R3b(72);
  Gaussian Gau(&R3,&R3b);
  ifstream from;
  ofstream to;
  Torus T;

  if(argc<13) {
    cerr << "Creates a stream given potential, progenitor J & th, initial velocity dispersion (code units), progenitor radius, "
	 << "time since first stripping (Myr) and number of stars\n";
    cerr << "Input: Potential J_r J_z Jphi th_r th_z th_phi del_J del_th tmax nstars output_file\n";
    cerr << "e.g. "<<argv[0]<<" pot/PJM11.Tpot 0.2 0.4 3.5 0.5 0.9 2 0.002 0.1 2000 400 tmp.st\n";
    exit(0);
  }
  int nstars;
  Actions J0, delJ, J;
  Angles A0, delA, A;
  Frequencies Om0, Om;
  double tmax, t, rp, ra, zmax, delv, delx;

  // Read in parameters
  my_open(from,argv[1]);
  my_open(to,argv[12]);
  GalaxyPotential Phi(from);
  from.close();

  for(int i=0;i!=3;i++) {
    J0[i] = atof(argv[i+2]);
    A0[i] = atof(argv[i+5]);
  }

  T.AutoFit(J0,&Phi);
  Om0 =  T.omega();
  rp = T.minR();
  ra = T.maxR();
  zmax = T.maxz();

  delv = atof(argv[8]);
  delx = atof(argv[9]);

  delJ[0] = delv*(ra-rp)/Pi;
  delJ[1] = 2*delv*zmax/Pi;
  delJ[2] = delv*rp;

  for(int i=0;i!=3;i++) delA[i] = delx*delv/delJ[i];

  for(int i=0;i!=3;i++) cerr << delJ[i] << '\t';
  cerr << '\n';
  for(int i=0;i!=3;i++) cerr << delA[i] << '\t';
  cerr << '\n';

  tmax = atof(argv[10]);
  nstars = atoi(argv[11]);


  // Create stars
  for(int i=0;i!=nstars;i++) {
    for(int j=0;j!=3;j++) {
      do{
	J[j] = J0[j] + delJ[j]*Gau();
      } while(J[j]<0 && j<2);
      A[j] = A0[j] + delA[j]*Gau();
    }
    t = tmax * R3();

    T.AutoFit(J,&Phi);
    Om = T.omega();
    
    // Add the drift away in angle
    for(int j=0;j!=3;j++) {
      A[j] += (Om[j]-Om0[j])*t; 
    }

    to << T.FullMap(A) << '\n';

  }



}


