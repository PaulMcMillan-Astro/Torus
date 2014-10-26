//
// C++ code by Paul McMillan
//
// Using the Torus machinery to populate a very simple stream model
//

/**
\file SimpleStreamMaker.cc
\brief Make a stream, very simple prescription

Creates a stream. Parameters input are potential, progenitor actions & angles,
spread in J (isotropic), initial spread in theta (also isotropic),
time since first stripping (Myr) and number of stars wanted.

We assume that stars have been stripped at a constant rate since the
first one was stripped, and that the actions and initial angles are
independant of the time of stripping.

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
    cerr << "Creates a stream given potential, progenitor J & th, initial spread in J & th, "
	 << "time since first stripping (Myr) and number of stars\n";
    cerr << "Input: Potential J_r J_z Jphi th_r th_z th_phi del_J del_th tmax nstars output_file\n";
    cerr << "e.g. "<<argv[0]<<" pot/PJM11.Tpot 0.2 0.4 3.5 0.5 0.9 2 0.01 0.005 2000 400 tmp.st\n";
    exit(0);
  }
  int nstars;
  Actions J0, J;
  Angles A0, A;
  Frequencies Om0;
  double delJ, delA, tmax, t;

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

  delJ = atof(argv[8]);
  delA = atof(argv[9]);
  tmax = atof(argv[10]);
  nstars = atoi(argv[11]);


  // Create stars
  for(int i=0;i!=nstars;i++) {
    for(int j=0;j!=3;j++) {
      do{
	J[j] = J0[j] + delJ*Gau();
      } while(J[j]<0 && j<2);
      A[j] = A0[j] + delA*Gau();
    }
    t = tmax * R3();

    T.AutoFit(J,&Phi);
    Frequencies Om = T.omega();
    
    // Add the drift away in angle
    for(int j=0;j!=3;j++) {
      A[j] += (Om[j]-Om0[j])*t; 
    }

    to << T.FullMap(A) << '\n';

  }



}


