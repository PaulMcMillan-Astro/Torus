/***************************************************************************//**
\file Orbit_basic.cc
\brief Should delete this before release

*                                                                              *
* Orbit_basic.cc                                                               *
*                                                                              *
* C++ code written by Paul McMillan, 2013                                      *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <vector>
//#include "PJMsm.h"
#include "Orb.h"
//#include "Torus.h"
#include "falPot.h"
#include "LogPot.h"
//#include "Random.h"
#include "PJMCoords.h"
#include "PJM_cline.h"
using std::cout;
using std::cerr;
using std::cin;


int main(int argc,char *argv[])
{

  if(argc<5) {
    cerr << "Input: potential input output_root output_time "
	 << "(output_time_2... maximum of 5)\n"
	 << argv[0] << " h for help\n";
    if(argc>1) if(argv[1][0] == 'h') {
	cerr << "Input file must have x,v in Galactocentric Cartesian "
	     << "coordinates in the first 6 columns\n"
	     << "Units kpc, km/s. Output times in accending order, in Myr\n";
      }
    exit(0);
  }

  ifstream from;
  Potential *Phi;

  ofstream to0,to1,to2,to3,to4,to5;
  string root = string(argv[3]),out1=root+".1.tab",out2=root+".2.tab",
    out3=root+".3.tab",out4=root+".4.tab",out5=root+".5.tab";
  
  if(string(argv[1]) == "LogPotential_220") {
    Phi = new LogPotential(220.*Units::kms,0.8,0.);
  } else {
    my_open(from,argv[1]);
    Phi = new GalaxyPotential(from);
    from.close();
  }


  my_open(from,argv[2]);
  int n_per_line = entrys_in_line(from);
  bool need_dummy = (n_per_line > 6); 
  double times[5];

  for(int i=0;i<5 && i+4<argc;i++) {
    times[i] = atof(argv[i+4]);
    if(i) if(times[i] <= times[i-1]) { 
	cerr << "times in sequence please\n";
	exit(1);
      }
  }

  my_open(to1,out1);
  if(argc>5) my_open(to2,out2);
  if(argc>6) my_open(to3,out3);
  if(argc>7) my_open(to4,out4);
  if(argc>8) my_open(to5,out5);
  int ntimes = argc-4;

  OmniCoords OC;
  GCA sGCA; 
  GCY sGCY;
  PSPT Q;
  int n = how_many_lines(from);
  string dummy = "";

  // Main loop - for each particle find position after orbit integration
  for(int i=0;i!=n;i++) {
    double t=0., dt =1.e-5;
    std::vector<bool> got(ntimes,false), test_got(ntimes,false);
    from >> sGCA; 
    if(need_dummy) std::getline(from,dummy);

    for(int j=3;j!=6;j++) sGCA[j] *= Units::kms; // units
    sGCY = OC.GCYfromGCA(sGCA);                  // convert to cyl polars
    for(int j=0;j!=6;j++) Q[j] = sGCY[j];        // put into needed form
    Phi->set_Lz(Q(0)*Q(5));
    Record3D X(Q,Phi);                           // integration class
    X.set_tolerance(1.e-12);
    X.set_maxstep(0.1);
    while(t<=times[ntimes-1]) {
      X.stepRK_by(dt);                      // Runge Kutta step (up to maxstep)
      t += dt;
      for(int j=0;j!=ntimes;j++)    
	if(t>times[j] && ! got[j]) {                // check if desired timestep
	  for(int k=0;k!=6;k++) sGCY[k] = X.QP3D(k);
	  sGCA = OC.GCAfromGCY(sGCY);
	  for(int k=3;k!=6;k++) sGCA[k] *= Units::kms_i; // to original form
	  if(!j)   to1 << sGCA << " " << dummy << "\n";
	  if(j==1) to2 << sGCA << " " << dummy << "\n";
	  if(j==2) to3 << sGCA << " " << dummy << "\n";
	  if(j==3) to4 << sGCA << " " << dummy << "\n";
	  if(j==4) to5 << sGCA << " " << dummy << "\n";
	  got[j] = true;
	}
    }
  }
}







