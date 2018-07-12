/***************************************************************************//**
\file AUTOTORI.cc
\brief Automatically generates file with tori of many (given) actions

*                                                                              *
* AUTOTORI.cc                                                                  *
*                                                                              *
* C++ code written by Paul McMillan, 2014                                      *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "PJMebf.h"
#include "Torus.h"
#include "falPot.h"
//#include "LogPot.h"
//#include "MiyamotoNagaiPot.h"
#include "PJM_cline.h"

using std::cerr;



//int GridAndGrid(Potential *Phi, double ***inact, double *****outact 

int main(int argc,char *argv[])
{
  if(argc<2) {
    cerr << "Input: potential output_toruslist output_values\n";
    if(argc>1)
      if(argv[1][0] == 'h') {
	cerr << "\nCode that computes a grid of tori and outputs whenever " 
	     << "they go through a grid of points in R and z\n";
      }
    exit(0);
  }

  double dJbyJ = 0.003;  
  // parameter that tells the programme how accurately to try to fit the torus 

  ifstream from;
  double Jrmin,Jrmax,dJr,
    Jzmin,Jzmax,dJz,
    Jpmin,Jpmax,dJp;
  int nJr, nJz, nJp;

  double Rmin, Rmax, dR, zmin, zmax, dz;
  int nR, nz;

  // Action grid parameters
  nJr = 11;
  nJz = 11;
  nJp = 11;

  Jrmin = 0.001; Jrmax = 0.2001; dJr = (Jrmax-Jrmin)/double(nJr-1); 
  Jzmin = 0.001; Jzmax = 0.2001; dJz = (Jzmax-Jzmin)/double(nJz-1); 
  Jpmin = 0.1;  Jpmax = 2.1; dJp = (Jpmax-Jpmin)/double(nJp-1); 

  // spacial grid parameters
  nR = 21;
  nz = 21;
  Rmin = 0.1; Rmax = 10.1; dR = (Rmax-Rmin)/double(nR-1);
  zmin = 0.; zmax = 2.0;   dz = (zmax-zmin)/double(nz-1);

  // Potential  
  my_open(from,argv[1]);
  Potential *Phi = new GalaxyPotential(from);
  from.close();

  // Output files
  ebf::WriteString(string(argv[2]),"/log","Torus list","w");
  ofstream to;
  my_open(to,argv[3]);

  Torus T;
  Actions J;
  Position x=0.;
  Velocity v[4];
  int nt = 0;


  // Big loop. Very multi-layered.

  for(int i=0;i!=nJp;i++) {
    J[2] = Jpmin + i*dJp;
    for(int j=0;j!=nJz;j++) {
      J[1] = Jzmin + j*dJz;
      for(int k=0;k!=nJr;k++) {
	J[0] = Jrmin + k*dJr;
	
	// Fit torus
	T.AutoFit(J,Phi,dJbyJ);
	
	// Store torus
	std::stringstream tor; tor << "t" << ++nt;
	T.write_ebf(string(argv[2]),tor.str());
	
	// Check if torus goes through point. If so, with what velocity? 
	for(int l=0;l!=nR;l++) {
	  x[0] = Rmin + l*dR;
	  for(int m=0;m!=nz;m++) {
	    x[1] = zmin + m*dz;
	    
	    // Don't need to check any further if point is clearly outside torus
	    // minR, maxR and maxz give (approx) extent of torus
	    if(x[0]>0.99*T.minR() && x[0] <1.01*T.maxR() 
	       && x[1] <1.01*T.maxz()) {
	      
	      if(T.containsPoint_Vel(x,v[0],v[1],v[2],v[3])) {
		// Output if torus goes through point
		to << J << ' ' << x << ' ';
		for(int ii=0;ii!=4;ii++) to << v[ii] << ' ';
		to << '\n';
	      }

	    }
	  }
	}
      }
    }
  }




}
