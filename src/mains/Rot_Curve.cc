/**
\file Rot_Curve.cc
\brief Should give some file that demonstrates potentials, but not this.

*                                                                              *
* C++ code written by Paul McMillan, 2013                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
using std::cout;

#include "falPot.h"
#include "PJM_cline.h"


int main(int argc,char *argv[])
{
  if(argc<5) {
    cerr << "Gives vc(R) for a set of values R from Rmin to Rmax separated by dR\n"
	 << "Input: potential Rmin Rmax dR\n";
  }
  ifstream from;

  my_open(from,argv[1]);
  GalaxyPotential Phi(from);
  from.close();

  double Rmin = atof(argv[2]),
    Rmax = atof(argv[3]),
    dR = atof(argv[4]);
  
  if(Rmin<=0. || (Rmax-Rmin)/dR<0. || dR==0) {
    cerr << "Bad input parameters (e.g. Rmin=0  which is not OK)\n";
    exit(1);
  }

  for(double R=Rmin;R<=Rmax;R+=dR) {
    cout << R << ' ' << sqrt(Phi.vcsquare(R))*Units::kms_i << '\n';
  }

}
