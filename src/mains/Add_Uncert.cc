
/*! \file Add_Uncert.cc
  \brief Adds "observational uncertainty" in distance modulus, vr & proper motion to input data.
 *
 *  No further explanation needed.
 */




#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;


#include "PJMCoords.h"
#include "Random.h"
#include "PJM_cline.h"

int main(int argc,char *argv[])
{
  ifstream  from;
  ofstream  to;
  time_t    cpu=time(NULL);
  bool      d_dm=false;
  int seed=0;
  if(argc>6) seed = 12734*atoi(argv[6]);
  Random3   R3(7*int(cpu)+seed),R3b(123*int(cpu)+seed);
  Gaussian  Gau(&R3,&R3b);
 
  if(argc < 6) {
    cerr <<  argv[0] << ": given input in HGP or HEQ, supplies output with quoted errors in distance modulus, vr, proper motion\n";
    cerr << "Usage: "<< argv[0] << " input output d_dm d_vr d_pm\n";
    exit(1);
  }

  double    mag, deldm, delvr, delpm;
  HGP sHGP;
  Vector <double,3> unc;
  my_open(from,argv[1]);
  my_open(to,argv[2]);
  deldm = atof(argv[3]); 
  delvr  = atof(argv[4]);  
  delpm  = atof(argv[5]); 

  int nline = how_many_lines(from);
 
  for(int i=0;i!=nline;i++) {
    from >> sHGP;
    if(deldm) {
      double tmp = 5.*log10(sHGP(0)) + 10.; // distance modulus
      tmp += deldm*Gau();
      sHGP[0] = powf(10.,0.2*tmp-2.);
    }
    sHGP[3] += delvr*Gau();
    sHGP[4] += delpm*Gau();
    sHGP[5] += delpm*Gau();
    to << sHGP << ' ' << deldm << ' ' 
       << delvr << ' ' << delpm << "\n"<< std::flush;
  }


}
