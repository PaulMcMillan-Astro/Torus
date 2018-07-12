/***************************************************************************//**
\file Sample_list_limits.cc
\brief Generates points from a TorusList within given limits. Takes output from Create_DF_tori. Though I should fix it so that it doesn't have to.

*                                                                              *
* Sample_list_limits.cc                                                        *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2007                                      *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
using std::cout;


#include "Torus.h"
#include "Random.h"
#include "Numerics.h"
#include "PJMCoords.h"
#include "PJMebf.h"
#include "PJM_cline.h"

bool in_phi_range(double p, double pmin, double pmax) {
  if(p>=pmin && p<=pmax) return true;
  if(p>pmax) {
    double dp = p-pmin; int test=int(dp/TPi);
    if(p-TPi*test<pmax) return true;
    else return false;
  }
  
  // if p<pmin
  double dp = pmax-p; int test = int(dp/TPi);
  if(p+TPi*test>pmin) return true;
  else return false;
  
}


int main(int argc,char *argv[])
{

  if(argc<4) {
    cerr << argv[0] << " - Produces positions, velocities & actions "
	 << "from input tori in given limits\n"
	 << "Output: R,z,phi,vR,vz,vphi,Jr,Jz,Jphi,th_R,th_z,th_phi\n"
	 << "Input: TorusList outputfile (- for cout) "
	 << "n_out (Rmin=) (Rmax=) (zmin=) (zmax=) (phimin=) (phimax=)\n";
    exit(0);
  }

  ofstream  to;
  bool      to_cout=false;
  time_t    cpu=time(NULL);
  int       NT, n_torus, nout=atoi(argv[3]), inlimit;//, *useable;
  Random3   R3(7*int(cpu));
  Angles    A;
  Torus     *T = new Torus;
  PSPT      QP;
  double    Rmin=0., Rmax=1.e99, zmin=-1.e99, zmax=1.e99,phimin=0,phimax = 2*Pi,
    running_tot=0.,*rank;
  GCY sGCY;
  GCA sGCA;
  OmniCoords OC; // converts various desriptions of velocity

  for(int i=4;i<argc;i++) {
    bool understood=false;
    if(parse_comm_line(argv[i],"Rmin=",Rmin))        understood=true;
    if(parse_comm_line(argv[i],"Rmax=",Rmax))        understood=true;
    if(parse_comm_line(argv[i],"zmin=",zmin))        understood=true;
    if(parse_comm_line(argv[i],"zmax=",zmax))        understood=true;
    if(parse_comm_line(argv[i],"phimin=",phimin))    understood=true;
    if(parse_comm_line(argv[i],"phimax=",phimax))    understood=true;
    if(!understood) {
      cerr << "Input "<<argv[i]<<" not understood\n"; 
      exit(0);
    }
  }

  //TorusList Tlist(argv[1],0);                       // Input Tori
  std::string Tlist = string(argv[1]);
  std::vector<std::string> tori;
  if(Ebf_GetToriNames(Tlist,tori)) {
    cerr << "Not a torus list\n"; exit(1);
  }
  NT = tori.size();
  if(!NT) { cerr << "No tori in TorusList\n"; exit(1); } 
  if(argv[2][0] == '-') to_cout=true;               // output file/cout
  else to.open(argv[2]);                            //
  
  rank = new double[NT];

  for(int i=0;i!=NT;i++) {
    if(ebf::ContainsKey(Tlist,"/"+tori[i]+"/Rmin")) {// assume if one -> all

      if( ebf::ContainsKey(Tlist,"/"+tori[i]+"/auxil/nMCMC")) {
	double TRmin, TRmax, Tzmax;
	int nMCMC;
	ebf::Read(Tlist,"/"+tori[i]+"/Rmin",&TRmin,1);
	ebf::Read(Tlist,"/"+tori[i]+"/Rmax",&TRmax,1);
	ebf::Read(Tlist,"/"+tori[i]+"/zmax",&Tzmax,1);
	ebf::Read(Tlist,"/"+tori[i]+"/auxil/nMCMC",&nMCMC,1);
	if(!(TRmax<Rmin || TRmin>Rmax || Tzmax<zmin || -Tzmax>zmax))
	  running_tot += nMCMC; 
	rank[i] = running_tot;
      } else {
	double TRmin, TRmax, Tzmax;
	ebf::Read(Tlist,"/"+tori[i]+"/Rmin",&TRmin,1);
	ebf::Read(Tlist,"/"+tori[i]+"/Rmax",&TRmax,1);
	ebf::Read(Tlist,"/"+tori[i]+"/zmax",&Tzmax,1);
	if(!(TRmax<Rmin || TRmin>Rmax || Tzmax<zmin || -Tzmax>zmax))
	  running_tot += 1; 
	rank[i] = running_tot;
      }
    } else {
      cerr << "Torus " + tori[i] + " doesn't have all the info I need\n";
      exit(0);
    }
  }    


  for(int i=0;i!=nout;) {
    double pick = R3();
    int which = int(NT*pick);               // random with metropolis weighting
    pick   *= running_tot;                  // guess appropriate point in list
    n_torus = hunt(rank,NT,pick,which)+1;   // find appropriate point in list
    if(!(T->read_ebf(Tlist,tori[n_torus]))) {
      std::cerr << "Fail " + Tlist+ " " + tori[n_torus]+"\n"; 
    }
    for(int j=0;j!=3;j++) A[j] = TPi*R3();  // random angle
    QP = T->Map3D(A);
    if(QP(0)<Rmax && QP(0)>Rmin && QP(1)<zmax && QP(1)>zmin && 
       in_phi_range(QP(2),phimin,phimax)) {
      for(int j=0;j!=6;j++) sGCY[j] = QP(j); // change to appropriate container
      OC.take_GCY(sGCY);                     // put into converter
      if(to_cout)
	cout << OC.give_GCA() <<' '  // get Cartesian coords, kpc & km/s
	     <<T->actions()<< ' ' << A << '\n' << std::flush;
      else
	to   << OC.give_GCA() <<' ' // get Cartesian coords, kpc & km/s
	     << T->actions()<< ' ' << A << '\n' << std::flush;
      i++;
    }
  }

}
