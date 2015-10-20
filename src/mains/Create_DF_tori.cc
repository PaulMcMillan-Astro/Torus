/***************************************************************************//**
\file Create_df_tori.cc
\brief Generates tori given a potential and output from Choose_any_df

*                                                                              *
* Create_DF_tori.cc                                                            *
*                                                                              *
* C++ code written by Paul McMillan, 2011-                                     *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/


#include <iostream>
#include <iomanip>
#include <fstream>
#include "Torus.h"
#include "ebf.hpp"
#include "PJM_cline.h"
#include "falPot.h"
#include "LogPot.h"

int main(int argc,char *argv[])
{
  if(argc<5) {
    cerr << "Input: potfile J&df&n_file output nout (Rmin=) (Rmax=) (zmin=)\n";
    exit(0);
  }
  Actions   J,oJ;
  int       nout=atoi(argv[4]), nMC;
  double    Rmin=0., Rmax=1.e99, zmin=0., df;
  Torus     T;
  //TorusList Tlist(argv[3],1); //write only
  string Tlist = string(argv[3]);
  ebf::WriteString(Tlist,"/log","Torus list","w");
  ifstream  from;
  Potential *Phi;

  for(int i=5;i<argc;i++) {
    bool understood=false;
    if(parse_comm_line(argv[i],"Rmin=",Rmin))        understood=true;
    if(parse_comm_line(argv[i],"Rmax=",Rmax))        understood=true;
    if(parse_comm_line(argv[i],"zmin=",zmin))        understood=true;
    if(!understood) {
      cerr << "Input "<<argv[i]<<" not understood\n"; 
      exit(0);
    }
  }



  if(string(argv[1]) == "LogPotential_220") {
    Phi = new LogPotential(220.*Units::kms,0.8,0.);
  } else {
    my_open(from,argv[1]);
    Phi = new GalaxyPotential(from);
    from.close();
  }

  my_open(from,argv[2]);

  for(int i=0; i!=nout && !from.eof();) {
    from >> J >> df >> nMC;
    //J[2] = -J[2]; // uncomment if you want galaxy rotating in positive sense
    //cerr << J << '\n' << std::flush;
    T.AutoFit(J,Phi,0.003); // Fit Torus from scratch
    T.FindLimits();         // To check torus enters relevant region
    if(T.minR()<Rmax && T.maxR()>Rmin && T.maxz()>zmin) {
      std::stringstream tor; tor << "t" << i+1;
      T.write_ebf(Tlist,tor.str());
      ebf::Write(Tlist,"/"+tor.str()+"/auxil/nMCMC",&nMC,"a","",1);
      i++;
    }
  }

}
