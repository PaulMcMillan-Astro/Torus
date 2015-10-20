/***************************************************************************//**
\file AUTOTORUS.cc
\brief Fits a single torus from actions given at the command line

*                                                                              *
* AUTOTORUS.cc                                                                 *
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

#include "PJMebf.h"
#include "Torus.h"

#include "falPot.h"
#include "LogPot.h"
#include "MiyamotoNagaiPot.h"

#include "PJM_cline.h"

using std::cerr;



int main(int argc,char *argv[])
{
  int      flag,err=1;
  double   dJ=0.003;
  Actions  J,oJ;
  Torus    *T;
  ifstream from;
  ofstream to;
  string outname="";
  Potential *Phi;
  T = new Torus;
//---------------------------------------------------------------------------
  if(argc<5 || argc>8){
    cerr <<"Input: potential actions (out=) (dJ=) (err=)\n\n"
	 <<"This programme is a good way of finding tori.\n";
    exit(1);
  }
//---------------------------------------------------------------------------
//===========================================================================
// 1. Get target potential
  if(string(argv[1]) == "LogPotentialYST") {
    Phi = new LogPotential(220.*Units::kms,0.8,0.);
  } else {
    my_open(from,argv[1]);
    Phi = new GalaxyPotential(from);
    from.close();
  }
  //Phi = new MiyamotoNagaiPotential(1.*Units::G_i,1.,.45);
  //Phi = new LogPotential(244.5*Units::kms,0.7,0.,0.);
  
// input actions
  J[0] = atof(argv[2]); J[1] = atof(argv[3]); J[2] = atof(argv[4]);

  for(int i=5;i!=argc;i++) {
    bool understood=false;
    if(parse_comm_line(argv[i],"out=",outname)) understood=true;
    if(parse_comm_line(argv[i],"dJ=",dJ)) understood=true;
    if(parse_comm_line(argv[i],"err=",err)) understood=true;
    if(!understood) {
      cerr << "Input "<<argv[i]<<" not understood\n"; exit(1);
    }

  }
  
  // flag = T->AutoFit(J,Phi,dJ,700,300,15,5,24,200,24,err);
  flag = T->AutoFit(J,Phi,dJ,600,150,12,4,16,200,40,err);
  std::cout << J << ", flag = " << flag<< "\n\n" << std::flush;

  std::cout << T->minR() << ' '  << T->maxR() << ' ' << T->maxz() << '\n' 
	    << std::flush;

  // Output. 
  // flag = 0  -> success
  // flag = -1 -> catastrophic breakage
  // flag = -2 -> target < dH <2*target
  // flag = -3 -> 2*target < dH
  // flag = -4 -> angle fit failed

  if(outname != "") {
    ebf::WriteString(outname,"/log","Torus list","w");
    T->write_ebf(outname,"T1");
    //T->read_ebf(outname,"T1");
  }
    T->show(std::cout); // give ASCII output to screen.

}


