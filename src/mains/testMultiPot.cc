/***************************************************************************//**
\file testMultiPot.cc
\brief Tests for the MultiPot class


*                                                                              *
*  testMultiPot.cc                                                             *
*                                                                              *
*  C++ code written by Paul McMillan, 2014,                                    *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "MultiPot.h"
#include "MiyamotoNagaiPot.h"
#include "LogPot.h"
#include "falPot.h"
#include "StackelPot.h"
#include "IsochronePot.h"
#include <cmath>

using std::cout;
using std::cerr;
using std::ifstream;


int main(int argc,char *argv[])  
{
  double P, dPdR, dPdz, R=8., z=0.2;
  Frequencies KNO;
  
  Potential **PotList = new Potential* [3];
  
  //MiyamotoNagaiPotential MNPot1(1e10,5,0.5);
  //MiyamotoNagaiPotential MNPot2(2e10,2.,1.);
  //LogPotential LGPot(240*Units::kms, 0.8, 1., 0.);
  
  //PotList[0] = &MNPot1;
  //PotList[1] = &MNPot2;
  //PotList[2] = &LGPot;

  PotList[0] = new MiyamotoNagaiPotential(1e10,5,0.5);
  PotList[1] = new MiyamotoNagaiPotential(2e10,2.,1.);
  PotList[2] = new LogPotential(240*Units::kms, 0.8, 0.);
  
  MultiPotential Phi(PotList,3);

  cout<< "Object Phi, of type `MultiPotential' constructed. Now testing.\n";
  cout<< "It is made of 2 `MiyamotoNagaiPotential's and 1 LogPotential\n";
  cout<< "You can also add any GalaxyPotential(s), \n";
  cout<< "IsochronePotentials, or user defined potential(s)\n"; 
  cout<<"\n We work at an (approximate) solar position, R=8.0, z=0.007\n";
  cout<<" MultiPotential can return just the potential at that point\n";
  cout<<" P = Phi(R,z) = ";
  P=Phi(R,z);
  cout<< P <<"\n";
  cout<< "Or the potential and derivatives (not the forces), dP/dR & dP/dz\n";
  cout<<" P = Phi(R,z,dPdR,dPdz) = ";
  P=Phi(R,z,dPdR,dPdz);
  cout<< P <<"\n dPdR = " << dPdR << "\n dPdz = " << dPdz << '\n';
  cout<<"\n Other possible queries are\n";
  cout<<" R = Phi.RfromLc(Lz), ";
  cout<<"the radius of a circular orbit with Ang Mom Lz\n";
  cout<<" Lz = Phi.LfromRc(R), "<<Phi.LfromRc(R) <<", the inverse\n";
  cout<<" KNO = Phi.KapNuOm(R), epicycle frequencies, radial, vertical & azimuthal\n"; 
}
