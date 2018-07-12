
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "Orb.h"
#include "Torus.h"
#include "falPot.h"
#include "Random.h"
#include "PJMCoords.h"
#include "PJM_cline.h"

using std::cout;
using std::cerr;
using std::cin;
using std::fstream;


int main(int argc,char *argv[])
{
  clock_t cpu0=clock(); 
  int jact,typex,typey,sizeT,sizeO=0,F,NT=50, NO=50,
    Ndyn=15,Nran=10000,n=0,Nlarge=500000,flag;
  time_t second=time(NULL);
  long cp=second;
  float ti[Nlarge],Ri[Nlarge], zi[Nlarge],phii[Nlarge], xi[Nlarge], yi[Nlarge];
  double Oz;
  PSPD StPo, oJT, JT;
  GCY *tabGCY;
  Actions J;
  Torus    *T, T0;
  string basename;
  fstream sosin;
  ofstream sosout;
  ifstream from;
  ofstream to;
  Potential *Phi;
  Random3 r3(cp);
  bool done = false,curs = false, first = true, both=false, rando=false;

  T = new Torus;
 
//---------------------------------------------------------------------------
  if( argc !=6 && argc !=7  && argc !=8 ){
    cerr << "Outputs SOS from torus and 5 different orbit integrations.\n"
	 <<"Input: potential JR Jz Jphi output_file (dJ/J)\n";
    exit(1);
  }
//---------------------------------------------------------------------------
//===========================================================================
// 1. Get target potential
  my_open(from, argv[1]);
  Phi = new GalaxyPotential(from);
  from.close();

// input actions
  J[0] = atof(argv[2]); J[1] = atof(argv[3]); J[2] = atof(argv[4]);
  double dJ = (argc==7)? atof(argv[6]) : 0.003;
  bool err = false;

// Fit the Torus -------
  flag = T->AutoFit(J,Phi,dJ);
  Phi->set_Lz(J(2));

  cerr << T->omega() << '\n';

  
  const char* cfile1 ="/tmp/tmp.sos";
  sosout.open(cfile1);
  T->SOS(sosout,NT);
  sosout.close();
  sosin.open(cfile1,fstream::in);
  float zT,pzT,trT;
  for (int i=0;i!=NT;++i){
    sosin >> xi[i] >> zT >> yi[i] >> pzT >> trT;
  } 
  sosin.close();
  sizeT=NT;
  std::remove("/tmp/tmp.sos");
  float xT[sizeT],yT[sizeT];

  for(int i=0; i!=sizeT; i++) {
      xT[i] = xi[i]; yT[i] = yi[i];
  }

// The orbit integration part --------------------------------------------------
  const char* cfile2  ="/tmp/tmp.oss";  
  const char* cfile2a ="/tmp/tmp.ossb";
  //  for(double tr=Pi/60.;tr<Pi-0.01;tr+=Pi/60.) {
  for(double tr=Pi/30.;tr<Pi-0.0001;tr+=Pi/30.) {

    // cerr << tr << '\n';
    {
      //double d1=r3.RandomDouble(), d2=r3.RandomDouble();
      Angles Atmp=0.; Atmp[0] = tr;
      Atmp[0] = r3.RandomDouble()*TPi;
      Atmp[1] = r3.RandomDouble()*TPi;
      //double d1=tr , d2=0.;
      StPo = T->Map(Atmp); //T->StartPoint(d1,d2);
    }
    F = orbit(Phi,StPo,1.e-12,NO,1.e-3*StPo(0),cfile2a,cfile2,Oz);
    
    sosin.open(cfile2,fstream::in);
    for (int i=sizeO;i!=sizeO+NO;++i){
      sosin >> xi[i] >> yi[i];
    } 
    sizeO += NO;
    sosin.close();
   
    std::remove("/tmp/tmp.oss");
    std::remove("/tmp/tmp.ossb");
  }

  float xO[sizeO],yO[sizeO];
  for(int i=0; i!=sizeO; i++) {
    xO[i] = xi[i]; yO[i] = yi[i];
  }



  string TorOut = string(argv[5]) + ".torSOS",
    OrbOut = string(argv[5]) + ".intSOS";
  
  my_open(to,TorOut);
  
  for(int i=0;i!=NT;i++) {
    to << xT[i] << ' '<< yT[i] << '\n';
  }
  for(int i=0;i!=NT;i++) {
    to << xT[NT-1-i] << ' '<< -yT[NT-1-i] << '\n';
  }
  to.close();

  my_open(to,OrbOut); 
  for(int j=0;j!=29;j++) {
    std::stringstream tmp;
    tmp << OrbOut << '.' << j;
    string OUT;
    tmp >> OUT;
    //my_open(to,OUT);
    
    for(int i=0;i!=NO;i++) {
      to << xO[j*NO+i] << ' '<< yO[j*NO+i] << '\n';
    }
    //to.close();
  }
}
