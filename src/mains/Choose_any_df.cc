/***************************************************************************//**
\file Choose_any_df.cc
\brief Chooses points from a df (see DF.h) using a Metropolis algorithm. Output actions, df value & corresponding weight (from MCMC).

*                                                                              *
* Choose_any_df.cc                                                             *
*                                                                              *
* C++ code written by Paul McMillan, 2008                                      *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*       Produce tori with actions MCMC taken from simple df                    *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using std::cout;

#include "Torus.h"
#include "Random.h"
#include "falPot.h"
#include "LogPot.h"
#include "PJM_cline.h"
#include "DF.h"

class tunableMCMC {
public:
  DF *f;
  Potential *Phi;
  Actions oJ, oJ_pr, sigmas, proposal;  
  double stepsize,odf;
  Gaussian *Gau;
  Random3 *R3;
  bool step();
  void burn_in(int);
  void find_sigs(int);
  void tune(int,int);
  void output(int,ofstream&);
  void output(int,std::vector<Actions>*,std::vector<int>*);
  tunableMCMC(DF*, Potential*, Actions, Gaussian*, Random3*);
  ~tunableMCMC() {;}
};

tunableMCMC::tunableMCMC(DF *F, Potential *phi, Actions J, 
			 Gaussian *gau, Random3 *r3) : 
  f(F), Phi(phi), oJ(J), Gau(gau), R3(r3)
{
  oJ_pr    = oJ;
  oJ_pr[0] = sqrt(oJ[0]);
  oJ_pr[1] = sqrt(oJ[1]);
  odf = f->df(Phi,oJ)*fabs(oJ_pr(0)*oJ_pr(1));
  stepsize = 1.;
  proposal = 0.1; proposal[2] = 0.3;
}

bool tunableMCMC::step() {
  Actions J,J_pr;
  double df;
  for(int j=0;j!=3;j++) {
    J_pr[j] = oJ_pr[j] + stepsize*proposal[j]*(*Gau)();
    J[j] = (j<2)? J_pr(j)*J_pr(j) : J_pr(j);
  }
  Phi->set_Lz(fabs(J(2)));
  df = f->df(Phi,J)*fabs(J_pr(0)*J_pr(1));  
  if(df>odf*(*R3)()) {
    odf = df; oJ = J; oJ_pr = J_pr; return true;
  }
  return false;
}

void tunableMCMC::burn_in(int nstep) {
  for(int i=0;i!=nstep;)  if(step()) i++;
}

void tunableMCMC::find_sigs(int nstep) {
  Actions mean=0.;
  sigmas = 0.;
  int ntry=0;
  for(int i=0;i!=nstep;) {
    if(step()) i++;
    mean += oJ_pr; ntry++;
    for(int j=0;j!=3;j++) 
      sigmas[j] += oJ_pr[j]*oJ_pr[j];
  }
  mean   /= double(ntry);
  for(int j=0;j!=3;j++) 
    sigmas[j] = sqrt(sigmas[j]/double(ntry) - mean[j]*mean[j]);
  proposal = sigmas;
}

void tunableMCMC::tune(int nstep, int nit) {
  int ntry;
  for(int i=0;i!=nit;i++) {
    ntry=0;
    for(int j=0;j!=nstep;) {
      if(step()) j++;
      ntry++;
    }
    double frac = double(nstep)/double(ntry);
    if(3*frac<1) 
      stepsize *= 1-pow(3*frac-1.,2)*0.8;
    else
      stepsize *= 1+pow(3*frac-1.,2)*0.75;
  }
}

void tunableMCMC::output(int nout,ofstream &to) {
  int ntry=0;
  double Rtmp;
  to << oJ << ' ' << odf/fabs(oJ_pr(0)*oJ_pr(1)) << ' ';
  for(int i=0;i!=nout;) {
    ntry++;
    if(step()) {
      i++;
      to << ntry << "\n" << std::flush;
      ntry = 0;
      if(i != nout) 
	to << oJ << ' ' << odf/fabs(oJ_pr(0)*oJ_pr(1)) << ' ' << std::flush;
    }
  }
}

void tunableMCMC::output(int nout,std::vector<Actions> *Jout, std::vector<int> *Wout) {
  int ntry=0;
  double Rtmp;
  Jout->push_back(oJ);
  for(int i=0;i!=nout;) {
    ntry++;
    if(step()) {
      i++;
      Wout->push_back(ntry);
      ntry = 0;
      if(i != nout) Jout->push_back(oJ);
    }
  }
}


int MCMC(DF *f, Potential* Phi, Actions &oJ_pr, double& odf, 
	 Gaussian &Gau, Random3 &R3) 
{
  int ntry=0;
  Actions J,J_pr;
  double df;
  do {
    for(int j=0;j!=2;j++) {
      J_pr[j] = oJ_pr[j] + 0.1*Gau();
      J[j] = J_pr(j)*J_pr(j);
    }
    J[2] = oJ_pr(2) + 0.3*Gau();
    df = f->df(Phi,J)*fabs(J_pr(0)*J_pr(1));
    ntry++;
  } while (df<odf*R3.RandomDouble());
  odf = df; oJ_pr = J_pr; oJ_pr[2] = J[2]; 
  return ntry;
}


int main(int argc,char *argv[])
{

  if(argc<5) { 
    cerr << "Produces list of actions, df value and N_MCMC "
	 << "sampled using MCMC from df\n";
    cerr << "Input parameters: Potfile df_file out_rootname ntor (seed)\n";
    exit(0);
  }

  ifstream  from;
  ofstream  to;
  time_t    cpu=time(NULL);
  int       flag, ntor, seed=0;
  Actions   J,rtJ,ortJ,J_pr;
  Potential *Phi;
  if(argc>5) seed = 127949*atoi(argv[5]);
  
  Random3   R3(7*int(cpu) + seed),R3b(123*int(cpu)+seed);
  Gaussian  Gau(&R3,&R3b);
  double    odf,newdf;
  string logname = string(argv[3]) + ".CdJ_log";

  if(string(argv[1]) == "LogPotential_220") {
    Phi = new LogPotential(220.*Units::kms,0.8,0.);
  } else {
    my_open(from,argv[1]);
    Phi = new GalaxyPotential(from);
    from.close();
  }

  my_open(from,argv[2]);
  DF *distfunc = set_DF(from);  
  from.close();

  ntor = atoi(argv[4]);


  // Write the log
  struct tm * timeinfo = localtime(&cpu);
  to.open(logname.c_str());
  string line;
  to << "Output file: " + string(argv[3]) + "\n"
     << "Created at: "  << asctime(timeinfo)
     << "By program: " << argv[0] << "\n"
     << "Input parameters: potfile df_file ntor "
     << "Ran_param (x2)\n" 
     <<  argv[1] << " " <<  argv[2] << " " << ntor  <<" "
     << 7*int(cpu) << " " << 123*int(cpu) << "\n"
     << "with potfile contents:\n";
  if(string(argv[1]) == "LogPotential_220") {
    to <<argv[1]<< "\n";
  } else { 
    from.open(argv[1]);
    while(from && !from.eof()) {
      getline(from,line);
      to << line << '\n';
    }
    from.close();
  }
  to << "and df_file contents:\n";
  from.open(argv[2]);
  while(!from.eof()) {
    getline(from,line);
    to << line << '\n';
  }
  from.close();

  to.close();

  my_open(to,argv[3]);
  J[0] = J[1] = 0.001; // for example
  J[2] = -1.8; 

  tunableMCMC tMC(distfunc,Phi,J,&Gau,&R3);

  int ntune = ntor/20, nburn = ntor/3, nit = 4;
  if(ntune<100) ntune = 100;
  nburn -= (nit+1)*ntune;
  if(nburn<10) nburn = 10;
  tMC.find_sigs(ntune);
  tMC.tune(ntune,nit);
  tMC.burn_in(nburn);  
  tMC.output(ntor,to);  
  to.close();

}
