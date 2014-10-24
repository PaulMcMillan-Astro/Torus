/***************************************************************************//**
\file GiveListProps.cc
\brief Gives the properties of a TorusList file (e.g. actions of all tori)

*                                                                              *
* GiveListProps.cc                                                             *
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
#include <vector>
#include <string>

#include "PJM_cline.h"
#include "Torus.h"
#include "ebf.hpp"

using std::cout;
using std::ofstream;
using std::setw;


int Ebf_GetTagNames(const std::string& filename1, 
		    std::vector<std::string> &tagnames)
{
  using namespace std;
  ebf::ebfc::EbfHeader ebfh;
  tagnames.resize(0);
  FILE* fd = fopen(filename1.c_str(), "rb");
  if (fd != NULL)
    {
      fseek(fd, 0, SEEK_END);
      int64_t loc = ftell(fd);
      fseek(fd, 0, SEEK_SET);
      while (ftell(fd) < loc)
	{
	  ebf::ebfc::EbfHeader_read(&ebfh,fd);
	  fseek(fd,ebfh.capacity_,SEEK_CUR);
	  tagnames.push_back(ebfh.name);
	  if(ebfh.ecode!=0)
	    {
	      cout<<"Ebf error: in reading header"<<endl;
	      fclose(fd);
	      return 1;
	    }
	}
      fclose(fd);
    }
  else
    {
      cout<<"Ebf error: file not found"<<endl;
    }
  return 0;
}

int main(int argc,char *argv[])
{ 
  int i;
  Actions J;
  Torus T;
  
  ofstream to;
//---------------------------------------------------------------------------
  if(argc<4){
    cerr << "Input: TorusList_name output_name(- for stdout) "
	 << "output_type(s), options:\n";
    cerr << " E - energy\n";
    cerr << " J - actions\n";
    cerr << " O - Frequencies\n";
    cerr << " d - dJ/J\n";    
    cerr << " e - full error output\n";
    cerr << " I - Toy (Isochrone) map parameters\n";
    cerr << " n - Number of S_n terms\n";
    cerr << " T - the full Torus output\n";
    cerr << " p - peri- and apo-centre\n";
    cerr << " P - Parameters of the Point Transform\n";
    cerr << " f - Whatever\'s in the df container\n";
    cerr << " z - z_max\n";    
    cerr << "\nGiveListProps outputs properties of all tori in a given list\n";
    exit(1);
  }
//---------------------------------------------------------------------------

  std::vector<string> stuff, tori;
  ebf::EbfDataInfo dinfo;
  string filename = string(argv[1]);
  Ebf_GetTagNames(filename, stuff);
  for(int i=0;i!=stuff.size();i++) { //cerr << stuff[i] << '\n';
    if(stuff[i].size()>7)
      if(stuff[i].substr(stuff[i].size()-7,7)=="actions") 
	tori.push_back(stuff[i].substr(0,stuff[i].size()-7));
  }
  for(int i=0;i!=tori.size();i++) {//cerr << tori[i] << '\n';
    if(ebf::ContainsKey(filename,tori[i]+string(argv[3]),dinfo)) {
      double output[dinfo.elements];
      ebf::Read(filename,tori[i]+string(argv[3]),output,dinfo.elements);
      for(int i=0;i!=dinfo.elements;i++)
	cerr << output[i] << ' ';
      cerr << '\n';
    }
  }
}
