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
  bool tocout = false;
  ofstream to;
//---------------------------------------------------------------------------
  if(argc<4){
    cerr << "Input: TorusList_name output_name(- for stdout) "
	 << "output_type(s), options:\n";
    cerr << " E - energy\n";
    cerr << " actions - actions\n";
    cerr << " frequencies - Frequencies\n";
    cerr << " d - dJ/J\n";    
    cerr << " errors - full error output\n";
    cerr << " toymap/isochrone - Toy (Isochrone) map parameters\n";
    cerr << " n - Number of S_n terms\n";
    cerr << " Rmin - pericentre\n";
    cerr << " Rmax - apocentre\n";
    cerr << " zmax - maximium height above the plane\n";
    cerr << " P - Parameters of the Point Transform\n";    
    cerr << "\nGiveListProps outputs properties of all tori in a given list\n";
    exit(1);
  }
//---------------------------------------------------------------------------

  std::vector<string> stuff, tori;
  ebf::EbfDataInfo dinfo;
  string filename = string(argv[1]);

  if(string(argv[2])=="-") tocout = true; 
  else to.open(argv[2]);

  Ebf_GetTagNames(filename, stuff);

  for(int i=0;i!=stuff.size();i++) { //cerr << stuff[i] << '\n';
    if(stuff[i].size()>7)
      if(stuff[i].substr(stuff[i].size()-7,7)=="actions") 
	tori.push_back(stuff[i].substr(0,stuff[i].size()-7));
  }

  for(int i=0;i!=tori.size();i++) {//cerr << tori[i] << '\n';
    for(int j=3;j<argc;j++) {
      if(string(argv[j]) == "n") {
	ebf::ContainsKey(filename,tori[i]+"generatingfunction/N1",dinfo);
	if(tocout) std::cout << dinfo.elements << ' ';
	else to << dinfo.elements << ' ';
      } else if(string(argv[j]) == "d") {
	double dc[4];
	ebf::Read(filename,tori[i]+"errors",&dc[0],4);
	if(tocout) std::cout << dc[0] << ' ';
	else to << dc[0] << ' ';
      }
      if(ebf::ContainsKey(filename,tori[i]+string(argv[j]),dinfo)) {
	double output[dinfo.elements];
	ebf::Read(filename,tori[i]+string(argv[j]),output,dinfo.elements);
	for(int i=0;i!=dinfo.elements;i++) {
	  if(tocout) std::cout << output[i] << ' ';
	  else to << output[i] << ' ';
	}
      }
    }
    if(tocout) std::cout  << '\n';
    else to  << '\n';
  }

  if(string(argv[2])!="-") to.close();
}
