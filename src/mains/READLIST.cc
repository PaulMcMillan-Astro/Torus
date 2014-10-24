/***************************************************************************//**
\file READLIST.cc
\brief Gives full description of a torus on a list of tori

*                                                                              *
* READLIST.cc                                                                  *
*                                                                              *
* C++ code written by Paul McMillan, 2007-                                     *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>

#include "Torus.h"
#include "PJMebf.h"
#include "PJM_cline.h"


using std::cout;
using std::cin;

int main(int argc,char *argv[])
{
  int      temp;
  Torus    *T, T0;

//---------------------------------------------------------------------------
  if(argc<2){
    cerr << "Input: TorusList_file (Torus_name)\n"
	 << "READLIST gives full information about one Torus in a list\n";
    exit(1);
  }
//---------------------------------------------------------------------------

  std::string Tlist = string(argv[1]);
  std::vector<std::string> tori;
  if(Ebf_GetToriNames(Tlist,tori)) {
    cerr << "Not a torus list\n"; exit(1);
  }

  T = new Torus;
  if(argc >2) {
    for(int i = 2;i!=argc;i++) { 
      //temp = atoi(argv[i]);
      if(!(T->read_ebf(Tlist,string(argv[i])))) {
	cerr << "Torus "+string(argv[i])+"is not in the list.\n";
      } else { 
	//Tlist.ExTorus(temp, *T);
	T->show(cout);
      }    
    } 
  } else {
    for(;;){
      string tname;
      cout << "There are " << tori.size() 
	   << " tori in the list, probably t1 to t" << tori.size() 
	   << ", which one do you want to know about?\n";
      cin >> tname;
      if ( !(T->read_ebf(Tlist,tname))) break;
      T->show(cout);
    }
  }
}
