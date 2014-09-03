/***************************************************************************//**
\file Err.h
\brief Error handling code. Barely used. Very much antiquated. Should definitely be brought up to date.

*                                                                              *
*  Err.h                                                                       *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
/*

Barely used, in truth

*/

#ifndef _Torus_Err_def_
#define _Torus_Err_def_ 1

#include <cstdlib>
#include <string>
using std::string;

#undef _exception_handling_   // change this if your compiler does allow
			      // for C++ exception handling

struct TorusExcept {
  string msgs;
  int   index;
  //TorusExcept(char* m, int i) : msgs(string(m)), index(i) {}
  TorusExcept(string m, int i) : msgs(m), index(i) {}
};

#ifdef  _exception_handling_

//inline void TorusError(char* m, int i)
//{
//    throw TorusExcept(m,i);
//}

inline void TorusError(string m, int i)
{
    throw TorusExcept(m,i);
}

#else

#include <iostream>
using std::cerr;

#ifndef _toruserrno_
#define _toruserrno_ 1

extern int toruserrno;

// global variable indicating the error occured (not very elegant but most
// C++ compilers do currently not support exception handling ...)
// meaning of values:  0  anything but the following
//                     1  error in auxiliary functions (in Aux.h)
//                     2  error in classes IM_par or IsoMap (in Iso.h)
//                     3  error in classes PT_par or PointTrans (in Poi.h)
//                     4  error in classes GF_par, GenFunc, GenFuncFit,
//                        or AngMap (in Gen.h)
//                     5  error in potential

#endif

inline void TorusError(string m, int i)
{
    toruserrno = i;
    if(i<0) {
        cerr<<" ERROR: "<<m<<'\n';
        exit(1);
    } else
        cerr<<" WARNING: "<<m<<'\n';
}

#endif

#endif
