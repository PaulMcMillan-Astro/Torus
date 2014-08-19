/*******************************************************************************
*                                                                              *
* Types.h                                                                      *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 34P, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#ifndef _TorusTypes_
#define _TorusTypes_ 1

#include "Vector.h"

typedef Vector<double,3> Frequencies;
typedef Vector<double,3> Actions;
typedef Vector<double,3> Angles;
typedef Vector<double,3> Position;
typedef Vector<double,3> Velocity;
typedef Vector<double,4> Errors;

typedef Vector<double,2> vec2;
typedef Vector<double,3> vec3;
typedef Vector<double,4> vec4;
typedef Vector<double,6> vec6;

#include "PSP.h"

typedef double*                 Pdble;

#endif
