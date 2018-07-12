/***************************************************************************//**
*                                                                              *
\file Point_ClosedOrbitCheby.h
\brief Contains class PoiClosedOrbit.
Code for the point transform used when Tori have very low J_R -- i.e. are
nearly closed orbits in R-z plane.

*                                                                              *
* C++ code written by Paul McMillan, 2008, modified by James Binney 2017       *
* e-mail: paul@astro.lu.se                                                     *
* github: https://github.com/PaulMcMillan-Astro/Torus                          *
*                                                                              *
*******************************************************************************/
//
// Header file for a point transformation.
// Input (r,TH,{phi},pr,pTH,{pphi}), Output (R,z,{phi},vR,vz,{v_phi})
//


#ifndef _PoiCh_tr_
#define _PoiCh_tr_ 1

#include "Potential.h"
#include "Maps.h"
#include "Types.h"
#include "CHB.h"
#include "Toy_Isochrone.h"
////////////////////////////////////////////////////////////////////////////////

/**

\brief Point transform used when fitting orbits with J_r << J_l, which
cannot be fitted otherwise. The transform is defined by:

 r^T  = x(th)*r
 th^T = y(r)*z(th)

 This is done so that r^T = a is constant for a J_R = 0 orbit (as required), and
 then fixes th^T so that pth^T has the correct dependence on th^T

------------------------------------------------------------------------------

 The class PoiClosedOrbit contains routines which implement this point transform
 and ones which find it for a given J_l, J_phi and gravitational potential.

 1. Finding the transform:
    Done by the function set_parameters(Potential*, Actions, a), where
    a is the radius of the toy shell orbit. This finds the
 correct closed (J_R=0) orbit in the target potential, and determines the
 functions x, y, and z (above), storing them as Chebyshev polynomials.

 2. Performing the transform:
    Done in the usual way, with Forward and Backward. x,y and z (above) found
 from the Chebyshev polynomials (for th<thmax) or a quadratic for x and z if
 th>thmax - this ensures that x and z are accurately fitted, and don't trend
 off to extreme values for th>thmax.  Note that we actually store z' = z/th


*///////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
class PoiClosedOrbit : public PoiTra {
 private:
  // properties of the point transform
  double Jl,Lz,alt2,thmax,omz,
    ax,bx,cx, az,bz,cz;   // coefficients of quadratic for x and z for th>thmax
  Cheby xa,ya,za;         // Chebyshev polynomials which define transform
                          // N.B. za stores z' = z/th because that's ~constant
  double thmaxforactint;  // used to find the action
  Cheby vr2, drdth2, pth2;//         ''
  double actint(double) const;
  double get_x(double,double&,double&) const;
  double get_z(double,double&,double&) const;
// ROUTINES USED TO FIND THE TRANSFORM -----------------------------------------
  // Integrate orbit over Pi in th_z & store output
  void do_orbit       (PSPD, double, Potential*, double*, double*, double*,
		       double*, double*, double*, double*, int&, int);
  // Shift starting radius towards that of closed orbit (or note we're there)
  void set_Rstart     (double&,double,double&,double&,bool&,bool&,bool&);
  // Shift orbit energy so that closed orbit of that energy has J_l = target
  void set_E          (const double,double&,double&,double&,bool&,bool&,bool&);
  // Rewrite tables so I can fit to chebyshev (also find thmax)
  void RewriteTables  (const int, double *, double *, double *, double *,
		       double *, double *, double *, double *, int &);
  // Routines needed to find y and z
  vec2 chebderivs     (const double, const vec2);
  double stepper      (vec2&, vec2&, const double, const double );
  void yzrkstep       (vec2&, vec2&, const double, double&, const double,
		       double&, const double, const double) ;
  double yfn          (double, vec2 *, const int, double * , int);
  void plot_grid      (int,double,double,double,double,double*,double*,int);
  void plot_SoS(double*,double*,const int,double,ToyMap*);
//------------------------------------------------------------------------------
 public:
  // various constructors
  PoiClosedOrbit();
  PoiClosedOrbit(Cheby, Cheby, Cheby, double);
  PoiClosedOrbit(const double*);
  PoiClosedOrbit(Potential *, const Actions, ToyMap *);
  ~PoiClosedOrbit() {}
  void    set_parameters    (Cheby, Cheby, Cheby, double);
  void    set_parameters    (Potential *, const Actions, ToyMap *);
  void    set_parameters    (const double*);
  void    parameters	    (double *)	                 const;
  int     NumberofParameters()                           const;
  PSPD    Forward           (const PSPD&)                const;
  PSPD    Backward          (const PSPD&)                const;
  PSPT    Forward3D         (const PSPT&)                const;
  PSPT    Backward3D        (const PSPT&)                const;
  void    Derivatives       (double[4][4])               const;
};

#endif
