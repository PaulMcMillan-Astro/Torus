/***************************************************************************//**
\file Potential.h \brief Contains class Potential. Code for the base
for all user defined potentials -- including those supplied by the
original writer (i.e. me. Hello.). 

Note that Potential can provides the effective potential at a given Lz
as well as the true potential.

*                                                                              *
*  Potential.h                                                                 *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2006/07                                   *
* e-mail: paul@astro.lu.se                                                     *
* github: https://github.com/PaulMcMillan-Astro/Torus                          *
*                                                                              *
*******************************************************************************/

#ifndef _Pot_def_
#define _Pot_def_ 1

#include "Types.h"
#include <cmath>
#include <cstdlib>

/**
   \brief Base class for user defined potentials (or those provided).

   The user defines a class derived from Potential, which provides the
   operator() with 2 and 4 arguments (the latter for derivatives),
   RfromLc, LfromRc and KapNuOm. An example is the logarithmic
   potential as implemented in class LogPotential.
*/

class Potential {
protected:
    double Lzsq;

public:

    void set_Lz    (const double Lz)    { Lzsq = Lz*Lz; }

    Potential      (const double Lz=0.) { set_Lz(Lz); }

    double Lzsquare() const             { return Lzsq; }
    /** \brief returns Phi given (R,z)                                  */
    virtual double operator()(				// returns Phi(R,z)
			      const double,		// given R
			      const double)const=0;	// and z

    /** \brief returns Phi given (R,z), also computes and returns 
    dPhi/dR and dPhi/dz
    */
    virtual double operator()(				// returns Phi(R,z)
			      const double,		// given R
			      const double,		// and z
			      double&,			// also computes dPhi/dR
			      double&)const=0;		// and dPhi/dz
    /** \brief returns radius of a circular orbit (R) with given
	angular momentum (Lc). Optionally also gives dR/dLc
     */
    virtual double RfromLc   (				// returns Rc,
		      const double,		// given Lz. possibly
		      double* =0)const;	// returns dRc/dLz.

    /** \brief returns angular momentum of a circular orbit (L) with given
	 radius (Rc). Optionally also gives dL/dRc
     */
    virtual double LfromRc   (				// returns Lc,
		      const double R,		// given R. possibly
		      double* dR=0) const {        	// returns dLz/dRc.
      double dPR,dPz,P;
      P = (*this)(R,0.,dPR,dPz);
      return sqrt(R*R*R*dPR);  
    }

    /** \brief Gives epicycle frequencies (kappa, nu, Omega - i.e. radial, 
	vertical, azimuthal) at a given R.
     */
    virtual Frequencies KapNuOm(			// returns kappa,nu,Om
				const double)const=0;	// given R at z=0

    /** \brief effective potential at R,z - i.e. Phi + Lz^2/(2R^2) 
     */
    double eff(const double R, const double z) const
	{ if(Lzsq==0.) return (*this)(R,z);
          if(R==0.) {
	      cerr << " error in class Potential::eff: R=0 at non-zero Lz\n";
              exit(1); }
          return (*this)(R,z) + 0.5 * Lzsq/(R*R); }

    /** \brief effective potential at R,z - i.e. Phi + Lz^2/(2R^2), and 
	dPhi_eff/dR, dPhi_eff/dz
     */
    double eff(const double R, const double z, double& dPdR, double& dPdz) const
	{ 
	  if(Lzsq==0.) return (*this)(R,z,dPdR,dPdz);
    	  if(R==0.) {
              cerr << " error in class Potential::eff: R=0 at non-zero Lz\n";
              exit(1); }
    	  register double potential     = (*this)(R,z,dPdR,dPdz);
   	  register double Lzsq_over_Rsq = Lzsq/(R*R);
   	  dPdR                         -= Lzsq_over_Rsq / R;
    	  return potential + 0.5 * Lzsq_over_Rsq; }
};

inline double Potential::RfromLc(const double L, double* dR) const
{
  bool more=false;
  double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL;
  R=exp(lR);
  P= (*this)(R,0.,dPR,dPz);
  LcR=pow(R*R*R*dPR,0.5);
  if(LcR == L) return R;
  if(L>LcR) more=true;
  oldL=LcR;
  
  for( ; ; ){
    lR += (more)? dlR : -dlR;
    R=exp(lR);
    P= (*this)(R,0.,dPR,dPz);
    LcR=pow(R*R*R*dPR,0.5);
    if(LcR == L) return R;
    if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
	return R;}
    oldL=LcR;
  }
  
}


#endif
