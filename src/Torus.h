/***************************************************************************//**
\file Torus.h
\brief Contains class Torus. The code that puts it all together.

*//*
 Torus.h

 C++ code written by Walter Dehnen, 1995-97,
                     Paul McMillan, 2007-
                     James Binney, 2015-
 e-mail: paul@astro.lu.se
 github: https://github.com/PaulMcMillan-Astro/Torus

 class Torus         A Torus is defined by the Actions and its parameters, it
	              gives mapping   Angle variables -> cylindrical co-ords.
		      Among other things, a Torus can fit itself, can make a
		      surface of section, can compute whether a space point is
		      ever hit, and if so with which velocity and probability
		      (density).


*******************************************************************************/

#ifndef _Torus_
#define _Torus_ 1

#include <fstream>
#include <iostream>
#include "GeneratingFunction.h"
#include "CHB.h"
#include "Vector.h"
#include "Matrix.h"
#include "Potential.h"
#include "PJMCoords.h"
#include "Fit.h"
#include "Types.h"
#include "Toy_Isochrone.h"
/**
\brief Class combining everything. Can fit a Torus with Actions J in a
given Potential. Then gives position and velocity given angles. Lots
of other fun stuff too.

This class is where everything gets combined together. It's the main event.

Fundamentally a Torus is defined by four objects, one each of classes AngMap
GenFnc ToyMap and PoiTra. Torus contains these objects and has the ability to
fit them for a given Potential. It can then use them to do the complete transformation from angle-action coordinates to position & velocity.

Fitting tk make this a heading

Given a Potential and a set of Actions (3 values J_R, J_z, J_phi), AutoFit
can fit the corresponding Torus. Once this is done the whole orbit is known.

The normal fitting procedure is:

1) Member function SetTopPot chooses plausible values for the
parameters of ToyIsochrone, and sets up original parameters of GenFnc
and AngMap (an object of class GenPar) via the function MakeGeneric (a
member function of GenPar).

2)

 */
class Torus : public PhaseSpaceMap {
private:
  Actions	J;                            // Torus actions
  double	E, Fs, Rmin, Rmax, zmax;      // Energy, df_value
					      // approx limits in R, z
  double        scale;			      // 1 unless has been scaled for interpolation
  Frequencies   Om;                           // orbital frequencies
  Errors        dc;                           // errors (see above)
  PoiTra       *PT;                           // Point transform
  ToyMap       *TM;                           // Toy Map (e.g. isochrone)
  GenFnc	GF;                           // Generating function (J,thT->JT)
  AngMap	AM;                           // Angle Mapping (th->thT)
  double Jtr[4];

  Angles      mirror_Angles(Angles,double) const;
  // For given angles & coord phi, find angles giving same x, -vR, -vz, vphi
  Velocity    mirror_Velocity(Velocity) const;      // return -vR, -vz, vphi
  GCY         build_GCY(Position,Velocity) const;   // return w=x,v

  void        SetMaps(const double*, const vec4&,
		      const GenPar&, const AngPar&);
  // set maps within Torus from parameters of point transform, toymap,
  // generating function S(J), dS/dJ
  void        SetMaps(const vec4&, const GenPar&, const AngPar&);
  // As above, but no point transform
  void        DelMaps();
  // Delete maps in Torus (avoiding memory leaks)
  void 	LevCof (const PSPD&, const Position&, const double&,
		const double&, PSPD&,
		double&, Vector<double,2>&, Matrix<double,2,2>&,
		Matrix<double,2,2>&) const;
  void 	LevCof (const PSPD&, const PSPD&, const double&,
		const double&, const double&, PSPD&,
		double&, Vector<double,2>&, Matrix<double,2,2>&,
		Matrix<double,4,2>&) const;
  void 	LevCof (const PSPD&, const PSPD&, const Vector<double,4>&,
		PSPD&, double&, Vector<double,2>&, Matrix<double,2,2>&,
		Matrix<double,4,2>&) const;
  // Compute chi, where
  //		chi^2 = a^2 * (R0 - R[th])^2 + b^2 * (z0 - z[th])^2
  //                      + vscale^2 * ( (vR0 - vR[th])^2 + (vz0 - vz[th])^2 )
  // and its derivatives w.r.t. th
  // for use when doing a Levenberg-Marquand minimisation
  void        SOSroot(const double, double&, double&) 	  const;
  // Used when finding surface of section at z=0. Determines z & dz/dth_z.
  double t1p;
  void        resSOSroot(const double, double&, double&)   const;
  // Used when finding SoS for orbit trapped by 1:1 resonance
  double RforSOS;
  void        SOS_z_root(const double, double&, double&) 	  const;
  // Used when finding surface of section at R=RSOS. Determines z & dz/dth_z.
  double GM(double,double,double);

/*     int          ThinFit          (Potential*,		  // galactic potential  */
/* 				    const int     =0,     // Full(0)/Half(1)Fit */
/* 				    const double  =0.001, // goal for |dJ|/|J| */
/* 				    const int     =700,   // max. number of Sn */
/* 				    const int     =200,	  // max. iterations */
/* 				    const int     =14,    // max. SN tailorings  */
/* 				    const int     =24,    // min. # of theta */
/* 							  //    per dim   */
/* 				    const int     =0, 	  // error output?  */
/* 				    const int     =3,	  // overdetermination */
/* 				    const int     =24,	  // min. # of cells */
/* 							  //	for angle fit */
/* 				    const int     =200);  // max. # of steps */
/* 							  //    on av. per cell */
/* // performs a fit either fully or partly (depending on parameter "type") */
/* // specialized for case of Jz << Jr */

/*     int          PointFit	   (Potential*,		  // galactic potential */
/* 				    const double  =0.001, // goal for |dJ|/|J| */
/* 				    const int     =200,   // max. number of Sn */
/* 				    const int     =200,	  // max. iterations */
/* 				    const int     =5,     // max. SN tailorings  */
/* 				    const int     =3,	  // overdetermination */
/* 				    const int     =24,	  // min. # of cells */
/* 							  //	for angle fit */
/* 				    const int     =200,   // max. # of steps */
/* 							  //    on av. per cell */
/* 				    const int     =24,    // min. # of theta */
/* 							  //    per dim  */
/* 				    const int     =0);	  // error output? */
/* // performs a complete fit: Actions are fixed while the torus parameters are */
/* // used as initial guess and are changed in order to fit the torus in the given */
/* // potential. */


public:
	Torus();
	Torus(double);
    Torus(const Torus&);
    Torus&	operator= (const Torus&);
    Torus&	operator*= (const double&);
    Torus&	operator+= (const Torus&);
    Torus&	operator-= (const Torus&);
    const Torus operator* (const double&);
    const Torus	 operator+ (const Torus&);
    const Torus	 operator- (const Torus&);
  // Constructors
   ~Torus() { DelMaps(); }
  // Destructor

   bool write_ebf(const string,
		  const string,             // write torus details to ebf file
      const string mode = "a"); // default mode: append);

   bool read_ebf (const string,
		  const string); // read torus details from ebf file


   double       energy  () const { return E; }
   double       fsample () const { return Fs; }
   Actions      actions () const { return J; }
   Frequencies  omega   () const { return Om; }
   Errors       errors  () const { return dc; }
   vec4         TP      () const { return TM->parameters(); }
   GenPar       SN      () const { return GF.parameters(); }
   AngPar       AP      () const { return AM.parameters(); }
   double       minR    () const { return Rmin; }
   double       maxR    () const { return Rmax; }
   double       maxz    () const { return zmax; }

   double       omega	 (const int i)const{ return Om(i); }
   double       error	 (const int i)const{ return dc(i); }
   double       action   (const int i)const{ return J(i); }
   double       TP	 (const int i)const{ return TM->parameters()(i); }
   double       SN	 (const int i)const{ return GF.coeff(i); }
   double       dS1	 (const int i)const{ return AM.dSdJ1(i);}
   double       dS2	 (const int i)const{ return AM.dSdJ2(i);}
   double       dS3	 (const int i)const{ return AM.dSdJ3(i);}
   int          n1       (const int i)const{ return GF.n1(i); }
   int          n2       (const int i)const{ return GF.n2(i); }
   // Functions returning parameters of the Torus

   void         show    (ostream&) const;  // Gives human readable torus details

   PoiTra&	 canmap  ()		    { return *PT; }
   ToyMap&	 toymap  ()		    { return *TM; }
   // returns addresses of point transform, toy map
   void	 SetPP	 (Potential*, const Actions);
   // Set point transform from potential, Jphi, Jz.
   void	 SetPP	 (Cheby, Cheby, Cheby, double);
   void	 SetPP	 (double*);
   void	 SetPP	 ();
   void	 SetTP	 (const vec4& );
   void	 SetSN	 (const GenPar& sn) { GF.set_parameters(sn); }
   void	 SetAP	 (const AngPar& ap) { AM.set_parameters(ap); }

   void         SetActions	   (const Actions&, const double=0.);
   void         SetFs   	   (const double);
   void         SetFrequencies    (const Frequencies&);
   void         TailorAndCutSN    (const double, const double, const double,
				   const int);
  /*   int          get		   (istream&); */
/*     void         put		   (ostream&)    const; */
    int          NumberofParameters()            const;
    int          NumberofSn        ()            const;

    PSPD         Forward	   (const PSPD&) const;
    PSPT         Forward3D	   (const PSPT&) const;
// these routines are here because they must be. They translate
// (Jr,Jl[,Jp],Tr,Tt,[Tp]) into (R,z[,phi],vR,vz[,vphi]), but give a
// WARNING if the actions do not match those of the Torus

    PSPD         Map	           (const Angles&) const;
// returns (R,z, vR,vz) given (Tr, Tt, phi)

    PSPD         MapfromToy	   (const Angles&) const;
// returns (R,z, vR,vz) given toy angles (tr, tt, phi)

    PSPD         MapfromToy	   (const Angles&, PSPD&, double& ) const;
// as above, but also returns  - (Jr,Jz,Tr,Tt) and
//                             - the determinant |d(Tr,Tt)/d(tr,tt)|

    PSPT         Map3D	           (const Angles&) const;
// returns (R,z,phi, vR,vz,vphi) given (Tr, Tt, Tphi)
    PSPT         MapfromToy3D	           (const Angles&) const;
// returns (R,z,phi, vR,vz,vphi) given (Tr, Tt, Tphi)

    double       DToverDt	   (const Angles&) const;
// given toy angles (tr,tt, phi) returns the determinant |d(Tr,Tt)/d(tr,tt)|
    double       DToverDt	   (const Angles&,Angles&) const;
// given toy angles (tr,tt, phi) returns the determinant
// |d(Tr,Tt)/d(tr,tt)| and the true angles

    double       DToverDt	   (const Angles&, PSPD&) const;
// as above, but also returns (Jr,Jz,Tr,Tt) in the second argument

    GCY		 FullMap	   (const Angles&) const;
// returns (R,z,phi, vR,vz,vphi) given (Tr, Tt, Tphi)

    Position     PosMap		   (const Angles&) const;
// returns just (R,z,phi) given (Tr, Tt, phi)

    void         SOS		   (ostream&, const int=200);
// writes positions of a R-pR surface of section computed from the orbit
// by default 200 angles are used.
    void resSOS			   (ostream&, const double);
//writes a point in the R-pR sos for orbits trappped by 1:1 resonance
    void         SOS_z		   (ostream&, const double,const int=200);
// writes positions of a z-pz surface of section computed from the orbit
// by default 200 angles are used.
    void         RpSOS		   (ostream&,const int=200);
// writes positions of a R-pR at z=0 phi=0 surface of section computed from the orbit
// by default 200 angles are used.
    int          SOS_z		   (double*, double*,
				    const double,const int=200);
// writes positions of a z-pz surface of section computed from the orbit
// by default 200 angles are used.

    double       DistancetoPoint   (const Position&) const;
// Returns distance from given position to nearest point on torus
    double       DistancetoPoint   (const Position&, double&, double&) const;
    double       DistancetoPoint_Ang   (const Position&, Angles&) const;
// Returns distance from given position to nearest point on torus
// Input toy angles for first guess, output for closest point.
    double       DistancetoPoint   (const Position&, Position&) const;
// Returns distance from given position to nearest point on torus (also given)
    Vector<double,2> DistancetoPSP (const PSPD&, double&, Angles&) const;
// Returns distance from given pspd to nearest point on torus

    Vector<double,4> DistancetoPSP (const PSPD&, Vector<double,4>&,
				    Angles&) const;
// Returns distance from given pspd to nearest point on torus
    void LevCof3D(const Angles&,Matrix<double,3,3>&, const Position&, PSPT&,double&,
		  Vector<double,3>&, Matrix<double,3,3>&, Matrix<double,3,3>&) const;
    void get_derivs(const Angles&, Matrix<double,3,3>&, Matrix<double,3,3>&) const;
//compute
//   chi^2 = (R0-R(th))^2+(z0-z(th))^2+R0^2(phi0-phi(th))^2 and its
//   derivatives wrs (toy) theta
    void LevCof3DTrue(const Angles&,Matrix<double,3,3>&, const Position&, PSPT&,double&,
		  Vector<double,3>&, Matrix<double,3,3>&) const;
//compute
//   chi^2 = (R0-R(th))^2+(z0-z(th))^2+R0^2(phi0-phi(th))^2 and its
//   derivatives wrs (true) theta

    void CheckLevCof(PSPD, Angles);
    void CheckLevCof(Position, Angles);

    double       DistancetoRadius  (const double) const;
// Returns distance from given axisymmetric Radius to nearest point on torus

    bool          containsPoint	   (const Position&) const;
// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise.

    bool containsPoint (const Position&,Velocity&,Matrix<double,2,2>&,Angles&,
		       Velocity&,Matrix<double,2,2>&,Angles&,bool,bool,bool,
		       double)  const;
/* Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the
 torus passes through the point given, this happens four times, in each case
 with a different velocity. However, only two of these are independent, since
 changing the sign of both vR and vz simultaneously gives the same orbit,
 at (2Pi-TR,Pi-Tz). For each of these two the two possible (vR,vz) and the
 matrix d(R,z)/d(Tr,Tz) is returned. If the first bool is true, the angles
 (TR,Tz,Tphi) are also calculated in each case. If the second bool is true,
 the TOY angles are returned.

   This is the function that does the actual work, but the end user will
 (probably) only use the related functions which use this one (see below).

   The other function can take as parameters (and thus return) any sensible
 mixture of 'Velocity's, 'GCY's (R,z,phi,vR,vz,vphi), Angles, and either
 the matrix d(R,z)/d(Tr,Tz) or the determinant |d(x,y,z)/d(Tr,Tz,Tphi)|
 which vanishes on the edge of the orbit, such that its inverse, the density
 of the orbit, diverges there (that's the reason why the density itself is
 not returned).

   Note that because Velocity and Angles are both just 3-vectors I have to call
 the function containsPoint_Ang if it only returns Angles, and
 containsPoint_Vel if it only returns Velocities
*/


// Those other functions
    bool containsPoint (const Position&,Velocity&,Velocity&,Matrix<double,2,2>&,
		       Angles&,Angles&, Velocity&,Velocity&,Matrix<double,2,2>&,
		       Angles&,Angles&,bool=false,bool=false,double=0.) const;
    bool containsPoint (const Position&,Velocity&,Matrix<double,2,2>&,Angles&,
		       Velocity&,Matrix<double,2,2>&,Angles&,bool=false,
		       bool=false,double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,Matrix<double,2,2>&,
		       Angles&,Angles&,GCY&,GCY&,Matrix<double,2,2>&,
		       Angles&,Angles&,bool=false,bool=false,double=0.) const;
    bool containsPoint (const Position&,GCY&,Matrix<double,2,2>&,Angles&,
		       GCY&, Matrix<double,2,2>&,Angles&,bool=false,
		       bool=false,double=0.) const;
    bool containsPoint(const Position&,Velocity&,Velocity&,double&,//used!
		      Angles&,Angles&,Velocity&,Velocity&,double&,
		      Angles&,Angles&,bool=false,bool=false,double=0.)const;
    bool containsPoint(const Position&,Velocity&,double&,Angles&,
		      Velocity&,double&,Angles&,bool=false,bool=false,
		      double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,double&,Angles&,Angles&,
		       GCY&,GCY&,double&,Angles&,Angles&,bool=false,
		       bool=false,double=0.) const;
    bool containsPoint (const Position&,GCY&,double&,Angles&,
		       GCY&,double&,Angles&,bool=false,bool=false,
		       double=0.) const;
// Lots of them
    bool containsPoint(const Position&,Velocity&,Velocity&,Angles&,Angles&,
		      Velocity&,Velocity&,Angles&,Angles&,bool=false,
		      bool=false,double=0.) const;
    bool containsPoint(const Position&,Velocity&,Angles&,Velocity&,Angles&,
		      bool=false,bool=false,double=0.)  const;
    bool containsPoint (const Position&,GCY&,GCY&,Angles&,Angles&,
		       GCY&,GCY&,Angles&,Angles&,bool=false,bool=false,
		       double=0.) const;
    bool containsPoint (const Position&,GCY&,Angles&,
		       GCY&,Angles&,bool=false,bool=false,double=0.) const;
    bool containsPoint_Ang (const Position&,Angles&,Angles&,bool=false,
			   bool=false,double=0.) const;
    bool containsPoint_Ang (const Position&,Angles&,Angles&,
			   Angles&,Angles&,bool=false,bool=false,
			   double=0.) const;

    bool containsPoint (const Position&,Velocity&,Velocity&,double&,
			Velocity&,Velocity&,double&,double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,double&,
		       GCY&,GCY&,double&,double=0.) const;
    bool containsPoint (const Position&,Velocity&,double&,
		       Velocity&,double&,double=0.) const;
    bool containsPoint (const Position&,GCY&,double&,GCY&,double&,
		       double=0.) const;
 // Lots and lots and lots.
    bool containsPoint (const Position&,Velocity&,Velocity&,Matrix<double,2,2>&,
		       Velocity&,Velocity&,Matrix<double,2,2>&,double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,Matrix<double,2,2>&,
		       GCY&,GCY&,Matrix<double,2,2>&,double=0.) const;
    bool containsPoint (const Position&,Velocity&,Matrix<double,2,2>&,
		       Velocity&,Matrix<double,2,2>&,double=0.) const;
    bool containsPoint (const Position&,GCY&,Matrix<double,2,2>&,
		       GCY&,Matrix<double,2,2>&,double=0.) const;

    bool containsPoint_Vel (const Position&,Velocity&,Velocity&,
			   Velocity&,Velocity&,double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,GCY&,GCY&,double=0.) const;
    bool containsPoint_Vel (const Position&,Velocity&,Velocity&,double=0.) const;
    bool containsPoint (const Position&,GCY&,GCY&,double=0.) const;

    PSPD	 StartPoint	   (const double=Pi, const double=Pih)const;
// returns (R,z,pR,pz) from given toy-angles. So the points cannot be used
// as representative distribution but just as start point for orbit integration

 /*    int        ManualFit	   (Potential*,	  // galactic potential */
/* 			    const double  =0.001, // goal for |dJ|/|J| */
/* 			    const int     =700,   // max. number of Sn */
/* 			    const int     =200,	  // max. iterations */
/* 			    const int     =14,    // max. SN tailorings  */
/* 			    const int     =3,	  // overdetermination */
/* 			    const int     =24,	  // min. # of cells */
/* 						  //	for angle fit */
/* 			    const int     =200,   // max. # of steps */
/* 						  //    on av. per cell */
/* 			    const int     =24,    // min. # of theta */
/* 						  //    per dim  */
/* 			    const int     =0);	  // error output? */
/* // performs a complete fit: Actions are fixed while the torus
parameters are */
/* // used as initial guess and are changed in order to fit the torus
in the given */
/* // potential. */

    int          AutoFit  	   (Actions,              // Actions
				    Potential*,		  // galactic potential
				    const double  =0.003, // goal for |dJ|/|J|
				    const int     =600,   // max. number of Sn
				    const int     =150,	  // max. iterations
				    const int     =12,    // max. SN tailorings
				    const int     =3,	  // overdetermination
				    const int     =16,	  // min. # of cells
							  //	for angle fit
				    const int     =200,   // max. # of steps
							  //    on av. per cell
				    const int     =12,    // min. # of theta
							  //    per dim
				    const int     =0);	  // error output?
// performs a complete fit: Actions are fixed while the torus parameters are
// used as initial guess and are changed in order to fit the torus in the given
// potential.

    int		FixMapsFit	(Actions,
				 Potential*,
				 const double	=0.003,
				 const int	=600,
				 const int	=150,
				 const int	=12,
				 const int	=3,
				 const int	=16,
				 const int	=200,
				 const int	=12,
				 const int	=0);
// Performs fit using a non-trivial point transformation

    int          AutoFit_TPfirst  	   (Actions,              // Actions
				    Potential*,		  // galactic potential
				    const double  =0.003, // goal for |dJ|/|J|
				    const int     =600,   // max. number of Sn
				    const int     =150,	  // max. iterations
				    const int     =12,    // max. SN tailorings
				    const int     =3,	  // overdetermination
				    const int     =16,	  // min. # of cells
							  //	for angle fit
				    const int     =200,   // max. # of steps
							  //    on av. per cell
				    const int     =12,    // min. # of theta
							  //    per dim
				    const int     =0);	  // error output?
// performs a complete fit: Actions are fixed while the torus parameters are
// used as initial guess and are changed in order to fit the torus in the given
// potential.
//
// In this case we use the strategy of fitting the toy potential parameters first
// then only varying the generating function


    int          AutoFit_SNfirst  	   (Actions,              // Actions
				    Potential*,		  // galactic potential
				    const double  =0.003, // goal for |dJ|/|J|
				    const int     =600,   // max. number of Sn
				    const int     =150,	  // max. iterations
				    const int     =12,    // max. SN tailorings
				    const int     =3,	  // overdetermination
				    const int     =16,	  // min. # of cells
							  //	for angle fit
				    const int     =200,   // max. # of steps
							  //    on av. per cell
				    const int     =12,    // min. # of theta
							  //    per dim
				    const int     =0);	  // error output?
// performs a complete fit: Actions are fixed while the torus parameters are
// used as initial guess and are changed in order to fit the torus in the given
// potential.
//
// In this case we use the strategy of fitting the low order terms in the generating
// function first then varying both the toy potential and the generating function
// (including higher order terms)


    int          AutoFit_LowJzonly (Actions,              // Actions
				    Potential*,		  // galactic potential
				    const double  =0.003, // goal for |dJ|/|J|
				    const int     =600,   // max. number of Sn
				    const int     =150,	  // max. iterations
				    const int     =12,    // max. SN tailorings
				    const int     =3,	  // overdetermination
				    const int     =16,	  // min. # of cells
							  //	for angle fit
				    const int     =200,   // max. # of steps
							  //    on av. per cell
				    const int     =12,    // min. # of theta
							  //    per dim
				    const int     =0);	  // error output?
// performs a complete fit: Actions are fixed while the torus parameters are
// used as initial guess and are changed in order to fit the torus in the given
// potential.
//
// In this case we use the strategy of fitting the low order terms in the generating
// function first then varying both the toy potential and the generating function
// (including higher order terms)


    int      FitWithFixToyPot  	   (Actions,              // Actions
				    vec4,               // Toy parameters
				    Potential*,		  // galactic potential
				    const double  =0.001, // goal for |dJ|/|J|
				    const int     =700,   // max. number of Sn
				    const int     =200,	  // max. iterations
				    const int     =14,    // max. SN tailorings
				    const int     =3,	  // overdetermination
				    const int     =24,	  // min. # of cells
							  //	for angle fit
				    const int     =200,   // max. # of steps
							  //    on av. per cell
				    const int     =24,    // min. # of theta
							  //    per dim
				    const int     =0, 	  // error output?
				    const bool    =false);//spherical potential?
// performs a complete fit: Actions and toy parameters are fixed while the
// other torus parameters are used as initial guess and are changed in order
// to fit the torus in the given potential.
//
// Most useful when looking to interpolate between tori.




 /*    int          ManualHalfFit   (Potential*,  // galactic potential */
/* 			    const double  =0.001, // goal for |dJ|/|J| */
/* 			    const int     =700,   // max. number of Sn */
/* 			    const int     =200,	  // max. iterations */
/* 			    const int     =14,    // max. SN tailorings  */
/* 			    const int     =24,    // min. # of theta */
/* 						  //    per dim  */
/* 			    const int     =0);	  // error output? */
/* // performs a fit partly: as above but angle map is not fitted. */


    double     check              (Potential*);           // dH/H

    void  SetToyPot               (Potential*,            // galactic potential
				   const Actions,         // Actions to fit
				   const double=0 );      // Scale radius

    void AutoPTTorus             (Potential*,             // galactic potential
				  const Actions,          // Actions to fit
				  const IsoPar );	  // TM parameters

    void         FindLimits        ();
};

// **************************** inline functions **************************** //

inline Torus::Torus() : J(0.), E(0.), Fs(0.), Om(0.), dc(0.), PT(0), TM(0), scale(1) {}

inline Torus::Torus(double x) : J(0.), E(0.), Fs(0.), Om(0.), dc(0.), PT(0), TM(0), scale(1) {}

inline Torus::Torus(const Torus& T)
    : PhaseSpaceMap(), J(T.J), E(T.E), Fs(T.Fs), Om(T.Om), dc(T.dc), scale(T.scale)
{
  if((T.PT)->NumberofParameters()) {
    double tmp[(T.PT)->NumberofParameters()];
    (T.PT)->parameters(tmp);
    SetMaps(tmp,
	    (T.TM)->parameters(),
	    (T.GF).parameters(),
	    (T.AM).parameters());
  }  else
    SetMaps((T.TM)->parameters(),
	    (T.GF).parameters(),
	    (T.AM).parameters());
}

inline void Torus::SetActions(const Actions& j, const double fs)
{ J  = j;  Fs = fs; }
inline void Torus::SetFrequencies(const Frequencies& om)
{  Om = om; }
inline void Torus::SetFs(const double fs)
{  Fs = fs; }


inline int Torus::NumberofParameters() const
{
    return   PT->NumberofParameters()
	   + TM->NumberofParameters()
	   + GF.NumberofParameters()
	   + AM.NumberofParameters();
}

inline int Torus::NumberofSn() const
{
    return GF.NumberofParameters();
}

inline PSPD Torus::Forward(const PSPD &JT) const
{
    if((JT(0) != J(0)) || JT(1) != J(1))
	cerr<<" WARNING: Torus::Forward() called with different action(s)\n";
    return JT >> AM >> GF >> (*TM) >> (*PT);
}
inline PSPT Torus::Forward3D(const PSPT &JT) const
{
    if((JT(0) != J(0)) || JT(1) != J(1) || JT(2) != J(2))
	cerr<<" WARNING: Torus::Forward() called with different action(s)\n";
    return JT >> AM >> GF >> (*TM) >> (*PT);
}

inline PSPD Torus::Map(const Angles& A) const {
    return PSPD(J(0),J(1),A(0),A(1)) >> AM >> GF >> (*TM) >> (*PT);
}
inline PSPT Torus::Map3D(const Angles& A) const {
  return PSPT(J(0),J(1),J(2),A(0),A(1),A(2)) >> AM >> GF >> (*TM) >> (*PT);
}

inline PSPD Torus::MapfromToy(		// return:	(R,z,vR,vz)
	    const Angles& A) const	// input:       (tr,tt,phi)
{
    return    PSPD(J(0),J(1),A(0),A(1)) >> GF >> (*TM) >> (*PT);
}
inline PSPT Torus::MapfromToy3D(const Angles& A) const
{
  return PSPT(J(0),J(1),J(2),A(0),A(1),A(2))>> GF >> (*TM) >> (*PT);
}

inline PSPD Torus::MapfromToy(		// return:	(R,z,vR,vz)
	    const Angles& A,		// input:       (tr,tt,phi)
	    PSPD&	  JT,		// output:	(Jr,Jz,Tr,Tt)
	    double&       Det) const    // output:      |d(Tr,Tt)/d(tr,tt)|
{
    double dTdt[2][2];
    PSPD Jt = PSPD(J(0),J(1),A(0),A(1));
    JT      = AM.BackwardWithDerivs(Jt,dTdt);
    Det     = dTdt[0][0]*dTdt[1][1] - dTdt[0][1]*dTdt[1][0];
    return    Jt >> GF >> (*TM) >> (*PT);
}

inline double Torus::DToverDt(				//return |d(Tr,Tt)/d(tr,tt)|
			      const Angles& A,		// input:	(tr,tt,phi)
			      Angles &AT) const		//output: (Tr,Tt,Tp)
{
	double dTdt[2][2];
	PSPT JT = AM.Backward3DWithDerivs(PSPT(J(0),J(1),J(2),A(0),A(1),A(2)), dTdt );
	for(int i=0;i<3;i++) AT[i]=JT[i+3];
	return  dTdt[0][0]*dTdt[1][1] - dTdt[0][1]*dTdt[1][0];
}

inline double Torus::DToverDt(		// return: 	|d(Tr,Tt)/d(tr,tt)|
			      const Angles& A,		// input:	(tr,tt,phi)
			      PSPD& JT) const		// output:	(Jr,Jz,Tr,Tt)
{
	double dTdt[2][2];
	JT = AM.BackwardWithDerivs( PSPD(J(0),J(1),A(0),A(1)) , dTdt );
	return  dTdt[0][0]*dTdt[1][1] - dTdt[0][1]*dTdt[1][0];
}

inline double Torus::DToverDt(		// return: 	|d(Tr,Tt)/d(tr,tt)|
	      const Angles& A) const	// input:	(tr,tt,phi)
{
    double dTdt[2][2];
    AM.BackwardWithDerivs( PSPD(J(0),J(1),A(0),A(1)) , dTdt );
    return  dTdt[0][0]*dTdt[1][1] - dTdt[0][1]*dTdt[1][0];
}

inline GCY Torus::FullMap(const Angles& A) const
{
    GCY  gcy;
    PSPT qp = PSPT(J(0),J(1),J(2),A(0),A(1),A(2)) >> AM >> GF >> (*TM) >> (*PT);
    gcy[0] = qp(0);		// R
    gcy[1] = qp(1);		// z
    gcy[2] = qp(2);		// phi
    gcy[3] = qp(3);		// vR   = pR
    gcy[4] = qp(4);		// vz   = pz
    gcy[5] = J(2)/qp(0);	// vphi = Lz/R
    return gcy;
}

inline Position Torus::PosMap(const Angles& A) const
{
    Position X;
    PSPT qp = PSPT(J(0),J(1),J(2),A(0),A(1),A(2)) >> AM >> GF >> (*TM) >> (*PT);
    X[0] = qp(0);
    X[1] = qp(1);
    X[2] = qp(2);
    return X;
}

inline Angles Torus::mirror_Angles(Angles A, double phi) const {
  Angles out;
  out[0] = TPi-A[0];
  out[1] = Pi-A[1];
  out[2] = 2*phi-A[2];
  for(int i=0;i!=3;i++) {
    while(out[i]>TPi) out[i] -= TPi;
    while(out[i]<0)   out[i] += TPi;
  }
  return out;
}
inline Velocity Torus::mirror_Velocity(Velocity v) const {
  Velocity out=-v;
  out[2] *= -1;
  return out;
}

inline GCY Torus::build_GCY(Position x, Velocity v) const {
  GCY out;
  for(int i=0;i!=3;i++) {
    out[i] = x[i]; out[i+3] = v[i];
  }
  return out;
}

//////// Very simple containsPoint.
inline bool Torus::containsPoint(    // return:     is X ever hit ?
	   const Position& X)	    // input:      position X
	   const
{
    return (DistancetoPoint(X) == 0.)?  true : false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Various ways of calling containsPoints if you're only interested in
//  certain properties of those points.
//

inline bool Torus::containsPoint (const Position& X,Velocity& v11,
			  Matrix<double,2,2>& M_1, Angles& A11,
			  Velocity& v21, Matrix<double,2,2>& M_2,
			  Angles& A21, bool toy, bool usea, double delr) const
{
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,true,toy,usea,delr)) return false;
  return true;
}

inline bool Torus::containsPoint (const Position& X,Velocity& v11,Velocity& v12,
			  Matrix<double,2,2>& M_1, Angles& A11, Angles& A12,
			  Velocity& v21, Velocity& v22, Matrix<double,2,2>& M_2,
			  Angles& A21,Angles& A22, bool toy, bool usea,
			  double delr) const
{
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,true,toy,usea,delr)) return false;
  v12=mirror_Velocity(v11);      v22=mirror_Velocity(v21);
  A12=mirror_Angles(A11,X(2));   A22=mirror_Angles(A21,X(2));
  return true;
}
inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
			  Matrix<double,2,2>& M_1, Angles& A11, Angles& A12,
			  GCY& w21, GCY& w22, Matrix<double,2,2>& M_2,
			  Angles& A21,Angles& A22, bool toy, bool usea,
			  double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint(X,v11,v12,M_1,A11,A12,v21,v22,M_2,A21,A22,toy,usea,delr))
    return false;
  w11 = build_GCY(X,v11); w12 = build_GCY(X,v12);
  w21 = build_GCY(X,v21); w22 = build_GCY(X,v22);
  return true;
}
inline bool Torus::containsPoint (const Position& X,GCY& w11,
			  Matrix<double,2,2>& M_1, Angles& A11,
			  GCY& w21, Matrix<double,2,2>& M_2,
			  Angles& A21, bool toy, bool usea, double delr) const
{
  Velocity v11,v21;
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,true,toy,usea,delr)) return false;
  w11 = build_GCY(X,v11);  w21 = build_GCY(X,v21);
  return true;
}
///////

inline bool Torus::containsPoint (const Position& X,Velocity& v11,Velocity& v12,//used
			  double& d_1, Angles& A11, Angles& A12,
			  Velocity& v21, Velocity& v22, double& d_2,
			  Angles& A21,Angles& A22, bool toy, bool usea,
				 double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,v12,M_1,A11,A12,v21,v22,M_2,A21,A22,toy,usea,delr))
    return false;
  d_1  = X(0) * (M_1(0,0)*M_1(1,1) - M_1(0,1)*M_1(1,0));
  d_2  = X(0) * (M_2(0,0)*M_2(1,1) - M_2(0,1)*M_2(1,0));
  return true;
}
inline bool Torus::containsPoint (const Position& X,Velocity& v11,
			  double& d_1, Angles& A11,
			  Velocity& v21, double& d_2,
			  Angles& A21, bool toy, bool usea, double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,toy,usea,delr)) return false;
  d_1  = X(0) * (M_1(0,0)*M_1(1,1) - M_1(0,1)*M_1(1,0));
  d_2  = X(0) * (M_2(0,0)*M_2(1,1) - M_2(0,1)*M_2(1,0));
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
			  double& d_1, Angles& A11, Angles& A12,
			  GCY& w21, GCY& w22, double& d_2,
			  Angles& A21,Angles& A22, bool toy, bool usea,
				 double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint(X,v11,v12,d_1,A11,A12,v21,v22,d_2,A21,A22,toy,usea,delr))
    return false;
  w11 = build_GCY(X,v11); w12 = build_GCY(X,v12);
  w21 = build_GCY(X,v21); w22 = build_GCY(X,v22);
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11, double& d_1,
				 Angles& A11, GCY& w21, double& d_2,
				 Angles& A21, bool toy, bool usea,
				 double delr) const
{
  Velocity v11,v21;
  if(!containsPoint(X,v11,d_1,A11,v21,d_2,A21,toy,usea,delr)) return false;
  w11 = build_GCY(X,v11);  w21 = build_GCY(X,v21);
  return true;
}
///////

inline bool Torus::containsPoint (const Position& X,Velocity& v11,Velocity& v12,
			  Angles& A11, Angles& A12,Velocity& v21,Velocity& v22,
			  Angles& A21,Angles& A22, bool toy, bool usea,
				 double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,v12,M_1,A11,A12,v21,v22,M_2,A21,A22,toy,usea,delr)) return false;
  return true;
}
inline bool Torus::containsPoint (const Position& X,Velocity& v11, Angles& A11,
			  Velocity& v21, Angles& A21, bool toy, bool usea,
				 double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,toy,usea,delr)) return false;
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
			  Angles& A11, Angles& A12, GCY& w21, GCY& w22,
			  Angles& A21,Angles& A22, bool toy, bool usea,
				 double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,w11,w12,M_1,A11,A12,w21,w22,M_2,A21,A22,toy,usea,delr))
    return false;
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11, Angles& A11,
			  GCY& w21, Angles& A21, bool toy, bool usea,
				 double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,w11,M_1,A11,w21,M_2,A21,toy,usea,delr)) return false;
  return true;
}

inline bool Torus::containsPoint_Ang (const Position& X, Angles& A11,
			  Angles& A21, bool toy, bool usea, double delr) const
{
  Velocity v11,v21;
  if(!containsPoint(X,v11,A11,v21,A21,toy,usea,delr)) return false;
  return true;
}

inline bool Torus::containsPoint_Ang (const Position& X, Angles& A11,Angles& A12,
			  Angles& A21, Angles& A22, bool toy, bool usea,
				     double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint(X,v11,v12,A11,A12,v21,v22,A21,A22,toy,usea,delr)) return false;
  return true;
}

///////////////////     Now if you're not interested in the Angles

inline bool Torus::containsPoint (const Position& X,Velocity& v11,Velocity& v12,
				 Matrix<double,2,2>& M_1, Velocity& v21,
				 Velocity& v22, Matrix<double,2,2>& M_2,
				 double delr) const
{
  if(!containsPoint(X,v11,M_1,v21,M_2,delr)) return false;
  v12=mirror_Velocity(v11);      v22=mirror_Velocity(v21);
  return true;
}
inline bool Torus::containsPoint (const Position& X,Velocity& v11,
				 Matrix<double,2,2>& M_1,
				 Velocity& v21, Matrix<double,2,2>& M_2,
				 double delr) const
{
  Angles A11, A21;
  if(!containsPoint(X,v11,M_1,A11,v21,M_2,A21,false,false,false,delr)) return false;
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
				 Matrix<double,2,2>& M_1,
				 GCY& w21, GCY& w22, Matrix<double,2,2>& M_2,
				 double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint(X,v11,v12,M_1,v21,v22,M_2,delr)) return false;
  w11 = build_GCY(X,v11); w12 = build_GCY(X,v12);
  w21 = build_GCY(X,v21); w22 = build_GCY(X,v22);
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,
			  Matrix<double,2,2>& M_1,
			  GCY& w21, Matrix<double,2,2>& M_2, double delr) const
{
  Velocity v11,v21;
  if(!containsPoint(X,v11,M_1,v21,M_2,delr)) return false;
  w11 = build_GCY(X,v11);  w21 = build_GCY(X,v21);
  return true;
}
///////

inline bool Torus::containsPoint (const Position& X,Velocity& v11,Velocity& v12,
				 double& d_1,
				 Velocity& v21, Velocity& v22, double& d_2,
				 double delr) const
{
  if(!containsPoint(X,v11,d_1,v21,d_2,delr)) return false;
  v12=mirror_Velocity(v11);      v22=mirror_Velocity(v21);
  return true;
}
inline bool Torus::containsPoint (const Position& X,Velocity& v11,double& d_1,
			  Velocity& v21, double& d_2, double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,M_1,v21,M_2,delr)) return false;
  d_1  = X(0) * (M_1(0,0)*M_1(1,1) - M_1(0,1)*M_1(1,0));
  d_2  = X(0) * (M_2(0,0)*M_2(1,1) - M_2(0,1)*M_2(1,0));
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
				 double& d_1, GCY& w21, GCY& w22,
				 double& d_2, double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint(X,v11,v12,d_1,v21,v22,d_2,delr)) return false;
  w11 = build_GCY(X,v11); w12 = build_GCY(X,v12);
  w21 = build_GCY(X,v21); w22 = build_GCY(X,v22);
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11, double& d_1,
			  GCY& w21, double& d_2, double delr) const
{
  Velocity v11,v21;
  if(!containsPoint(X,v11,d_1,v21,d_2,delr)) return false;
  w11 = build_GCY(X,v11);  w21 = build_GCY(X,v21);
  return true;
}


///////

inline bool Torus::containsPoint_Vel (const Position& X,Velocity& v11,
			       Velocity& v12,Velocity& v21,Velocity& v22,
				     double delr) const
{
  if(!containsPoint_Vel(X,v11,v21,delr)) return false;
  v12=mirror_Velocity(v11);      v22=mirror_Velocity(v21);
  return true;
}
inline bool Torus::containsPoint_Vel(const Position& X,Velocity& v11,
				    Velocity& v21, double delr) const
{
  Matrix <double,2,2> M_1,M_2;
  if(!containsPoint(X,v11,M_1,v21,M_2,delr)) return false;
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11,GCY& w12,
				 GCY& w21, GCY& w22, double delr) const
{
  Velocity v11,v12,v21,v22;
  if(!containsPoint_Vel(X,v11,v12,v21,v22,delr)) return false;
  w11 = build_GCY(X,v11); w12 = build_GCY(X,v12);
  w21 = build_GCY(X,v21); w22 = build_GCY(X,v22);
  return true;
}

inline bool Torus::containsPoint (const Position& X,GCY& w11, GCY& w21,
				 double delr) const
{
  Velocity v11,v21;
  if(!containsPoint_Vel(X,v11,v21,delr)) return false;
  w11 = build_GCY(X,v11);  w21 = build_GCY(X,v21);
  return true;
}
// AND we're done with containsPoint, finally.
////////////////////////////////////////////////////////////////////////////////

inline PSPD Torus::StartPoint(const double th1, const double th2) const
{
    return PSPD(J(0),J(1),th1,th2) >> GF >> (*TM) >> (*PT);
}

/* inline bool Torus::ManualFit(Potential *Phi, const double tol, const int Max, */
/* 			    const int Mit, const int Nta, const int Over,  */
/* 			    const int Ncl, const int ipc, const int Nth,  */
/* 			    const int err) */
/* { */
/*     register int F; */
/*     GenPar SN=GF.parameters(); */
/*     AngPar AP=AM.parameters(); */
/*     F = AllFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP, */
/* 	       Om,E,dc,0,false,Nta,ipc,E,Nth,err);  */
/*     if(F && J(1)<0.005*J(0)) { */
/*       Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc; */
/*       vec4 tmpIP = TM->parameters(); */
/*       GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F; */
/*       SN=GF.parameters(); */
/*       F = LowJzFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP, */
/* 		   Om,E,dc,0,Nta,ipc,E,Nth,err); */
/*       if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) || (F==-4 || F== -1))) {  */
/* 	Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;  */
/* 	TM->set_parameters(tmpIP); F=oF; */
/*       } */
/*     }  */
/*     if(F && J(0)<0.05*J(1)) { // possibly something with number of terms used */
/*       Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc; */
/*       vec4 tmpIP = TM->parameters(); */
/*       GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F; */
/*       SN=GF.parameters(); */

/*       AutoPTTorus(Phi,J,3.); */
/*       F = PTFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP, */
/* 		Om,E,dc,0,Nta,ipc,E,Nth,err); */
/*       if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) || */
/* 	       ((F==-4 && oF!=-1) || F== -1))) { */
/* 	Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP; */
/* 	TM->set_parameters(tmpIP); F=oF; SetPP(); */
/*       } */
/*     } */
/*     //SN.write(cerr); */
/*     GF.set_parameters(SN); */
/*     AM.set_parameters(AP); */
/*     FindLimits(); */
/*     return F; */
/* } */

inline int Torus::AutoFit(Actions Jin, Potential *Phi, const double tol,
			  const int Max, const int Mit, const int Nta,
			  const int Over, const int Ncl,
			  const int ipc, const int Nth, const int err)
{
	J = Jin;
	SetToyPot(Phi,Jin);
  // show(cerr);
	register int F;
	GenPar SN=GF.parameters();
	AngPar AP=AM.parameters();
	F = AllFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		   Om,E,dc,0,0,Nta,ipc,E,Nth,err);
	double JzFactor = 0.5 *
			  (1 + tanh( (J(0)-0.5*fabs(J(2))) / (0.25*fabs(J(2)) ) ) ) ;
	double LowJzFitParam = 0.005 + 0.05*JzFactor;

	if(F && (J(1)<LowJzFitParam*J(0)) ) {
		Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
		vec4 tmpIP = TM->parameters();
		GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
		SN=GF.parameters();
		F = LowJzFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
			     Om,E,dc,0,Nta,ipc,E,Nth,err);
		if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) || (F==-4 || F== -1))) {
			Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
			TM->set_parameters(tmpIP); F=oF;
		}
	}

	if(F && J(0)<0.05*J(1)) { // possibly something with number of terms used
		Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
		//SetToyPot(Phi,J,0);//By now TM screwed up; start over
		IsoPar tmpIP = TM->parameters();
		printf("IsoPar: %f %f %f %f\n",tmpIP[0],tmpIP[1],tmpIP[2],tmpIP[3]);
		GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
		SN=GF.parameters();
		if(err) cerr << "using AutoPTTorus\n";
		AutoPTTorus(Phi,J,tmpIP);
		if(err) cerr << "using PTFit\n";
		F = PTFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
			  Om,E,dc,0,Nta,ipc,E,Nth,err);
		/*if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) ||
			 ((F==-4 && oF!=-1) || F== -1))) {
			Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
			TM->set_parameters(tmpIP); F=oF; SetPP();
		}*/
	}
	GF.set_parameters(SN);
	AM.set_parameters(AP);

	FindLimits();
	return F;
}

inline int Torus::FixMapsFit(Actions Jin, Potential *Phi, const double tol,
			     const int Max, const int Mit, const int Nta,
			     const int Over, const int Ncl,
			     const int ipc, const int Nth, const int err)
{
	J = Jin;
	SetToyPot(Phi,Jin);
  // show(cerr);
	register int F;
	GenPar SN=GF.parameters();
	AngPar AP=AM.parameters();
	Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
	vec4 tmpIP = TM->parameters();
	GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
	SN=GF.parameters();
	if(err) cerr << "using AutoPTTorus\n";
	AutoPTTorus(Phi,J,tmpIP);
	if(err) cerr << "using PTFit\n";
	F = PTFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		  Om,E,dc,0,Nta,ipc,E,Nth,err);
	if(F && ((dc(0) > tmpdc(0)))) {
		Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
		TM->set_parameters(tmpIP); F=oF; SetPP();
	}
	GF.set_parameters(SN);
	AM.set_parameters(AP);

	FindLimits();
	return F;
}


inline int Torus::AutoFit_TPfirst(Actions Jin, Potential *Phi, const double tol,
				  const int Max, const int Mit, const int Nta,
				  const int Over, const int Ncl,
				  const int ipc, const int Nth, const int err)
{
  J = Jin;
  SetToyPot(Phi,Jin);
  // show(cerr);
  register int F;
  GenPar SN=GF.parameters();
  AngPar AP=AM.parameters();
  F = AllFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
	       Om,E,dc,0,1,Nta,ipc,E,Nth,err);
  double JzFactor = 0.5 *
    (1 + tanh( (J(0)-0.5*fabs(J(2))) / (0.25*fabs(J(2)) ) ) ) ;
  double LowJzFitParam = 0.005 + 0.05*JzFactor;

  if(F && (J(1)<LowJzFitParam*J(0)) ) {
    Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
    vec4 tmpIP = TM->parameters();
    GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
    SN=GF.parameters();
    F = LowJzFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		 Om,E,dc,0,Nta,ipc,E,Nth,err);
    if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) || (F==-4 || F== -1))) {
      Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
      TM->set_parameters(tmpIP); F=oF;
    }
  }

  if(F && J(0)<0.05*J(1)) { // possibly something with number of terms used
    Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
    vec4 tmpIP = TM->parameters();
    GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
    SN=GF.parameters();
    if(err) cerr << "using AutoPTTorus\n";
    AutoPTTorus(Phi,J,tmpIP);
    if(err) cerr << "using PTFit\n";
    F = PTFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
	      Om,E,dc,0,Nta,ipc,E,Nth,err);
    if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) ||
	     ((F==-4 && oF!=-1) || F== -1))) {
      Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
      TM->set_parameters(tmpIP); F=oF; SetPP();
    }
  }
    GF.set_parameters(SN);
    AM.set_parameters(AP);

    FindLimits();

    return F;
}


inline int Torus::AutoFit_SNfirst(Actions Jin, Potential *Phi, const double tol,
				  const int Max, const int Mit, const int Nta,
				  const int Over, const int Ncl,
				  const int ipc, const int Nth, const int err)
{
  J = Jin;
  SetToyPot(Phi,Jin);
  // show(cerr);
  register int F;
  GenPar SN=GF.parameters();
  AngPar AP=AM.parameters();


  F = AllFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
	     Om,E,dc,0,2,Nta,ipc,E,Nth,err);
  double JzFactor = 0.5 *
    (1 + tanh( (J(0)-0.5*fabs(J(2))) / (0.25*fabs(J(2)) ) ) ) ;
  double LowJzFitParam = 0.005 + 0.05*JzFactor;

   if(F && (J(1)<LowJzFitParam*J(0)) ) {

    Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
    vec4 tmpIP = TM->parameters();
    GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
    SN=GF.parameters();
    SetToyPot(Phi,Jin);
    F = LowJzFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		 Om,E,dc,0,Nta,ipc,E,Nth,err);
    if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) || (F==-4 || F== -1))) {
      Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
      TM->set_parameters(tmpIP); F=oF;
    }
  }

  if(F && J(0)<0.05*J(1)) { // possibly something with number of terms used
    Frequencies tmpOm = Om; double tmpE = E; Errors tmpdc = dc;
    IsoPar tmpIP = TM->parameters();
    GenPar tmpSN = SN; AngPar tmpAP = AP; int oF = F;
    SN=GF.parameters();
    if(err) cerr << "using AutoPTTorus\n";
    AutoPTTorus(Phi,J,tmpIP);
    if(err) cerr << "using PTFit\n";
    F = PTFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
	      Om,E,dc,0,Nta,ipc,E,Nth,err);
    if(F && ((dc(0) > tmpdc(0) && oF!=-4 && oF!=-1) ||
	     ((F==-4 && oF!=-1) || F== -1))) {
      Om = tmpOm; E = tmpE; dc = tmpdc; SN = tmpSN; AP = tmpAP;
      TM->set_parameters(tmpIP); F=oF; SetPP();
    }
  }
    GF.set_parameters(SN);
    AM.set_parameters(AP);

    FindLimits();

    return F;
}


inline int Torus::AutoFit_LowJzonly(Actions Jin, Potential *Phi, const double tol,
				  const int Max, const int Mit, const int Nta,
				  const int Over, const int Ncl,
				  const int ipc, const int Nth, const int err)
{
  J = Jin;
  SetToyPot(Phi,Jin);
  // show(cerr);
  register int F;
  GenPar SN=GF.parameters();
  AngPar AP=AM.parameters();

    F = LowJzFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		 Om,E,dc,0,Nta,ipc,E,Nth,err);

    GF.set_parameters(SN);
    AM.set_parameters(AP);

    FindLimits();

    return F;
}



inline int Torus::FitWithFixToyPot(Actions Jin, vec4 TPin,
			   Potential *Phi, const double tol,
			   const int Max, const int Mit, const int Nta,
			   const int Over, const int Ncl,
			   const int ipc, const int Nth, const int err,
			   const bool sph_pot)
{
  J = Jin;
  SetToyPot(Phi,Jin);
  TM->set_parameters(TPin);
  //  cerr << TP() << '\n';IP
  register int F;
  GenPar SN=GF.parameters();
  AngPar AP=AM.parameters();
  printf("FitWithFixToyPot: %f %f %f %f\n",TPin[0],TPin[1],TPin[2],TPin[3]);
  F = ToySetFit(J,Phi,tol,Max,Mit,Over,Ncl,*PT,*TM,SN,AP,
		Om,E,dc,Nta,ipc,E,Nth,err,false);

  GF.set_parameters(SN);
  AM.set_parameters(AP);

  FindLimits();

  return F;
}

inline double Torus::check(Potential *Phi) {

  double H, dH,lambda,tmp=0.;
  int nbad;
  GenPar SN=GF.parameters();
  if(SbyLevMar(J,Phi,7,std::max(24, 6*(SN.NumberofN1()/4+1)),
	       std::max(24, 6*(SN.NumberofN2())/4+1),1,
	       0.,0.,SN,*PT,*TM,lambda,H,dH,nbad,tmp,0)) {
    cerr << "Got some negative actions\n"; return 0.;
  }
  return dH;
}
Torus InterpTorus(Torus*,double*,int,double);
Torus InterpTorus2(Torus*,double*,int,double);
Torus InterpTorus_n2(Torus**,int,Actions,Actions,Actions);
Torus InterpTorus_2n(Torus**,int,Actions,Actions,Actions);
Torus InterpTorus(Torus***,Actions,Actions,Actions);

#include "resTorus.h"
#endif
