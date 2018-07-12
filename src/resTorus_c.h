/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
/* Class that provides tori of orbits trapped at corotation computed from
   * p-theory (Binney 2017) */
#ifndef _resTorus_c_
#define _resTorus_c_
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "Torus.h"
#include "Numerics.templates"
#include "Numerics.h"
#include "Random.h"
#include "jjb_utils.h"
#include "PJMebf.h"
#include "PJMCoords.h"
#include "PJM_cline.h"
#include "falPot.h"
#include "bar_pot.h"
#include "eTorus.h"

class resTorus_c{
	private:
		eTorus eT;
		Torus **Tgrid; //Pointer to grid of fitted un-trapped tori
		int np;//# of values of 1st index of Tgrid
		int inout;// +-1 if circilating, otherwise 0
		Actions Jbar_grid,dJ_grid;
		double lJ;// libration action
		double lOmega; //Frequency of libration
		double *t1grd,*tpgrd,*Dgrd;//pointers to stores of t1 and corresponding t_true
		double G;//dOmega1/dJ1
		vec8 hres;// resonant term and its 1st 2 derivatives
		int3 resN;// integers defining resonance
		Actions resJp; //Actions of underlying resonant torus
		double t1c,Dt1;// t1=resN[0]*t_r+resN[2]*t_p librates by +-Dt1 about t1c
		bool set_Dtheta1(void);// set Dt1
		void Dtheta1_fn(double,double&,double&) const;
		double theta_dot(double,double);// thetadot as f(t1,Jr-resJp[0])
		void store_theta(void);//Constructs table of t1 versus libration angle
		int nangle;// size of above table
		double from_librationAngle(double,double&);//Get t1 and Jr-resJp[0] from libration angle
		void getImin_max(void);// Set Imin,max
		void getDelta(double,double*);// find 2 values of Delta=Jr-resJp[0] & their difference
		double dDelta(double,double);// computes dJ/dtp
		double get_H(eTorus&);// extract mean H from eTorus
		double get_driver(eTorus&);// extract resonant Fourier term from eTorus
		double get_OLR(eTorus&);// extract OLR Fourier term from eTorus
		double librationAngle(double,double);//return libration angle from (t1,Jr-resJp[0]).
		void findJp(Actions&,Potential*,bar_pot*,double,double);//find Jphi of resonance given Jr
		bool fix_derivs(Actions&,Actions&,Actions&,Potential*,bar_pot*,double,double);
		Actions fixJ(Actions,const Angles&,double);//Add contributions of non-resonant terms
		Actions fixJ(Actions,const Angles&,double,double,Matrix<double,3,3>&);//Add contributions of non-resonant terms
		double resCond(Frequencies Om){
			return .5*(resN[0]*Om[0]+resN[1]*Om[1]+resN[2]*Om[2]);
		}
		bool oldV(Velocity*,Angles*,int,const Velocity&,const Angles& A,double&);
	public:
		resTorus_c(Torus**,int,Actions,Potential*,bar_pot*,double,double);
		~resTorus_c(){delete[] t1grd; delete[] tpgrd; delete[] Dgrd;};
		double I,Imin,Imax;
		void check_angle(void);
		int setI(double);
		int setI(double,int);
		void LevMstart(Angles,const Position,PSPT&,Vector<double,3>&,
			       Matrix<double,3,3>&,Matrix<double,3,3>&,double&,double&);
		void LevMstep(Angles&,const Position,PSPT&,Vector<double,3>&,
			      Matrix<double,3,3>&,Matrix<double,3,3>&,double&,double&,double&);
		double librationAction(void);// return resJ
		double librationOmega(void);// return libration frequency
		double from_librationAngle(double,double&,double&);//Get t1, dt1/dla and Jr-resJp[0] from libration angle
		GCY FullMap(const Angles&);// (R,..,v_phi) a f(theta)
		void SOS(ostream&,int);// write z=0 SoS to file
		int containsPoint(const Position&,Velocity*,Angles*,double*);//returns number of angles at which torus reaches given position
		void getJs(Actions&,Actions&,double);//returns extremes of unperturbed actions
		Actions get_resJp(void){return resJp;}
		Frequencies omega(void){
			Frequencies Om=eT.omega(); Om[2]=librationOmega(); return Om;
		}
};
#endif
