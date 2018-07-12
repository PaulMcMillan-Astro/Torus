#ifndef _resTorus_
#define _resTorus_
/* Class that provides tori of 1:1 resonant orbits computed from
   * p-theory (Binney 2016).  The creator requires at least two untrapped
   * orbits of a common E, one either side of the region of traping. The other arguments
   * are the actions of the underlying resonant torus, for this torus d Omega_N/d J1,
   * and the driving term in H, h_{2,-2) and its 1st 2 derivatives wrt J1=J_r-J_z.
   * Currently it is assumed that h_{2,-2)>0 so libration is around
   * theta_r-theta_z=0 or PI. Tori assume libraton around 0, so clockwise
   * circulation in (R,z) plane, but from these counter-clockwise circulating tori
   * are derived by theta_r -> TPI-theta_r, theta_z -> PI-theta_z.
   * Coded by James Binney March 2016 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "Torus.h"
#include "Potential.h"
#include "Numerics.templates"
#include "jjb_utils.h"

class resTorus{
	private:
		Torus *Tgrid; //Pointer to grid of fitted un-trapped tori
		double *Jrgrid;//Pointer to actions of above
		bool O;// Do trapped orbits circulate in Rz plane
		double lJ;// libration action
		double lOmega; //Frequency of libration
		double *t1grd,*trgrd,*Dgrd;//pointers to stores of t1 and corresponding t_true
		double omp,*hres;//dOmega1/dJ1 and pointers to h_{2,-2} & its 1st 2 derivatives
		int norb; // number of tori in input grid
		Actions Jb; //Actions of underlying resonant torus
		double t1c,Dt1;// t1=t_r-t_z librates by +-Dt1 about t1c
		void set_Dtheta1(void);// set Dt1
		void Dtheta1_fn(double,double&,double&) const;
		double theta_dot(double,double);// thetadot as f(t1,Jr-Jbr)
		void store_theta(void);//Constructs table of t1 versus libration angle
		int nangle;// size of above table
		double from_librationAngle(double,double&);//Get t1 and Jr-Jbr from libration angle
		bool search_eval(Position&,double,int,double&,double&,Angles&,Angles&);
		Angles lockit(Position&,double,double,double&,double&,int,int,double*);
		void getImin_max(void);// Set Imin,max
		void getDelta(double,double*);// find 2 values of Delta=Jr-Jbr & their difference
		double librationAngle(double,double);//return libration angle from (t1,Jr-Jbr).
		double rk_step(double&,double,double&,int);
	public:
		resTorus(Torus*,double*,int,Actions,double,double*);
		~resTorus(){delete[] t1grd; delete[] trgrd; delete[] Dgrd;};
		double I,Imin,Imax;
		void reset(Actions,double);
		int setI(double);
		double librationAction(void);// return resJ
		double librationOmega(void);// return libration frequency
		GCY FullMap(const Angles&);// (R,..,v_phi) a f(theta)
		void SOS(ostream&,const int);// write z=0 SoS to file
		void SOSr(ostream&,const int);// write z=0 SoS to file
		int containsPoint_Ang(Position&,Angles*,double*);//returns number of angles at which torus reaches given position
		int containsPoint(Position&,Velocity*,Angles*,double*);//last argument has values of t1
		void getJs(Actions&,Actions&);//returns extremes of unperturbed actions
		Actions actions(const Angles&);// return unperturbed actions at given point
};
#endif
