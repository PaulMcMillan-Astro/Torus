/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
#ifndef _eTorus_
#define _eTorus_
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "Torus.h"
#include "Numerics.templates"
#include "PJMebf.h"
#include "PJMCoords.h"
#include "PJM_cline.h"
#include "jjb_utils.h"
#include "falPot.h"
#include "bar_pot.h"
#include "ebf.hpp"

class eTorus{
	private:
		Torus T;
		int tmax;
		double Omp;
		vec10 hN;
		int10 iN,jN,kN;
		double H0(Potential*,bar_pot*,Angles&);
		void getH0Hn(Potential*,bar_pot*,int,int);
		int fq(int,int);
		double insert(double*,int&,double,int,int,int);
	public:
		eTorus(void){};
		eTorus(Torus&,Potential*,bar_pot*,const double);
		eTorus(Actions,Potential*,bar_pot*,const double,const double);
		int AutoFit(Actions,Potential*,bar_pot*,const double,const double);
		void reset(Torus&,Potential*,bar_pot*);
		Frequencies omega(void) const;
		Actions actions(void) const{return T.actions();}
		int10 i1(void) const{return iN;}
		int10 i2(void) const{return jN;}
		int10 i3(void) const{return kN;}
		vec10 hn(void) const{return hN;}
		GCY FullMap(Angles A) const{return T.FullMap(A);}
		double DToverDt(Angles A,Angles &AT) const{return T.DToverDt(A,AT);}
		void LevCof3D(const Angles&,Matrix<double,3,3>&,const Position&,PSPT&,double&,
			      Vector<double,3>&,Matrix<double,3,3>&,Matrix<double,3,3>&) const;
		void LevCof3DTrue(const Angles&,Matrix<double,3,3>&,const Position&,PSPT&,double&,
			      Vector<double,3>&,Matrix<double,3,3>&) const;
		double DistancetoPoint(const Position&,Angles&) const;
		bool containsPoint (const Position&,Velocity&,Velocity&,double&,Angles&,Angles&,
				    Velocity&,Velocity&,double&,Angles&,Angles&)
				    const;
		bool write_ebf(const string,const string);
		bool read_ebf(const string,const string);
		~eTorus(){};
		eTorus&	operator= (const eTorus&);
		eTorus&	operator*= (const double&);
		eTorus&	operator+= (const eTorus&);
		eTorus&	operator-= (const eTorus&);
		const eTorus operator* (const double&);
		const eTorus	 operator+ (const eTorus&);
		const eTorus	 operator- (const eTorus&);
};
eTorus Interp_eTorus_mn(eTorus**,int,int,Actions,Actions,Actions);

#endif