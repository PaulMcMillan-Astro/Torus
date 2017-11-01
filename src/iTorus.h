/* Code written by J Binney 2017 as an upgrade to the TorusModeller
 * MNRAS 456 1982 (2016).  Described in Binney (2017) */
#ifndef _iTorus_
#define _iTorus_
#include "eTorus.h"

class iTorus{
	private:
		eTorus **Tgrid;
		int np,nr;
		Actions J,Jbar,dJ;
		eTorus eT;
		double refine(const Position&,double&,double&);
		eTorus fixJ(const Angles&) const;
		bool oldV(Velocity*,Angles*,int,const Velocity&,const Angles&,double&) const;
	public:
		iTorus(Actions,eTorus**,int,int,Actions,Actions);
		void setJ(Actions);
		eTorus fixJ(const Angles&,Matrix<double,3,3>&) const;
		eTorus eT1(Actions) const;
		vec10 hn(void) const{return eT.hn();}
		int10 i1(void) const{return eT.i1();}
		int10 i2(void) const{return eT.i2();}
		int10 i3(void) const{return eT.i3();}
		Actions actions(void) const{return J;}
		Frequencies omega(void) const{return eT.omega();}
		void get_crit_Jp(const Position&,double&,double&,double,double);
		GCY FullMap(const Angles&) const;
		void SOS(ostream&,const int);
		bool InOrbit(const Position&);
		void LevMstart (Angles,const Position,PSPT&,Vector<double,3>&,
			       Matrix<double,3,3>&,Matrix<double,3,3>&,double&) const;
		double LevMstep (Angles&,const Position,PSPT&,Vector<double,3>&,
			      Matrix<double,3,3>&,Matrix<double,3,3>&,double&,double&) const;
		int containsPoint (const Position&,Velocity*,Angles*,double*,int) const;
};

#endif