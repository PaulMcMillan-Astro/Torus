#ifndef _bar_pot_
#define _bar_pot_

class bar_pot{
	private:
		double fourPiG;
		double K;
		double q;
		double Rb,Rb2;
		inline double m2(double R,double z){return R*R+pow(z/q,2);};
	public:
		bar_pot(double qin,double Rbin,double rho0in){
			q=qin; Rb=Rbin; Rb2=Rbin*Rbin;
			fourPiG=5.65318158e-11;
			K=fourPiG*rho0in*pow(Rb,5);
		}
		bar_pot(double qin,double Rbin,double Vc,double A){
			q=qin; Rb=Rbin; Rb2=Rbin*Rbin;
			fourPiG=5.65318158e-11;
			K=pow(Vc,2)*pow(Rb,3)*A;
		}
		void reset_K(double);
		void Sormani(void);
		double Phi(double,double);
		double Phi(double,double,double&,double&);
		double dPhidR(double,double);
		double oRddRRdPhidR(double,double);
		double dPhidz(double,double);
		double d2Phidz2(double,double);
		double rho(double,double);
};

#endif