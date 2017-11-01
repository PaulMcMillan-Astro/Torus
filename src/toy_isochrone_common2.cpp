derivs_ok = true;

// Extract and scale the actions and angles.
jr = double(JT(0)) / sMb;
jt = double(JT(1)) / sMb;
tr = double(JT(2));
tt = double(JT(3));
if(std::isnan(tr) || std::isinf(tr) || fabs(tr)>INT_MAX) 
	tr = 0.; // just in case  
	while(tr<0.)  tr+=TPi;
	while(tr>TPi) tr-=TPi;
	if(jr<0. || jt<0.|| this->jp<0.) {
		TorusError("ToyIsochrone: negative action(s)",2); 
		return fail;
	}
// Set up auxilliary variables independent of the angles tr and tt.
	at  = this->jp+jt;
	sGam= (this->jp==0.) ? 1. : sqrt(1.-pow(this->jp/at,2));
	sq  = hypot(2.,at);
	double fac = at/sq;
	HH  = 2./(2.*jr+at+sq);
	H2  = HH*HH;
	H   =-0.5*H2;
	a   = 1./H2-1.;
	if(at==0.) e = 1.;
	else {
		double e2    = 1. - pow(at/a,2)/H2;
		e     = (e2>0) ? sqrt(e2) : 0;
		if(e2<-1.e-3) TorusError("ToyIsochrone: bad e in ForwardWithDerivs",2);
	}
	ae  = a*e;
	eps = ae*H2;
	wr  = H2*HH;
	wt0r= 0.5*(1.+fac);
	wt  = wt0r*wr;
// Solve for the relevant other variables.
	psisolve();           // gives psi, cpsi, spsi
	double dw;
	wh  = wfundw(fac,dw); // fac is input, dw is output
	chi = tt-wt0r*tr+wh;
	u   = a*(1.-e*cpsi);
// Calculate co-ordinates and the conjugate momenta.
	double schi= sin(chi), cchi= cos(chi);
	r   = sqrt(u*(u+2.));
	double csth= cos(th = asin(sGam*schi));
	pr  = HH*ae*spsi/r;
	pt  = (at==0) ? 0. : at*sGam*cchi/csth;
// Calculate derivatives of the coordinates
	dQdT[0][0] = ae*spsi/r/H2*b;
	dQdT[0][1] = 0.;
	dQdT[1][1] = sGam*cchi/csth;
	double dchidtr = dw/(1.-eps*cpsi)-wt0r;
	dQdT[1][0] = dQdT[1][1] * dchidtr;
