//into ForwardWithDerivatives

	if(!derivs_ok)
		TorusError(" `ToyIsochrone::Derivatives()' called without For/Backward",-2);

	double schi=sin(chi), cchi=cos(chi), sith=sin(th), csth=cos(th);
	double fac = at/HH, cGam= sqrt(1.-sGam*sGam);
// w. r. t. Jr:
	double Hx = wr, HxH = 0.5*Hx/H, wrx = 1.5*wr/H, wtx = wt*wrx;
	wrx*= wr;
	double ax = HxH/H, axa = ax/a, ex = pow(at/a,2)/(e*H2)*(HxH+axa);
	double exe = ex/e, axe = ax*e, exa = ex*a, psix=(exa+H2*axe)*spsi/(u+1.);
	double ux = u*axa-cpsi*exa+ae*spsi*psix;
	double wx =-wh*HxH + fac*(u+1.)/(r*r)*psix
		   -0.5*fac*(  Fint(a,-ae,ax,-axe-exa) + Fint(a+2.,-ae,ax,-axe-exa));
	double chix= tr/wr * (wrx*wt/wr - wtx) + wx;
	dQPdJ[0][0] = (u+1.)/r*ux;
	dQPdJ[1][0] = sGam*cchi*chix/csth;
	dQPdJ[2][0] = pr*(HxH -dQPdJ[0][0]/r +axa +exe) + HH*ae*cpsi/r*psix;
	dQPdJ[3][0] = at*sGam*(-schi/csth*chix +cchi*sith/pow(csth,2)*dQPdJ[1][0]);
	dQPdJ[0][0]/= sMob;
	dQPdJ[1][0]/= sMb;
	dQPdJ[2][0]/= b;
// w.r.t. Jt:
	Hx  = wt; HxH = 0.5*Hx/H; wrx = wtx; wtx = wt0r*wrx+2.*wr/pow(sq,3);
	ax  = HxH/H; axa = ax/a; ex  = pow(at/a,2) / (e*H2) * (HxH+axa-1./at);
	exe = ex/e; axe = ax*e; exa = ex*a; psix= (exa+H2*axe)*spsi/(u+1.);
	ux  = u*axa-cpsi*exa+ae*spsi*psix;
	wx  = wh*(1./at-HxH) + fac*(u+1.)/(r*r)*psix
	      -0.5*fac*(  Fint(a,-ae,ax,-axe-exa) + Fint(a+2.,-ae,ax,-axe-exa));
	chix= tr/wr * (wrx*wt/wr - wtx) + wx;
	if(jt==0.) {
		dQPdJ[1][1] = 0.;
		dQPdJ[3][1] = 100.;		// actually, it's infinite
	} else {
		double xGam= cGam/(sGam*at);
		dQPdJ[1][1] = (sGam*cchi*chix+schi*cGam*xGam)/csth;
		dQPdJ[3][1] = pt/at + at *( cchi*cGam/csth*xGam + sGam 
					    * (-schi/csth*chix +cchi*sith/pow(csth,2)*dQPdJ[1][1]));
	}
	dQPdJ[0][1] = (u+1.)/r*ux;
	dQPdJ[2][1] = pr*(HxH -dQPdJ[0][1]/r +axa +exe) + HH*ae*cpsi/r*psix;
	dQPdJ[0][1]/= sMob;
	dQPdJ[1][1]/= sMb;
	dQPdJ[2][1]/= b;
