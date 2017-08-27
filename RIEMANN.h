//-------------------------------------------------
// Include approximative Riemann-Solver from
// Y. Akamatsu et al., arXiv.org (2013) 34, 1302.1665
//-------------------------------------------------



//-------------------------------------------------
// Class for creation of left and right initial states
//-------------------------------------------------

class RSOLVER 
{ 
	EOS * eos;
	STATE IS;	
	public: 

	RSOLVER(EOS * eqn, STATE s): eos(eqn), IS(s)
	{
		iP=IS.pressure;
		iE=IS.energy;
		iCS=IS.soundvel;
		iVEL=IS.velx;
		iW=IS.gamma;
		iNB=IS.baryon;
		side=IS.side;

	}		

	double initial_pressure()
	{
		return iP;
	}

	double velocity(double T);
	double velocityderivative(double T); 
	double baryondensity(double T);
	double shockvel(double T);
	double wavespeed_min();
	double wavespeed_max();

	private:
	double iP,iE,iCS,iVEL,iW,iNB;
	int side;
	double baryonflux(double P, double E, double CS);
	double zeta(double P, double E, double CS);
};


//-------------------------------------------------
// Compute velocity
// Akamatsu et. al (2013) Eq: (30)
//-------------------------------------------------

double RSOLVER::velocity(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	double CS=eos->get_soundvel(T);

	double C=((iE+iP)*iW*iW);
	double A=((C*iVEL)+((P-iP)*zeta(P,E,CS)));
	double B=(C+((P-iP)*((iVEL*zeta(P,E,CS))+1)));
	return A/B;
}


//-------------------------------------------------
// Compute derivative of velocity
// Akamatsu et. al (2013) Eq: (35)
//-------------------------------------------------

double RSOLVER::velocityderivative(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	double CS=eos->get_soundvel(T);
	double zeta2=0;
	double A=(iE+iP)/(E+iP);
	double B=(iE+P)/(E+iP);
	double taub=A*(1-(B*(1/CS)));
	double D=(iE+iP)*iW*iW;
	zeta2=-0.5*iW*iW*((taub+(1/baryonflux(P,E,CS)))/((zeta(P,E,CS)*(1-(iVEL*iVEL)))-iVEL));
return (((zeta(P,E,CS)+zeta2)*(1-(iVEL*velocity(T))))-velocity(T))/(D+((P-iP)*((iVEL*zeta(P,E,CS))+1)));
}


//-------------------------------------------------
// Compute velocity
// Akamatsu et. al (2013) Eq: (27)
//-------------------------------------------------

double RSOLVER::baryondensity(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	return iNB*sqrt(((E+iP)*(E+P))/((iE+iP)*(iE+P)));
}


//-------------------------------------------------
// Compute square of Baryon flux
// p-pS and e-eS --> 0 leads to numerical incorrect
// values for J.
// In this case, switch to analytical value J(cS) with
// cS as the sound velocity
// Akamatsu et. al (2013) Eq: (29)
//-------------------------------------------------

double RSOLVER::baryonflux(double P, double E, double CS)
{
	double J=((E+iP)/(iE+iP))*((P-iP)/((E-iE)-(P-iP)));
	if(fabs(P-iP)<=1/fermi3 || fabs(E-iE)<=1/fermi3) 
	{
		return CS/(1-CS);
	}
	else
	{
        // Break, if J gets negative
        // Print values for error analysis
		if(J<0){
			cout << "NEGATIVE J" << endl;
			cout << "(E-iE): " << (E-iE)<< endl; 
			cout << "(P-iP): "<< (P-iP) << endl;
			cout << "iE: " << iE<< endl; 
			cout << "E: " << E<< endl; 
			cout << "iP: "<< iP << endl;
			cout << "P: "<< P << endl;
			exit(0);
		}
		return J;
	}

}


//-------------------------------------------------
// Compute zeta
// Akamatsu et. al (2013) Eq: (31)
// THe sign is chosen that
// + for R(1) and - L(0) 
//-------------------------------------------------

double RSOLVER::zeta(double P, double E, double CS)
{
	double sign=1.0;
	if(side==0) sign=-1.0; 
	else sign=1.0;
	return	(iVEL + (sign*sqrt(1+1/baryonflux(P,E,CS))))/(1-(iVEL*iVEL));
}


//-------------------------------------------------
// Compute velocity of shockwave
// A. Mignone, T. Plewa and G. Bodo, arXiv.org (2005), astro-ph/0505200v1
// Eq: (88)
//-------------------------------------------------

double RSOLVER::shockvel(double T) 
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	double CS=eos->get_soundvel(T);
	return (1/zeta(P,E,CS))+iVEL;	
}


//-------------------------------------------------
// Compute minimal wave speed
//-------------------------------------------------

double RSOLVER::wavespeed_min()
{
	return (iVEL-sqrt(iCS))/(1-(iVEL*sqrt(iCS)));
}


//-------------------------------------------------
// Compute maximal wave speed
//-------------------------------------------------

double RSOLVER::wavespeed_max()
{
	return (iVEL+sqrt(iCS))/(1+(iVEL*sqrt(iCS)));
}


//-------------------------------------------------
// Class for Riemann Solver that
// samples results on x/t
//-------------------------------------------------

class RIEMANN
{	
	EOS * eofstate;
	STATE ILEFT;
	STATE IRIGHT;
	RSOLVER * L, * R;

	public:

	RIEMANN (EOS * a, STATE l, STATE r) : eofstate(a), ILEFT(l), IRIGHT(r), L(new RSOLVER(a,l)), R(new RSOLVER(a,r)) {}

	~RIEMANN () 
	{
		delete L;
		delete R;
	}

	double pstar();
	void getflux(double * flux);
	private:

	void fill_star(RSOLVER * EX, STATE * S, double PS);
};


//-------------------------------------------------
// Search for pressure pstar in intermediate state
//-------------------------------------------------

double RIEMANN::pstar()
{
	int MAX=100;
	double TOLERANZ=10e-8;		
	double abb2=0;
	double piterate=0;
	double solution=0;
	double dif=0;
	double abb=0;
	double T;


    // Guess initial pressure
	piterate=0.5*(L->initial_pressure()+R->initial_pressure());

	T=eofstate->get_temperature(piterate);
	for(int i=0; i<= MAX; i++)
	{
		dif=((L->velocity(T)-R->velocity(T))/((L->velocityderivative(T)-R->velocityderivative(T))));
		solution=piterate-dif;				


		//-------------------------------------------------
        // Stop if difference of relative pressure smaller 
        // than error
		//-------------------------------------------------

		abb=(solution-piterate)/(0.5*(piterate+solution));
		if(fabs(abb)<=TOLERANZ) 
		{
			if (L->velocity(T)-R->velocity(T)>TOLERANZ) 
			{
				cout << "PSTAR: TEST FAILED!" << endl;
			}
			return piterate; 
		}
		else 
		{

            // Small routine that covers some convergence issues
			if(abb2<fabs(abb)+TOLERANZ && abb2>fabs(abb)-TOLERANZ)
			{
				return piterate;
			}

			else
			{
				abb2=abb;
				piterate=fabs(solution);
				T=eofstate->get_temperature(piterate);
			}	
		}
	}
	cout << "PSTAR: Reached maximal number of iterations!" << endl;
	cout << "PLEFT: " << L->initial_pressure() << endl;
	cout << "PRIGHT: " << R->initial_pressure() << endl;
	exit(0);
}



//-------------------------------------------------
// Fill structure of intermediate state
//-------------------------------------------------

void RIEMANN::fill_star(RSOLVER * EX,STATE * S, double PS)
{
	double T=0;
	T=S->temperature=eofstate->get_temperature(PS);
	S->pressure=PS;	
	S->energy=eofstate->get_energy(T);
	S->baryon=EX->baryondensity(T); 
	S->velx=EX->velocity(T);
	S->soundvel=eofstate->get_soundvel(T);	
	//Gamma-Faktor 1 Dim
	S->gamma=1/sqrt(1-(pow(S->velx,2)));
}


//-------------------------------------------------
// Sample solutions on x/t to get numerical flux
// A. Mignone, T. Plewa and G. Bodo, arXiv.org (2005), astro-ph/0505200v1
//-------------------------------------------------

void RIEMANN::getflux(double *flux)
{
	STATE SAMPLE, SSAMPLE, SOLUTION;	 
	RSOLVER * SAMPEX, * SSAMPEX;
	double PS=pstar();
	double LV,LVS,sig,m=0;

	//-------------------------------------------------
	// sigma= -sign(v*)
	//-------------------------------------------------

	if((L->velocity(eofstate->get_temperature(PS))+R->velocity(eofstate->get_temperature(PS)))/2 < 0) sig=1;
	else sig=-1;


	//-------------------------------------------------
    // Compute lambda, depending on the case if the wave
    // between S and S* is shock or rarefaction wave
	//-------------------------------------------------

	if(sig < 0) 
	{
		SOLUTION.side=0;
		SAMPLE=ILEFT;
		SSAMPLE.side=0;
		fill_star(L,&SSAMPLE,PS);
	}		
	else
	{
		SOLUTION.side=1;
		SAMPLE=IRIGHT;
		SSAMPLE.side=1;
		fill_star(R,&SSAMPLE,PS);
	}
	SAMPEX = new RSOLVER(eofstate, SAMPLE);	
	SSAMPEX = new RSOLVER(eofstate, SSAMPLE);	

	if(PS > SAMPLE.pressure)
	{
		LV=LVS=SAMPEX->shockvel(SSAMPLE.temperature);

	}
	else
	{
		if(sig <0) 
		{
			LVS=SSAMPEX->wavespeed_min();
			LV=min(SAMPEX->wavespeed_min(),LVS);
		}
		else 
		{
			LVS=SSAMPEX->wavespeed_max();
			LV=max(SAMPEX->wavespeed_max(),LVS);
		}
	}	


	if(sig*LVS >0)
	{
		SOLUTION=SSAMPLE;
	}
	if(sig*LV <0)
	{
		SOLUTION=SAMPLE;
	}
	if(sig*LVS <0 && 0< sig*LV)
	{
		SOLUTION.pressure=(LV*SSAMPLE.pressure-LVS*SAMPLE.pressure)/(LV-LVS);
		SOLUTION.energy=(LV*SSAMPLE.energy-LVS*SAMPLE.energy)/(LV-LVS);
		SOLUTION.velx=(LV*SSAMPLE.velx-LVS*SAMPLE.velx)/(LV-LVS);
		SOLUTION.baryon=(LV*SSAMPLE.baryon-LVS*SAMPLE.baryon)/(LV-LVS);
		SOLUTION.gamma=1/sqrt(1-pow(SOLUTION.velx,2));
	}

	//-------------------------------------------------
    // Compute numerical flux in conservative variables
	//1. D*v
	//2. m*v+p
	//3. m
	//-------------------------------------------------

	m=(SOLUTION.energy+SOLUTION.pressure)*pow(SOLUTION.gamma,2)*SOLUTION.velx;
	flux[0]=SOLUTION.gamma*SOLUTION.baryon*SOLUTION.velx;
	flux[1]=(m*SOLUTION.velx)+SOLUTION.pressure;	
	flux[2]=m;
}
