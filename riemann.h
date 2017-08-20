class EXACT 
{ 
	EOS * eos;
	STATE IS;	
	public: 

    	 EXACT(EOS * eqn, STATE s): eos(eqn), IS(s)
	{
		//cout << "T: " << IS.temperature << endl; 
		iP=IS.pressure;
		iE=IS.energy;
		iCS=IS.soundvel;
		iVEL=IS.velx;
		//optional for 3D Euler
		/*
		iVELy=IS.vely;
		iVELz=IS.velz;
		*/
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
	//double tang_vel_y(double T);
	//double tang_vel_z(double T);
	double shockvel(double T);
	double wavespeed_min();
	double wavespeed_max();

	private:
	double iP,iE,iCS,iVEL,iW,iNB;
	//double iVELy,iVELz;
	int side;
	double baryonflux(double P, double E, double CS);
	double zeta(double P, double E, double CS);
};

double EXACT::velocity(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	double CS=eos->get_soundvel(T);

	double C=((iE+iP)*iW*iW);
	double A=((C*iVEL)+((P-iP)*zeta(P,E,CS)));
	double B=(C+((P-iP)*((iVEL*zeta(P,E,CS))+1)));
	return A/B;
}

double EXACT::velocityderivative(double T)
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

double EXACT::baryondensity(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	return iNB*sqrt(((E+iP)*(E+P))/((iE+iP)*(iE+P)));
}

double EXACT::baryonflux(double P, double E, double CS)
		//for p-pS and e-eS --> 0 J becomes numerically inaccurate. so switch to analytical falue J(cS) what is small??
		{
		double J=((E+iP)/(iE+iP))*((P-iP)/((E-iE)-(P-iP)));
		if(fabs(P-iP)<=1/fermi3 || fabs(E-iE)<=1/fermi3) 
		{
		//	cout << "BFLUX->CS" << endl;
			return CS/(1-CS);
		}
		//THER SHOULD NOT BE A FABS()
			else
			{
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

double EXACT::zeta(double P, double E, double CS)
{
	// sign is choosen so that + for R(1) and - for L(0)
 	double sign=1.0;
	if(side==0) sign=-1.0; 
	else sign=1.0;
//cout << "BFLUX: " << baryonflux(P,E,CS) << endl;
//cout << "DEN: " << (iVEL + (sign*sqrt(1+1/baryonflux(P,E,CS)))) << endl;
//cout << "ENUM: " << (1-(iVEL*iVEL)) << endl;
//	return	(iVEL + (sign*sqrt(1+(((1-(iVEL*iVEL))*iW*iW)/baryonflux(P,E,CS)))))/(1-(iVEL*iVEL));
//	

return	(iVEL + (sign*sqrt(1+1/baryonflux(P,E,CS))))/(1-(iVEL*iVEL));

}

/*
double EXACT::tang_vel_y(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_pressure(T);
	double CS=eos->get_soundvel(T);
	double A=iW*(iE+iP);
	double D=(iW*iNB);
	
	return (A*iVELy)/(A+((P-iP)*((zeta(P,E,CS)*iVEL)+D)));
}
double EXACT::tang_vel_z(double T)
{
	double P=eos->get_pressure(T);
	double E=eos->get_pressure(T);
	double CS=eos->get_soundvel(T);
	double A=iW*(iE+iP);
	double D=(iW*iNB);
	
	return (A*iVELz)/(A+((P-iP)*((zeta(P,E,CS)*iVEL)+D)));
}
*/

double EXACT::shockvel(double T) 
{
	double P=eos->get_pressure(T);
	double E=eos->get_energy(T);
	double CS=eos->get_soundvel(T);
	//double D=(iW*iNB);
//cout << "SHOCKVEL:" << (D/zeta(P,E,CS))+iVEL<< endl;
//cout << "iV: " << iVEL << endl;
//	cout << "iNB: " << iNB << endl;
//	cout << "iW: " << iW << endl;
//	cout << "zeta: " << zeta(P,E,CS) << endl;
//
	return (1/zeta(P,E,CS))+iVEL;	//return (D/zeta(P,E,CS))+iVEL;
}

double EXACT::wavespeed_min()
{
	//Eigenvalues for 3D Euler
/*
	double DELTA=1-((pow(iVEL,2)+pow(iVELy,2)+pow(iVELz,2))*iCS);
	double ETHA=1-pow(iVEL,2)-(iCS*(pow(iVELy,2)+pow(iVELz,2)));
	return ((iW*iVEL*(1-iCS))-(sqrt(iCS*ETHA)))/(iW*DELTA);
*/
//	cout << "iVEL: " << iVEL << endl;
//	cout << "iCS: " << iCS << endl;
//	cout << "wavespeed_min: " << (iVEL-sqrt(iCS))/(1-(iVEL*sqrt(iCS))) << endl;
	return (iVEL-sqrt(iCS))/(1-(iVEL*sqrt(iCS)));
}

double EXACT::wavespeed_max()
{
/*
 	double DELTA=1-((pow(iVEL,2)+pow(iVELy,2)+pow(iVELz,2))*iCS);
	double ETHA=1-pow(iVEL,2)-(iCS*(pow(iVELy,2)+pow(iVELz,2)));
	return ((iW*iVEL*(1-iCS))+(sqrt(iCS*ETHA)))/(iW*DELTA);
*/
//	cout << "wavespeed_max: " << (iVEL+sqrt(iCS))/(1+(iVEL*sqrt(iCS))) << endl;
	return (iVEL+sqrt(iCS))/(1+(iVEL*sqrt(iCS)));
}


class RIEMANN
{	
	EOS * eofstate;
	STATE ILEFT;
	STATE IRIGHT;
	EXACT * L, * R;

	public:
	
	RIEMANN (EOS * a, STATE l, STATE r) : eofstate(a), ILEFT(l), IRIGHT(r), L(new EXACT(a,l)), R(new EXACT(a,r)) {}

	~RIEMANN () 
	{
		delete L;
		delete R;
	}
			
double pstar();
void getflux(double * flux);
	private:

void fill_star(EXACT * EX, STATE * S, double PS);
};


double RIEMANN::pstar()
{
#define MAX 100
#define TOLERANZ 10e-8		
	
	double abb2=0;
		
	double piterate=0;
	double solution=0;
	double dif=0;
	double abb=0;
	double T;


	//Guess pressure:
	
	piterate=0.5*(L->initial_pressure()+R->initial_pressure());
	T=eofstate->get_temperature(piterate);
	for(int i=0; i<= MAX; i++)
	{
		dif=((L->velocity(T)-R->velocity(T))/((L->velocityderivative(T)-R->velocityderivative(T))));
		solution=piterate-dif;				//stop iteration if relative pressure change is smaller than tolerance (Toro P. 127)
	abb=(solution-piterate)/(0.5*(piterate+solution));
	//	abb=solution-piterate;

		//different possibility for tolerance:
		//abb=L.velocity(T)-R.velocity(T);
	//	cout << fabs(abb) << "	" << i << endl;
		if(fabs(abb)<=TOLERANZ) 
		{
			if (L->velocity(T)-R->velocity(T)>TOLERANZ) 
			{
				cout << "PSTAR: TEST FAILED!" << endl;
			}
			return piterate; 
		}
		//fabs ist zu viel
		else 
		{
		
			//Additional Routine if TOLERANZ is too big for some Pressure Values on grid points	
			
			if(abb2<fabs(abb)+TOLERANZ && abb2>fabs(abb)-TOLERANZ)
			{
				//cout << "OUT OF BOUNDS!" << endl;
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
	return 0;
}

void RIEMANN::fill_star(EXACT * EX,STATE * S, double PS)
{
	double T=0;
	T=S->temperature=eofstate->get_temperature(PS);
	S->pressure=PS;	
	S->energy=eofstate->get_energy(T);
	S->baryon=EX->baryondensity(T); //Evtl. is das Falsch und ich brauche das
//nB aus den Rankin Hugenott Conditions
	S->velx=EX->velocity(T);
	S->soundvel=eofstate->get_soundvel(T);	
	//for RS with tangential velocities
	//S->vely=EX.tang_vel_y(T);
	//S->velz=EX.tang_vel_z(T);
	////gamma factor for 1Dim
	S->gamma=1/sqrt(1-(pow(S->velx,2)));
}
void RIEMANN::getflux(double *flux)
{
	//-----------
	//sample like Mignone et. al.
	//If that fails, try sampling like 
	//Marti & Mueller (1996)
	//-----------
	STATE SAMPLE, SSAMPLE, SOLUTION;	 
	EXACT * SAMPEX, * SSAMPEX;
	double PS=pstar();
	double LV,LVS,sig,m=0;
	//	fill_star(L,&LSTAR,PS);	
	//	fill_star(R,&RSTAR,PS);	

	//	printoutdebug(LSTAR,RSTAR);	


	//------------
	//sigma from Mignone et. al. in sampling procedure
	//sigma=- sign(v*)
	//------------
	if((L->velocity(eofstate->get_temperature(PS))+R->velocity(eofstate->get_temperature(PS)))/2 < 0) sig=1;
	else sig=-1;

	//-----------
	//compure lambda according whether the 
	//wave seperating S and S* is shock or rarefraction
	//-----------
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
	SAMPEX = new EXACT(eofstate, SAMPLE);	
	SSAMPEX = new EXACT(eofstate, SSAMPLE);	

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

//	cout << "LVS " << LVS << endl;
//	cout << "LV " << LV << endl;

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
//	cout << SOLUTION.pressure << endl;

	//-----------
	//Compute the Flux in conservative variables
	//1. D*v
	//2. m*v+p
	//3. m
	//-----------
	m=(SOLUTION.energy+SOLUTION.pressure)*pow(SOLUTION.gamma,2)*SOLUTION.velx;
	flux[0]=SOLUTION.gamma*SOLUTION.baryon*SOLUTION.velx;
	flux[1]=(m*SOLUTION.velx)+SOLUTION.pressure;	
	flux[2]=m;
}
