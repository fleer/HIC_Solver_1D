//-------------------------------------------------
// MUSCL-HANCOCK SCHEME
//-------------------------------------------------


#include <stdlib.h>

class MUSCL: public NUMSCHEME
{
	EOS * EOFSTATE;
	double CHEMPOT;

	int GRIDN;
	double *X,  *XL,  *XR, *DX, **GVAR;
	double **VARL, **VARR;
	double *LAMBD1, *LAMBD2, *LAMBD3;
	STATE *EFLEFT, *EFRIGHT; 
	string SLIMITER;
	public:
	MUSCL(EOS * F, int N, double S, STATE iLEFT1, STATE iLEFT2, STATE iRIGHT1, STATE iRIGHT2, STATE PERT, double **VAR, double *POSITION, double *XSTEP, string SL): EOFSTATE(F), GRIDN(N), SLIMITER(SL), SYSTEMLENGTH(S) 
	{
		
    //-------------------------------------------------
    // Parameters for Landau-Model
    // Define width of collisionsurface for AU-Nucelus
    // with R=6.5fm
    // Lorentz-Factor: 100
    // RHIC Energy --> sqrt(SNN)=200 MeV
    //-------------------------------------------------

		double NUCL=0.065;

    //-------------------------------------------------
    // 1/2 width of pertubation
    //-------------------------------------------------
		
		double PERTL=0.005;


    //-------------------------------------------------
    // Initialize variables
    // Fill the grid
    //-------------------------------------------------

		ZERO=4;
		STOP=GRIDN+4;
		X = new double [GRIDN+8];
		XL = new double [GRIDN+8];
		XR = new double [GRIDN+8];
		DX = new double [GRIDN+8];
		LAMBD1 = new double [GRIDN+8];
		LAMBD2 = new double [GRIDN+8];
		LAMBD3 = new double [GRIDN+8];
		EFLEFT= new STATE [GRIDN+8];
		EFRIGHT= new STATE [GRIDN+8];

		GVAR = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) GVAR[I] = new double [3];
		VARL = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) VARL[I] = new double [3];
		VARR = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) VARR[I] = new double [3];
		
		fill_grid();

		
    //-------------------------------------------------
    // Create start configuration
    //-------------------------------------------------

		for(int I=ZERO; I<STOP; I++)
		{
			if(X[I] < (SYSTEMLENGTH/2.-NUCL))
			{
				VAR[I][0]=iLEFT1.baryon;
				VAR[I][1]=iLEFT1.velx;
				VAR[I][2]=iLEFT1.pressure;
			}
			else
			if(X[I] > (SYSTEMLENGTH/2.-NUCL) && (X[I] < SYSTEMLENGTH/2.-PERTL || X[I] == SYSTEMLENGTH/2.-PERTL))
			{
				VAR[I][0]=iLEFT2.baryon;
				VAR[I][1]=iLEFT2.velx;
				VAR[I][2]=iLEFT2.pressure;
			}
			else
			if(X[I] > (SYSTEMLENGTH/2.-PERTL) && (X[I] <  SYSTEMLENGTH/2.+PERTL || X[I] <  SYSTEMLENGTH/2.+PERTL))
			{
				VAR[I][0]=PERT.baryon;
				VAR[I][1]=PERT.velx;
				VAR[I][2]=PERT.pressure;
			}
			else
			if((X[I] > SYSTEMLENGTH/2.+PERTL || X[I] > SYSTEMLENGTH/2.+PERTL) && X[I] < SYSTEMLENGTH/2.+NUCL)
			{
				VAR[I][0]=iRIGHT1.baryon;
				VAR[I][1]=iRIGHT1.velx;
				VAR[I][2]=iRIGHT1.pressure;
			}
			else
			{
				VAR[I][0]=iRIGHT2.baryon;
				VAR[I][1]=iRIGHT2.velx;
				VAR[I][2]=iRIGHT2.pressure;
			}
			POSITION[I]=X[I];
			XSTEP[I]=DX[I];
		}


//-------------------------------------------------
// Include boundary conditions
//-------------------------------------------------

		boundary(VAR);


//-------------------------------------------------
// Read in chemical potential and test that it is
// equal for both sides
//-------------------------------------------------

		if(iLEFT1.chempot != iRIGHT1.chempot)
		{
			cout << "SHOULD BE THE SAME ON BOTH SIDES OF PROBLEM!" << endl;
			exit(0);
		}
		CHEMPOT=iLEFT1.chempot;


	}


//-------------------------------------------------
// Destructor
//-------------------------------------------------

	~MUSCL() 
	{
		delete [] X;
		delete [] XL;
		delete [] XR;
		delete [] DX;

		delete [] *GVAR;
		delete [] GVAR;
		delete [] *VARR;
		delete [] VARR;
		delete [] *VARL;
		delete [] VARL;

	}

	void evaluate(double **VAR,double TIMESTEP);
	void boundary(double **VAR);
	double get_timestep();

	private:
	double SYSTEMLENGTH;
	int ZERO;
	int STOP;
	double DT;
	void fill_grid();
	void boundary();
	void get_lambda(double *LAMBD1, double *LAMBD2, double *LAMBD3);
	double timestep(double CFL, double iDT, double *LAMBD1, double *LAMBD3);
	void primrecovery(double *VAR, double *SCONS);
	double sign(double A);
	double minmod(double A, double B);
	double maxmod(double A, double B);
	double mc(double A, double B, double C);
	double superbee(double A, double B);
};



//-------------------------------------------------
// Compute the 3 eigenvalues for the equations of 
// ideal relativistic hydrodynamics in 1+1 dim
//-------------------------------------------------

void MUSCL::get_lambda(double *LAMBD1, double *LAMBD2, double *LAMBD3)
{


	for(int I=0; I<GRIDN+8; I++)
	{
		double T=EOFSTATE->get_temperature(GVAR[I][2]);

		LAMBD1[I] = (GVAR[I][1] - EOFSTATE->get_soundvel(T))/(1. - GVAR[I][1]*EOFSTATE->get_soundvel(T));
		LAMBD2[I] = GVAR[I][1];
		LAMBD3[I] = (GVAR[I][1] + EOFSTATE->get_soundvel(T))/(1. + GVAR[I][1]*EOFSTATE->get_soundvel(T));
	}	
}


//-------------------------------------------------
// Fill the grid
//-------------------------------------------------


void MUSCL::fill_grid()
{

	double DLX=SYSTEMLENGTH/((double) GRIDN);
	X[ZERO]=XL[ZERO]=XR[ZERO]=DX[ZERO]=0;

	for(int i=ZERO+1; i<STOP+1;i++)
	{
		XL[i]=XL[i-1]+DLX;
	}

	for(int i=ZERO; i<STOP;i++)
	{
		XR[i]=XL[i+1];
		X[i]=0.5*(XL[i]+XR[i]);
		DX[i]=XR[i]-XL[i];
	}
}


//-------------------------------------------------
// Compute next time step (CFL-Condition)
//-------------------------------------------------


double MUSCL::timestep(double CFL, double iDT, double *LAMBD1, double *LAMBD3)
{

	double DTCC=0;
	double DTEST=0;

	for(int I=ZERO; I < STOP; I++) 
	{
		double LAMBDA   = max(fabs(LAMBD1[I]),fabs(LAMBD3[I]));
		DTEST = LAMBDA/(XR[I] - XL[I]);
		if(DTEST > DTCC)
		{
			DTCC = DTEST;
		}
	}
	return min(CFL/DTCC, iDT);
}


//-------------------------------------------------
// Boundary conditions
//-------------------------------------------------

void MUSCL::boundary(double **VAR)
{
//-------------------------------------------------
//FLOW OUT
//-------------------------------------------------
	
		for(int I=1; I<=4; I++)
		for(int K=0; K<3; K++)
		{
			VAR[ZERO-I][K]=VAR[ZERO][K];
			VAR[STOP-1+I][K]=VAR[STOP-1][K];
		}
	
	for(int I=1; I<=4; I++)
	{
		DX[ZERO-I]=DX[ZERO+1];
		DX[STOP-1+I]=DX[STOP-1];
	}


}


double MUSCL::get_timestep()
{
	return DT; 
}	



//-------------------------------------------------
// Compute next time step
//-------------------------------------------------

void MUSCL::evaluate(double **VAR, double TIMESTEP)
{
	double GFLUX[GRIDN+8][3];
	double SCONS[GRIDN+8][3];
	double DTX[GRIDN+8];
	double FLUXL[GRIDN+8][3];
	double FLUXR[GRIDN+8][3];

	for(int I=0; I<GRIDN+8; I++)
	for(int K=0; K<3; K++)
	GVAR[I][K]=VAR[I][K];

	boundary(GVAR);

	get_lambda(LAMBD1,LAMBD2,LAMBD3);

	//CFL-Number= 0.4  
	DT=timestep(0.4,TIMESTEP,LAMBD1,LAMBD3);

	//-------------------------------------------------
    // Compute conservative variables
	//-------------------------------------------------

	for(int I=ZERO-2; I<STOP+2 ;I++) 
	{
		double W=1/sqrt(1-pow(GVAR[I][1],2));
		double E=EOFSTATE->get_energy(EOFSTATE->get_temperature(GVAR[I][2]));	
		DTX[I]=DT/DX[I];
		SCONS[I][0]=GVAR[I][0]*W;
		SCONS[I][1]=(E+GVAR[I][2])*pow(W,2)*GVAR[I][1];
		SCONS[I][2]=(E+GVAR[I][2])*pow(W,2)-GVAR[I][2];

	}

	//-------------------------------------------------
	// Slope-Limiter
	//-------------------------------------------------

	for (int I=ZERO-1; I<STOP+1; I++)
		for(int K=0; K<3 ;K++)
		{

			double DSCONSL,DSCONSR;
			double SLOPE=0;

			DSCONSL=SCONS[I][K]-SCONS[I-1][K];
			DSCONSR=SCONS[I+1][K]-SCONS[I][K];

			//minmod slope limiter
			if(SLIMITER== "MINMOD") SLOPE=minmod(DSCONSL/DX[I],DSCONSR/DX[I]); 
			//GODUNOV METHOD WITH SLOPE=0;
			else if(SLIMITER== "GODUNOV") SLOPE=0;	
			else
			{
				cout << "WRONG SLOPE LIMITER" << endl;
				exit(0);
			}


			//-------------------------------------------------
			// Reconstruction of UL and UR
			//-------------------------------------------------

			VARL[I][K]=SCONS[I][K]-SLOPE*0.5*DX[I];
			VARR[I][K]=SCONS[I][K]+SLOPE*0.5*DX[I];


		}







	//-------------------------------------------------
    // Reconstruction of primitive variables for VARL
	//-------------------------------------------------
#pragma omp parallel for
	for(int I=ZERO-1; I<STOP+1; I++) 
	{
		primrecovery(GVAR[I], VARL[I]);
	}


	//-------------------------------------------------
    // Compute regular flux FLUXL
	//-------------------------------------------------

	for (int I=ZERO-1; I<STOP+1; I++)
	{
		FLUXL[I][0]=VARL[I][0]*GVAR[I][1];
		FLUXL[I][1]=VARL[I][1]*GVAR[I][1]+GVAR[I][2];
		FLUXL[I][2]=VARL[I][1];
	}

	//-------------------------------------------------
    // Reconstruction of primitive variables for VARR
	//-------------------------------------------------

#pragma omp parallel for
	for(int I=ZERO-1; I<STOP+1; I++) 
	{
		primrecovery(GVAR[I], VARR[I]);
	}


	//-------------------------------------------------
    // Compute regular flux FLUXL
	//-------------------------------------------------

	for (int I=ZERO-1; I<STOP+1; I++)
	{
		FLUXR[I][0]=VARR[I][0]*GVAR[I][1];
		FLUXR[I][1]=VARR[I][1]*GVAR[I][1]+GVAR[I][2];
		FLUXR[I][2]=VARR[I][1];
	}


	//-------------------------------------------------
    // Evolution of UL and UR about delta t/2
	//-------------------------------------------------

	for (int I=ZERO-1; I<STOP+1; I++)
		for(int K=0; K<3 ;K++)
		{
			double F=DTX[I]*(FLUXL[I][K]-FLUXR[I][K])/2;
			VARL[I][K]=VARL[I][K]+F;	
			VARR[I][K]=VARR[I][K]+F;	
		}


	//-------------------------------------------------
    // Effective left states
	//-------------------------------------------------

#pragma omp parallel for
	for(int I=ZERO-1; I<STOP; I++) 
	{
		primrecovery(GVAR[I], VARR[I]);
	}
	for(int I=ZERO-1; I<STOP; I++)
	{	

		EFLEFT[I].side=0;
		EFLEFT[I].pressure = GVAR[I][2];
		EFLEFT[I].velx= GVAR[I][1];
		EFLEFT[I].baryon= GVAR[I][0]; 
		EFLEFT[I].temperature=EOFSTATE->get_temperature(EFLEFT[I].pressure);
		EFLEFT[I].energy=EOFSTATE->get_energy(EFLEFT[I].temperature); 
		EFLEFT[I].soundvel=EOFSTATE->get_soundvel(EFLEFT[I].temperature);
		EFLEFT[I].gamma=1/sqrt(1-pow(EFLEFT[I].velx,2));
		EFLEFT[I].chempot=CHEMPOT; 
	}


	//-------------------------------------------------
    // Effective right states
	//-------------------------------------------------

#pragma omp parallel for
	for(int I=ZERO-1; I<STOP+1; I++) 
	{
		primrecovery(GVAR[I], VARL[I]);
	}
	for(int I=ZERO-1; I<STOP; I++)
	{
		EFRIGHT[I].side=1;
		EFRIGHT[I].pressure = GVAR[I+1][2];  
		EFRIGHT[I].velx= GVAR[I+1][1];	
		EFRIGHT[I].baryon= GVAR[I+1][0];
		EFRIGHT[I].temperature=EOFSTATE->get_temperature(EFRIGHT[I].pressure);
		EFRIGHT[I].energy=EOFSTATE->get_energy(EFRIGHT[I].temperature); 
		EFRIGHT[I].soundvel=EOFSTATE->get_soundvel(EFRIGHT[I].temperature);
		EFRIGHT[I].gamma=1/sqrt(1-pow(EFRIGHT[I].velx,2));
		EFRIGHT[I].chempot=CHEMPOT; 
	}


#pragma omp parallel for
	for(int I=ZERO-1; I < STOP; I++) 
	{
		RIEMANN NFLUX(EOFSTATE, EFLEFT[I], EFRIGHT[I]);
		NFLUX.getflux(GFLUX[I]);
	}


	//-------------------------------------------------
    // Temporal evolution of effective variables
	//-------------------------------------------------

	for(int I=ZERO; I<STOP; I++) 
		for(int K=0; K<3; K++)
		{
			SCONS[I][K]=SCONS[I][K]+DTX[I]*(GFLUX[I-1][K]-GFLUX[I][K]);
		}


	//-------------------------------------------------
    // If baryon density --> 0
	//-------------------------------------------------

	for(int I=ZERO; I<STOP; I++) 
	{
		SCONS[I][0]=max(10e-12,SCONS[I][0]);
	}


	//-------------------------------------------------
    // Reconstruct primitive variables
	//-------------------------------------------------		

#pragma omp parallel for
	for(int I=ZERO; I<STOP; I++) 
	{
		primrecovery(GVAR[I], SCONS[I]);
	}

	for(int I=0; I<GRIDN+8; I++)
		for(int K=0; K<3; K++)
			VAR[I][K]=GVAR[I][K];
}


//-------------------------------------------------
// Algorithm for recovery of primitive variables
//-------------------------------------------------

	void MUSCL::primrecovery(double *VAR, double *SCONS)
{
	double T;
	VAR[2]=EOFSTATE->recovery(SCONS);

	T=EOFSTATE->get_temperature(VAR[2]);	
	VAR[0]=EOFSTATE->get_baryondensity(T,24);

	
//-------------------------------------------------
// Test if baryon density is unphysical
//-------------------------------------------------

	if(VAR[0] > 1) 
	{
		cout << "NB NOT PHYSICAL: " << VAR[0] << endl;
		VAR[0]=10e-10;
	}
	VAR[1]=SCONS[1]/(SCONS[2]+VAR[2]);
}
	double MUSCL::sign(double A)
{
	if(A < 0.) return -1.;
	else if(A > 0.) return 1.;
	else return 0.;
}


//---------
//MINMOD SLOPE LIMITER (KOLGAN 1972; VAN LEER 1979)
//REZZOLLA P. 428
//---------

	double MUSCL::minmod(double A, double B)
{
	double SIGNA;
	SIGNA=sign(A);
	return SIGNA*max(0.,min(fabs(A),B*SIGNA));
}

