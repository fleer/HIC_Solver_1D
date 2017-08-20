//NEU EINGEFUEGT AM 11. Juni: get_lambda unter boundary() in evaluationsschritt

#include <stdlib.h>

class PPM : public NUMSCHEME
{
	EOS * EOFSTATE;
	double CHEMPOT;

	int GRIDN;
	double *X,  *XL,  *XR, *DX, **GVAR;
	double **VARL, **VARR;
	double ***EFVARL, ***EFVARR;
	double **DELU;
	double *LAMBD1, *LAMBD2, *LAMBD3;
	STATE *EFLEFT, *EFRIGHT; 
	public:
	PPM(EOS * F, int N, double S, STATE iLEFT, STATE iRIGHT, double **VAR, double *POSITION, double *XSTEP): EOFSTATE(F), GRIDN(N), SYSTEMLENGTH(S)
	{
		//-----------
		//VARIABLE WHICH SETS 0-POINT BECAUSE OF DIFFERENT BOUNDARY CONDITIONS 
		//-3 To the Left +3 to the Right 
		//----------
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

		DELU= new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8 ;I++) DELU[I] = new double [3];
		GVAR = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) GVAR[I] = new double [3];
		VARL = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) VARL[I] = new double [3];
		VARR = new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) VARR[I] = new double [3];
		EFVARL = new double ** [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++)
		{
			EFVARL [I] = new double * [3];
			for(int K=0; K<3;K++) EFVARL[I][K] = new double [3];
		}
		EFVARR = new double ** [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++)
		{
			EFVARR [I] = new double * [3];
			for(int K=0; K<3;K++) EFVARR[I][K] = new double [3];
		}
		fill_grid();

		for(int I=ZERO; I<STOP; I++)
		{
			if(X[I] < SYSTEMLENGTH/2.)
			{
				VAR[I][0]=iLEFT.baryon;
				VAR[I][1]=iLEFT.velx;
				VAR[I][2]=iLEFT.pressure;
			}
			else
			{
				VAR[I][0]=iRIGHT.baryon;
				VAR[I][1]=iRIGHT.velx;
				VAR[I][2]=iRIGHT.pressure;
			}
			POSITION[I]=X[I];
			XSTEP[I]=DX[I];
		}

		boundary(VAR);


		//GET BARYON NUMBER --> SHOULD BE THE SAME ON BOTH SIDES OF PROBLEM!!!!
		if(iLEFT.chempot != iRIGHT.chempot)
		{
			cout << "SHOULD BE THE SAME ON BOTH SIDES OF PROBLEM!" << endl;
			exit(0);
		}
		CHEMPOT=iLEFT.chempot;


	}


	~PPM() 
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
	void interpolation();
	void detection(int K);
	void flaten();
	void monotonization();
	void get_average(int charact, double *LAMB,double DT);

	void get_effective(double *LAMB1,double *LAMB2, double *LAMB3);
};

void PPM::fill_grid()
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

void PPM::get_lambda(double *LAMBD1, double *LAMBD2, double *LAMBD3)
{


	for(int I=0; I<GRIDN+8; I++)
	{
		double T=EOFSTATE->get_temperature(GVAR[I][2]);

		LAMBD1[I] = (GVAR[I][1] - EOFSTATE->get_soundvel(T))/(1. - GVAR[I][1]*EOFSTATE->get_soundvel(T));
		LAMBD2[I] = GVAR[I][1];
		LAMBD3[I] = (GVAR[I][1] + EOFSTATE->get_soundvel(T))/(1. + GVAR[I][1]*EOFSTATE->get_soundvel(T));
	}	
}


/*
*     -------- 
*    NAME: timestep 
*     --------
*
*    PURPOSE: 
*    COMPUTES THE NEW TIMESTEP VALUE FROM COURANT CONDITION
*
*/

double PPM::timestep(double CFL, double iDT, double *LAMBD1, double *LAMBD3)
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


/*
 * RECONSTRUCTION PROCEDURE TO GET FIRST GUESS FOR A(R,J) & A(L,J)
 * STEP 1
 *INTERPOLATION TO GET A(J=1/2) 
 */

void PPM::interpolation()
{
	double COEFF1[GRIDN+8],COEFF2[GRIDN+8],COEFF3[GRIDN+8],COEFF4[GRIDN+8],COEFF5[GRIDN+8];
	double SCRCH[GRIDN+8],SCRCH1[GRIDN+8],SCRCH2[GRIDN+8],SCRCH3[GRIDN+8],SCRCH4[GRIDN+8];

	for(int I=ZERO-1; I<STOP+2 ; I++)
	{
		SCRCH1[I] = DX[I] + DX[I-1];
		SCRCH2[I] = SCRCH1[I] + DX[I];
		SCRCH3[I] = SCRCH1[I] + DX[I-1];
	} 

	for(int I=ZERO-1; I < STOP+1; I++) 
	{
		SCRCH4[I] = DX[I]/(SCRCH1[I] + DX[I+1]);
		COEFF1[I] = SCRCH4[I]*SCRCH3[I]/SCRCH1[I+1];
		COEFF2[I] = SCRCH4[I]*SCRCH2[I+1]/SCRCH1[I];
	}
	for(int I=ZERO-1; I<STOP; I++)
	{	
		SCRCH4[I] = 1./(SCRCH1[I] + SCRCH1[I+2]);
		COEFF3[I] = -SCRCH4[I]*DX[I]*SCRCH1[I]/SCRCH3[I+1];
		COEFF4[I] = SCRCH4[I]*DX[I+1]*SCRCH1[I+2]/SCRCH2[I+1];
		COEFF5[I] = DX[I] - 2.*(DX[I+1]*COEFF3[I] + DX[I]*COEFF4[I]);
		COEFF5[I] = COEFF5[I]/SCRCH1[I+1];
	}

	for(int K=0; K<3; K++) 
	{


		for(int I=ZERO-3; I<STOP+4; I++) SCRCH[I] = GVAR[I][K] - GVAR[I-1][K]; 


		//     ------------------------------------------------------
		//     DELU[I] AS IN EQ.61 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
		//     ------------------------------------------------------

		for(int I=ZERO-1; I<STOP+1; I++) 
			DELU[I][K] = COEFF1[I]*SCRCH[I+1] + COEFF2[I]*SCRCH[I];


		//     -------------------------------------------------------
		//     DELU(I) AS IN EQ.60 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
		//     -------------------------------------------------------

		for(int I=ZERO-1; I<STOP+1; I++) 
		{
			if(SCRCH[I+1]*SCRCH[I]>0)
			{
				double SDELU     = DELU[I][K]/fabs(DELU[I][K]);
				DELU[I][K] = min(min(fabs(DELU[I][K]),2.*fabs(SCRCH[I])),fabs(SCRCH[I+1]))*SDELU;
			}
			else DELU[I][K]   = 0;
		}

		//     ----------------------------------------------
		//     INTERFACE VALUES AS IN  EQ.59 OF MARTI AND MUELLER (1996), JCP, 
		//     VOL. 123, 1-14
		//     ----------------------------------------------


		for(int I=ZERO-1; I<STOP; I++) 
		{
			double	VP = GVAR[I][K]  + COEFF5[I]*SCRCH[I+1] + COEFF3[I]*DELU[I+1][K];
			VP   = VP + COEFF4[I]*DELU[I][K];
			VARR[I][K]=VP; 
			VARL[I+1][K]=VP; 
		}
	}

}

/*
 * RECONSTRUCTION PROCEDURE TO GET FIRST GUESS FOR A(R,J) & A(L,J)
 * STEP 2
 * DETECTION OF DISCONTINUITY & STEEPENING 
 */
void PPM::detection(int K)
{
	double SCRCH1[GRIDN+8],SCRCH2[GRIDN+8],SCRCH3[GRIDN+8],SCRCH4[GRIDN+8];
	double ETA[GRIDN+8],ETATIL[GRIDN+8];
	double D2U[GRIDN+8];

	//-----------
	//CHOOSE PARAMETERS FROM TABLE V IN MARTI & MUELLER 1996
	//-----------

	//double EPSLN=0.01, ETA1=6.4, ETA2=0.06;//, AK0=1.;
	//double EPSLN=0.01, ETA1=7, ETA2=0.06;//, AK0=1.;
	double EPSLN=0.1, ETA1=5, ETA2=0.05;//, AK0=1.;



	for(int I=ZERO-2; I<STOP+3; I++)
	{
		double A;
		SCRCH1[I] = DX[I] + DX[I-1];
		SCRCH2[I] = SCRCH1[I] + DX[I+1];
		A= GVAR[I][K] - GVAR[I-1][K];
		SCRCH1[I] = A/SCRCH1[I];
	}


	for(int I=ZERO-2; I<STOP+2; I++)
	{
		D2U[I]    = (SCRCH1[I+1] - SCRCH1[I])/SCRCH2[I];
		SCRCH4[I] = pow(X[I] - X[I-1],3);
	}



	for(int I=ZERO-1; I<STOP+1; I++)
	{
		SCRCH1[I] = D2U[I+1]*D2U[I-1];
		SCRCH3[I] = fabs(GVAR[I+1][K] - GVAR[I-1][K]) - EPSLN*min(fabs(GVAR[I+1][K]),fabs(GVAR[I-1][K]));
	}

	//     ------------------------------------------------------
	//     ETATIL(I) AS IN EQ.67 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
	//     -------------------------------------------------------

	for(int I=ZERO-1; I<STOP+1; I++)
	{	 
		// so that VARIABLE don't become 0
		if((GVAR[I+1][K]-GVAR[I-1][K]) == 0) SCRCH2[I]=10.e-5;
		else SCRCH2[I]=GVAR[I+1][K]-GVAR[I-1][K];

		if (SCRCH1[1] > 0 || SCRCH3[I] < 0) ETATIL[I]=0;
		else
			ETATIL[I] =( (D2U[I-1] - D2U[I+1])*(SCRCH4[I] + SCRCH4[I+1]))/(X[I+1] - X[I-1])/SCRCH2[I];
	}


	//     -------------------------------------------------------
	//     ETA(I) AS IN EQ.66 OF MARTI AND MUELLER (1996), JCP, VOL. 123, 1-14
	//     ONLY FOR ZONES VERIFYING EQ.63 (OTHERWISE, ZERO)
	//     -------------------------------------------------------

	for(int I=ZERO-1; I<STOP+1; I++)
	{
		ETA[I]    = max(0.,min(ETA1*(ETATIL[I] - ETA2),1.));
		//-----------------
		//CONDITION IS NORMALLY ONLY FOR GAS DYNAMICS
		//-----------------
		/*
		   double AK0=1;
		   double A,B;
		   A = fabs(GVAR[I+1][2] - GVAR[I-1][2])/min(GVAR[I+1][2],GVAR[I-1][2]);
		   B = fabs(GVAR[I+1][0] - GVAR[I-1][0])/min(GVAR[I+1][0],GVAR[I-1][0]);

		   if((AK0 * B )< A) ETA[I] = 0.;
		   */
	}

	//     ----------------------------
	//     NEW RECONSTRUCTED VALUES (EQ.65)
	//     ----------------------------

	for(int I=ZERO-1; I<STOP+1; I++)
	{
		double A,B;
		A = GVAR[I-1][K] + 0.5*DELU[I-1][K];
		B = GVAR[I+1][K] - 0.5*DELU[I+1][K];
		VARR[I][K]     = VARR[I][K]*(1.-ETA[I]) + A*ETA[I];
		VARL[I][K]     = VARL[I][K]*(1.-ETA[I]) + B*ETA[I];
	}
}

//     -------- 
//    NAME: F L A T E N  
//     --------

//    PURPOSE: 
//    THIS SUBROUTINE FLATTENS ZONE STRUCTURE IN REGIONS WHERE SHOCKS
//    ARE TOO THIN    


//    COMMENTS: 
//    STEP 3 IN THE RECONSTRUCTION PROCEDURE (SEE APPENDIX I IN MARTI 
//    & MUELLER 1996)
void PPM::flaten()
{
	double DP[GRIDN+8], DVEL[GRIDN+8], SCRCH1[GRIDN+8],SCRCH2[GRIDN+8],SCRCH3[GRIDN+8];
	//double EPSILN=9.0, w1=0.47, w2=4.8;
	double EPSILN=1.0, w1=0.52, w2=10.;


	for(int I=ZERO-3; I<STOP+3; I++)
	{
		DP[I]     = GVAR[I+1][2]   - GVAR[I-1][2];
		DVEL[I]   = GVAR[I+1][1] - GVAR[I-1][1];
		SCRCH1[I] = EPSILN*min(GVAR[I+1][2],GVAR[I-1][2]) - fabs(DP[I]);

		if(SCRCH1[I] < 0 && DVEL[I]<0) SCRCH1[I]=1.;
		else SCRCH1[I]=0;
	} 

	for(int I=ZERO-2; I<STOP+2; I++)
	{
		double DP2 = GVAR[I+2][2] - GVAR[I-2][2];

		if(DP2 == 0)
		{
			if(DP[I] ==0)
				SCRCH2[I]= - w1;
			else
				SCRCH2[I]=1.-w1;
		}
		else
			SCRCH2[I]= DP[I]/DP2 - w1;

		SCRCH3[I] = SCRCH1[I]*max(0.,SCRCH2[I]*w2);
	}



	for(int I=ZERO-1; I<STOP+1; I++)
	{
		if(DP[I] < 0)
			SCRCH2[I]=SCRCH3[I+1];
		else
			SCRCH2[I]=SCRCH3[I-1];
	}


	for(int I=ZERO-1; I<STOP+1; I++)
	{
		double FLATN,FLATN1;
		FLATN  = max(SCRCH2[I],SCRCH3[I]);
		FLATN  = max(0.,min(1.,FLATN));
		FLATN1 = 1. - FLATN;
		for(int K=0; K<3; K++)
		{
			VARL[I][K]=FLATN*GVAR[I][K] + FLATN1*VARL[I][K];
			VARR[I][K]=FLATN*GVAR[I][K] + FLATN1*VARR[I][K];
		}
	}
}

/*
 * RECONSTRUCTION PROCEDURE TO GET FIRST GUESS FOR A(R,J) & A(L,J)
 * STEP 4
 * MONOTONIZATION TO ENSURE THAT VELOCITIES ARE SMALLER THAN C 
 */

 void PPM::monotonization()
{

/*  
*     -----------------------------------------------------
*     NEW RECONSRUCTED VALUES IF CONDITION IN EQ.73 OF MARTI & MUELLER
*     1996 HOLDS
*     -----------------------------------------------------
*/

	for(int I=ZERO-1; I<STOP+1; I++)
	for(int K=0; K<3; K++)
	{
		double DVAR= VARR[I][K]-VARL[I][K];
		double SCRCH= VARR[I][K]-GVAR[I][K];
		double SCRCH1= GVAR[I][K]-VARL[I][K];
		double PVAR= VARL[I][K] + VARR[I][K];

		if(SCRCH*SCRCH1<0 || SCRCH*SCRCH1==0)
		{
			VARL[I][K]=GVAR[I][K];
			VARR[I][K]=GVAR[I][K];
		}
			if((DVAR*(GVAR[I][K]-PVAR/2))>pow(DVAR,2)/6)
			{
				VARL[I][K]=3*GVAR[I][K]-2*VARR[I][K];
			}
			if(-(DVAR*(GVAR[I][K]-PVAR/2))>pow(DVAR,2)/6)
			{
				VARR[I][K]=3*GVAR[I][K]-2*VARL[I][K];
			}

	}

}


/*     -------- 
*    NAME: A V R G 1 D 
*     --------

*    PURPOSE: 
*    THIS SUBROUTINE CALCULATES AVERAGES OF QUANTITIES P,RHO,VEL, OVER
*    THE PART OF THE DOMAIN OF DEPENDENCE FOR THE LAMBDA 
*    CHARACTERISTIC OF RADM(I) FOR THE TIME INTERVAL (T(N),T(N+1)).
*

*    COMMENTS: 
*    THIS ROUTINE CLOSELY FOLLOWS THE ANALYTICAL DEVELOPMENTS DESCRIBED IN 
*    MARTI & MUELLER, JCP, 1996
*/

void PPM::get_average(int charact, double *LAMB,double DT)
{

	double SCRCH1[GRIDN+8];
	double VAR6;

	for(int I=ZERO-1; I<STOP+1; I++)
	{
		SCRCH1[I] = max(0.,DT*LAMB[I]/DX[I]);
	}

		for(int I=ZERO ; I<STOP+1; I++)
		for(int K=0; K<3; K++)
		{ 
			VAR6=6.*GVAR[I-1][K]- 3.*(VARL[I-1][K]+VARR[I-1][K]);	
			EFVARL[I][charact][K]   = VARR[I-1][K]   - SCRCH1[I-1]/2.*(VARR[I-1][K]-VARL[I-1][K]  - (1. - 2.*SCRCH1[I-1]/3.)*VAR6 );
		}

	for(int I=ZERO-1; I<STOP+1; I++)
	{
		SCRCH1[I] = max(0.,-DT*LAMB[I]/DX[I]);
	}

		for(int I=ZERO ; I<STOP+1; I++)
		for(int K=0; K<3; K++)
		{ 
			VAR6=6.*GVAR[I][K]- 3.*(VARL[I][K]+VARR[I][K]);	
			EFVARR[I][charact][K]   = VARR[I][K]   + SCRCH1[I]/2.*(VARR[I][K]-VARL[I][K]  + (1. - 2.*SCRCH1[I-1]/3.)*VAR6 );
		}
}


/*
*     -------- 
*    NAME: S T A T 1 D 
*     --------
*
*    PURPOSE: 
*    THIS SUBROUTINE CALCULATES EFFECTIVE SECOND-ORDER-ACCURATE LEFT 
*    AND RIGHT STATES FOR RIEMANN PROBLEMS IN ONE DIMENSIONAL 
*    CALCULATIONS.
*
*
*    COMMENTS: 
*    THIS ROUTINE CLOSELY FOLLOWS THE ANALYTICAL DEVELOPMENTS DESCRIBED IN 
*    MARTI & MUELLER, JCP, 1996
*/

void PPM::get_effective(double *LAMB1,double *LAMB2, double *LAMB3)
{
	double BETAL[GRIDN+8][3], BETAR[GRIDN+8][3]; 
	double C[GRIDN+8],R[GRIDN+8];
	//     -------------
	//     EFFECTIVE LEFT STATES
	//     -------------

	for(int I=ZERO; I<STOP+1; I++)
	{
		double T=EOFSTATE->get_temperature(EFVARL[I][2][2]);

		double E=EOFSTATE->get_energy(T);
		double CS=EOFSTATE->get_soundvel(T);
		double NB=EFVARL[I][2][0];
		double P=EFVARL[I][2][2];
		double V=EFVARL[I][2][1];
		double GAMMAS=1/(1-pow(V,2));
		C[I]=(E+P)*GAMMAS*sqrt(CS);
		R[I]=NB*pow(GAMMAS,2)*(E+P);

		BETAL[I][0] = 0.5*(V - EFVARL[I][0][1] - (P - EFVARL[I][0][2])/C[I]);
		BETAL[I][1] =((NB-EFVARL[I][1][0])/R[I])-((P-EFVARL[I][1][2])/pow(C[I],2));  
		BETAL[I][2] = 0; 
	}

	for(int I=ZERO; I<STOP+1; I++)
	{
		if(LAMB1[I-1] <= 0)
			BETAL[I][0] = 0;

		if(LAMB2[I-1] <= 0)
			BETAL[I][1] = 0;
	}


	for(int I=ZERO; I<STOP+1; I++)
	{
		double T=EOFSTATE->get_temperature(EFVARL[I][2][2]);

		double E=EOFSTATE->get_energy(T);
		double CS=EOFSTATE->get_soundvel(T);
		double NB=EFVARL[I][2][0];
		double P=EFVARL[I][2][2];
		double V=EFVARL[I][2][1];
		double GAMMAS=1/(1-pow(V,2));
		C[I]=(E+P)*GAMMAS*sqrt(CS);
		R[I]=NB*pow(GAMMAS,2)*(E+P);
		EFLEFT[I].side=0;
		EFLEFT[I].pressure = P + C[I]*BETAL[I][0];  
		// PL(I)     = DMAX1(SMALLP,PL(I))

		EFLEFT[I].velx= V - BETAL[I][0];	
		//	cout << "V " << V << "	BETA " << BETAL[0][I] << endl; 
		EFLEFT[I].baryon= NB + R[I]*(BETAL[I][0]-BETAL[I][1]);
		EFLEFT[I].temperature=EOFSTATE->get_temperature(EFLEFT[I].pressure);
		EFLEFT[I].energy=EOFSTATE->get_energy(EFLEFT[I].temperature); 
		EFLEFT[I].soundvel=EOFSTATE->get_soundvel(EFLEFT[I].temperature);
		EFLEFT[I].gamma=1/sqrt(1-pow(EFLEFT[I].velx,2));
		EFLEFT[I].chempot=CHEMPOT; //Here I have to Work 
	}

	//     ----------
	//     EFFECTIVE RIGHT STATES
	//     ----------

	for(int I=ZERO; I<STOP+1; I++)
	{
		double T=EOFSTATE->get_temperature(EFVARR[I][0][2]);

		double E=EOFSTATE->get_energy(T);
		double CS=EOFSTATE->get_soundvel(T);
		double NB=EFVARR[I][0][0];
		double P=EFVARR[I][0][2];
		double V=EFVARR[I][0][1];
		double GAMMAS=1/(1-pow(V,2));
		C[I]=(E+P)*GAMMAS*sqrt(CS);
		R[I]=NB*pow(GAMMAS,2)*(E+P);

		BETAR[I][0] = 0; 
		BETAR[I][1] =((NB-EFVARR[I][1][0])/R[I])-((P-EFVARR[I][1][2])/pow(C[I],2));  
		BETAR[I][2] = -0.5*(V - EFVARR[I][0][1] + (P - EFVARR[I][0][2])/C[I]);
	}

	for(int I=ZERO; I<STOP+1; I++)
	{
		if(LAMB3[I] >= 0 )
			BETAR[I][2] = 0;

		if(LAMB2[I] >= 0 )
			BETAR[I][1] = 0;
	}


	for(int I=ZERO; I<STOP+1; I++)
	{
		double T=EOFSTATE->get_temperature(EFVARR[I][0][2]);

		double E=EOFSTATE->get_energy(T);
		double CS=EOFSTATE->get_soundvel(T);
		double NB=EFVARR[I][0][0];
		double P=EFVARR[I][0][2];
		double V=EFVARR[I][0][1];
		double GAMMAS=1/(1-pow(V,2));
		C[I]=(E+P)*GAMMAS*sqrt(CS);
		R[I]=NB*pow(GAMMAS,2)*(E+P);
		EFRIGHT[I].side=1;
		EFRIGHT[I].pressure = P + C[I]*BETAR[I][2];  
		// PL(I)     = DMAX1(SMALLP,PL(I))

		//	cout << "V " << V << "	BETA " << BETAL[0][I] << endl; 
		EFRIGHT[I].velx= V + BETAR[I][2];	
		EFRIGHT[I].baryon= NB + R[I]*(BETAR[I][2]+BETAR[I][1]);
		EFRIGHT[I].temperature=EOFSTATE->get_temperature(EFRIGHT[I].pressure);
		EFRIGHT[I].energy=EOFSTATE->get_energy(EFRIGHT[I].temperature); 
		EFRIGHT[I].soundvel=EOFSTATE->get_soundvel(EFRIGHT[I].temperature);
		EFRIGHT[I].gamma=1/sqrt(1-pow(EFRIGHT[I].velx,2));
		EFRIGHT[I].chempot=CHEMPOT; //Here I have to Work 

	}

}


void PPM::evaluate(double **VAR, double TIMESTEP)
{
	double GFLUX[GRIDN+8][3];
	double SCONS[GRIDN+8][3];
	double DTX[GRIDN+8];

		for(int I=0; I<GRIDN+8; I++)
		for(int K=0; K<3; K++)
			GVAR[I][K]=VAR[I][K];

		boundary(GVAR);

		//NEU EINGEFUEGT!!!!!!
		get_lambda(LAMBD1,LAMBD2,LAMBD3);

	for(int I=ZERO; I<STOP ;I++) 
	{
		double W=1/sqrt(1-pow(GVAR[I][1],2));
		double E=EOFSTATE->get_energy(EOFSTATE->get_temperature(GVAR[I][2]));	
		SCONS[I][0]=GVAR[I][0]*W;
		SCONS[I][1]=(E+GVAR[I][2])*pow(W,2)*GVAR[I][1];
		SCONS[I][2]=(E+GVAR[I][2])*pow(W,2)-GVAR[I][2];

	}
		//CFL Number is 0.4  accourding to details for FIG 8.
		DT=timestep(0.4,TIMESTEP,LAMBD1,LAMBD3);
		interpolation();
		detection(0);
		detection(1);
		detection(2);
		flaten();
		monotonization();
		get_average(0,LAMBD1,DT);
		get_average(1,LAMBD2,DT);
		get_average(2,LAMBD3,DT);
		get_effective(LAMBD1,LAMBD2,LAMBD3);


		for(int I=ZERO; I<STOP ;I++) 
			DTX[I]=DT/DX[I];

#pragma omp parallel for
		for(int I=ZERO; I < STOP+1; I++) 
		{
			//	cout << "I: " << I << endl; 
			//	cout << "PLEFT: " << EFLEFT[I].pressure*fermi3 << "	PRIGHT: " << EFRIGHT[I].pressure*fermi3 << endl;
			//		cout << "VELOCITYL: "<< EFLEFT[I].velx << endl;
			//		cout << "VELOCITYR: " << EFRIGHT[I].velx << endl;
			RIEMANN NFLUX(EOFSTATE, EFLEFT[I], EFRIGHT[I]);
			NFLUX.getflux(GFLUX[I]);
		}

		//-------
		//TIME ADVANCE
		//------

			for(int I=ZERO; I<STOP; I++) //HER SHOUD BE ZERO!
			for(int K=0; K<3; K++)
			{
				SCONS[I][K]=SCONS[I][K]-DTX[I]*(GFLUX[I+1][K]-GFLUX[I][K]);
			}
		//THEORETICALLU FOR NB =max(smallnb,NB)
		for(int I=ZERO; I<STOP; I++) 
		{
			SCONS[I][0]=max(10e-12,SCONS[I][0]);
		}

		//-------
		//PRIMITIVE RECOVERY
		//-------
#pragma omp parallel for
		for(int I=ZERO; I<STOP; I++) 
		{
			double T;
		//	double TEMP[3];
//			for(int J=0; J<3; J++) TEMP[J]=SCONS[I][J];
			//	cout << "I " << I-ZERO <<endl;
			//       cout << "D " << TEMP[0] << endl;
			//        cout << "M " << TEMP[1] << endl;
			//      cout << "E " << TEMP[2] << endl;
			GVAR[I][2]=EOFSTATE->recovery(SCONS[I]);
			T=EOFSTATE->get_temperature(GVAR[I][2]);	
			GVAR[I][0]=EOFSTATE->get_baryondensity(T,24);

			//NEE FCT for DETECT DISCONTINUITIES!!!!!!!!!!!!!!!!!!!
			if(GVAR[I][0] > 1) 
			{
				cout << "NB NOT PHYSICAL: " << GVAR[I][0] << endl;
				GVAR[I][0]=10e-10;
			}
			GVAR[I][1]=SCONS[I][1]/(SCONS[I][2]+GVAR[I][2]);

		}

	for(int I=0; I<GRIDN+8; I++)
		for(int K=0; K<3; K++)
			VAR[I][K]=GVAR[I][K];
}


void PPM::boundary(double **VAR)
{

	//--------
	//BOUNDARY CONDITIONS
	//-------
	//
	//-------
	//PERIODIC BOUNDARY
	//-------
	/*	
		for(int K=0; K<3; K++)
		{
		for(int I=1; I<=3; I++)
		{
		GVAR[K][ZERO-I]=GVAR[K][STOP-I];
		GVAR[K][STOP-1+I]=GVAR[K][ZERO-1+I];
		}

		}
		for(int I=1; I<=3; I++)
		{
		DX[ZERO-I]=DX[STOP-I];
		DX[STOP-1+I]=DX[ZERO-1+I];
		}
		*/
	//-------
	//FLOW OUT
	//-------
	
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


double PPM::get_timestep()
{
//CFL Number is 0.4  accourding to details for FIG 8.
		return DT; 
}	
