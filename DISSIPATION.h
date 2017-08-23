//-------------------------------------------------
// Implementation dissipativer Effekte
//-------------------------------------------------

class DISSIPATION
{
	EOS *EOFSTATE;
	double **VHALF;	
	double *TAU;
	double ETA;
	double *VEL, *VELRIEMANN;
	int GRIDN;
	double *DX, *PXXEVOLVE, *BULKEVOLVE;
	double  DT;
	double *BULKCONST;
	public:
	DISSIPATION(EOS *e, double *t, double shear, double *V, double *VRIEMANN,  int GRIDPOINTS, double *X, double TIME, double *B): EOFSTATE(e), TAU(t), ETA(shear), VEL(V), VELRIEMANN(VRIEMANN) ,GRIDN(GRIDPOINTS), DX(X), DT(TIME), BULKCONST(B)
	{
		ZERO=4;
		STOP=GRIDN+4;

		VHALF= new double * [GRIDN+8];
		for(int I=0; I<GRIDN+8;I++) VHALF[I] = new double [3];

		PXXEVOLVE= new double [GRIDN+8];
		BULKEVOLVE= new double [GRIDN+8];

//-------------------------------------------------
// Randbedingungen
//-------------------------------------------------

		boundary(VEL);
		boundary(VELRIEMANN);
		boundary(TAU);
		boundary(DX);
		boundary(BULKCONST);
	}

	void get_PXX(double *PXX)
	{
		for(int I=ZERO; I<STOP; I++)
			PXX[I]=PXXEVOLVE[I];
	}
	void get_bulk(double *BULK)
	{
		for(int I=ZERO; I<STOP; I++)
			BULK[I]=BULKEVOLVE[I];
	}

	void evolve_dissip(double *PXX, double *BULK);

	void evolve_complete(double **VAR, double *PXX, double *BULK);

	private:
	int ZERO;
	int STOP;
	void boundary(double *VAR);

	double	DTIME(double EVELX, double V);
	double DPOS(double LVELX, double RVELX, double DEX);
	double PIXX(double DERIVX, double DERIVT, double V);
	double BULKNS(double DERIVX, double DERIVT, double V);
	double	UPWIND(double UL, double U, double UR, double VELX,double DELTAX);
};


//-------------------------------------------------
// Zeitliche Entwicklung der disipativen Groessen &
// Erzeugen des Flusses fuer die Entwicklung des gesamten Systems
//-------------------------------------------------

void DISSIPATION::evolve_dissip(double *PXX, double *BULK)
{
	boundary(PXX);
	boundary(BULK);
	//Variablen fuer Scherviskositaet
		double RELAX[GRIDN+8];
		double NSPIXX[GRIDN+8];

	//Variablen fuer Dehnviskositaet
		double BULKRELAX[GRIDN+8];
		double NSBULK[GRIDN+8];

		for(int I=ZERO-2; I<STOP+2;I++)
			{
				double DZ=DTIME(VELRIEMANN[I], VEL[I]);
				double DO=DPOS(VEL[I-1],VEL[I+1], DX[I]);	
				double BNS=BULKNS(DO,DZ,VEL[I]);
				NSPIXX[I]=PIXX(DO,DZ,VEL[I]);
				NSBULK[I]=-BULKCONST[I]*BNS;
			}

		for(int I=ZERO-2; I<STOP+2; I++)
		{
			double GAMMA=1/(1-pow(VELRIEMANN[I],2));
			RELAX[I]=UPWIND(PXX[I-1],PXX[I], PXX[I+1], VEL[I], DX[I]);
			RELAX[I]=((RELAX[I]-NSPIXX[I])*exp(-DT/TAU[I]))+NSPIXX[I];	

			BULKRELAX[I]=UPWIND(BULK[I-1],BULK[I], BULK[I+1], VEL[I], DX[I]);
			BULKRELAX[I]=((BULKRELAX[I]-NSBULK[I])*exp(-DT/TAU[I]))+NSBULK[I];


//-------------------------------------------------
// physikalische Fluesse
//-------------------------------------------------

			VHALF[I][0]=0;
			VHALF[I][1]=(BULKRELAX[I]*(1+GAMMA*pow(VELRIEMANN[I],2)))+RELAX[I];
			VHALF[I][2]=(BULKRELAX[I]*GAMMA*VELRIEMANN[I])+(RELAX[I]*VELRIEMANN[I]);
		}


//-------------------------------------------------
// Zeitliche Entwicklung der dissip Groessen
//-------------------------------------------------

		for(int I=ZERO-1; I<STOP+1; I++)
		{
			PXXEVOLVE[I]=UPWIND(RELAX[I-1],RELAX[I], RELAX[I+1], VELRIEMANN[I], DX[I]);
			BULKEVOLVE[I]=UPWIND(BULKRELAX[I-1],BULKRELAX[I], BULKRELAX[I+1], VELRIEMANN[I], DX[I]);
			
			
// Vermeidung von unphysikalischen Werten 
// (Zu kleine Werte fuehren zu schlechten Ergebnissen im Upwind-Scheme)
				if(fabs(PXXEVOLVE[I])<10e-12)
				PXXEVOLVE[I]=0;
		}	
}



//-------------------------------------------------
// Entwicklung des gesamten Systems
//-------------------------------------------------

void DISSIPATION::evolve_complete(double **VAR,double *PXX, double *BULK)
{
	double V [GRIDN+8];
	double UINITIAL[GRIDN+8][3];

	
//-------------------------------------------------
// Berechne konservative Variablen und die Geschwindigkeit 
// vor dem Anwenden des Riemann-Solvers
//-------------------------------------------------

	for(int I=ZERO; I<STOP ;I++) 
	{
		double W=1/sqrt(1-pow(VAR[I][1],2));
		double T=EOFSTATE->get_temperature(VAR[I][2]);
		double E=EOFSTATE->get_energy(T);	
		V[I]=VEL[I];
		UINITIAL[I][0]=VAR[I][0]*W;
		UINITIAL[I][1]=(E+VAR[I][2])*pow(W,2)*VAR[I][1];
		UINITIAL[I][2]=(E+VAR[I][2])*pow(W,2)-VAR[I][2];
	}

	for(int I=ZERO; I<STOP; I++)
	{
		double GAMMA=1/(1-pow(V[I],2));
		UINITIAL[I][0]=UINITIAL[I][0];
		UINITIAL[I][1]=UINITIAL[I][1]+BULK[I]*GAMMA*V[I]+PXX[I]*V[I];
		UINITIAL[I][2]=UINITIAL[I][2]+(GAMMA-1)*BULK[I]+PXX[I];
	}

	
//-------------------------------------------------
// Zeitliche Entwicklung des gesamten Systems
//-------------------------------------------------

	for(int I=ZERO; I<STOP; I++)
	{
		double DELTA=DT/(DX[I]);

		for(int K=0; K<3; K++)
		{
			double FMIN=0.5*(VHALF[I-1][K]+VHALF[I][K]);
			double FPLU=0.5*(VHALF[I+1][K]+VHALF[I][K]);
			UINITIAL[I][K]=(UINITIAL[I][K])+DELTA*(FMIN-FPLU);
		}
	}

	for(int I=ZERO; I<STOP; I++)
	{
		V[I]=VELRIEMANN[I];
	}


//-------------------------------------------------
// Wiederherstellung der primitiven Variablen mit
// Beruecksichtigung der dissipativen Effekte
//-------------------------------------------------

#pragma omp parallel for
	for(int I=ZERO; I<STOP; I++) 
	{
		int COND=0;
		double UDISSIPATION[3];

		double T1=0;
		while(COND == 0)
		{
		double GAMMA=1/(1-pow(V[I],2));
			UDISSIPATION[0]=0;  
			UDISSIPATION[1]=BULK[I]*GAMMA*V[I]+PXX[I]*V[I];
			UDISSIPATION[2]=(GAMMA-1)*BULK[I]+PXX[I];
			double T;
			double TEMP[3];

			for(int J=0; J<3; J++) 
			{
				if(fabs(UINITIAL[I][J])< 0.8*fabs(UDISSIPATION[J]))
				{
					cout << "UNPHYSICAL DISSIPATION!!!" << endl;
					TEMP[J]=UINITIAL[I][J];
				}
				else
					TEMP[J]=UINITIAL[I][J]-UDISSIPATION[J];
			}
			VAR[I][2]=EOFSTATE->recovery(TEMP);
			T=EOFSTATE->get_temperature(VAR[I][2]);	
			VAR[I][0]=EOFSTATE->get_baryondensity(T,24);

			//Vermeidung von zu kleinen Druecken
			if(VAR[I][2] < 10e-10) VAR[I][2]=0; 


//-------------------------------------------------
// Testen, ob Baryondichte oder Geschwindigkeit unphysikalisch ist
//-------------------------------------------------

			if(VAR[I][0] > 1) 
			{
				cout << "NB NOT PHYSICAL: " << VAR[I][0] << endl;
				VAR[I][0]=10e-10;
			}
			VAR[I][1]=TEMP[1]/(TEMP[2]+VAR[I][2]);

			if(VAR[I][1] > 1) cout << "UNPHYSICAL VELOCITY" << endl;

			if(fabs(T-T1) <10e-8) COND=1; 
			else 
			{
				T1=T;
				V[I]=VAR[I][1];
			}
		}
	}
}



//-------------------------------------------------
// Realisierung des upwind-schmes
//-------------------------------------------------

double DISSIPATION::UPWIND(double UL, double U, double UR, double VELX,double DELTAX)
		{
		//Wegen Strang Splitting nur halber Zeitschritt!
			double DELTA=(DT*VELX)/(DELTAX*2);
			if(VELX >0)
			{	
				return U-(DELTA*(U-UL));
			}
			else 
			{
			if(VELX <0)
			{	
				return U-(DELTA*(UR-U));
			}
			else
				return U;
			}
		}

 
//-------------------------------------------------
// Realisierung der Orts- und Zeitableitung fuer 
// die Navier-Stokes Groessen
//-------------------------------------------------

	double	DISSIPATION::DTIME(double EVELX, double V)
{
	return (EVELX-V)/DT;
}

	double DISSIPATION::DPOS(double LVELX, double RVELX, double DEX)
{
	return (RVELX-LVELX)/(2*DEX);
}

//-------------------------------------------------
// Navier-Stokes Groessen 
//-------------------------------------------------

	double DISSIPATION::PIXX(double DERIVX, double DERIVT, double V)
{
	double GAMMA=1/sqrt(1-pow(V,2));
	double A=(pow(GAMMA,3)*pow(V,2)+GAMMA)*DERIVX;
	double B=(pow(GAMMA,3)*pow(V,2)+GAMMA)*DERIVT;
	double C=pow(GAMMA,3)*V*DERIVT;
	double D= (1+pow(GAMMA,2)*pow(V,2))/3;
	return 2*ETA*(-A-GAMMA*V*(GAMMA*B+GAMMA*V*A)+D*(C+A));
}

	double DISSIPATION::BULKNS(double DERIVX, double DERIVT, double V)
{
	double GAMMA=1/sqrt(1-pow(V,2));
	double A=(pow(GAMMA,3)*pow(V,2)+GAMMA)*DERIVX;
	double C=pow(GAMMA,3)*V*DERIVT;
	return (C+A);
}	




//-------------------------------------------------
// Randbedingungen
//-------------------------------------------------

		void DISSIPATION::boundary(double *VAR)
{
	//-------------------------------------------------
	//FLOW OUT
	//-------------------------------------------------
	for(int I=1; I<=4; I++)
	{

		VAR[ZERO-I]=VAR[ZERO];
		VAR[STOP-1+I]=VAR[STOP-1];
	}
}
