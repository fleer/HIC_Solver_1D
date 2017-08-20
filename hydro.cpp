#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>
#include<limits>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
//#include<time.h>
#include "CPUTIME.h"
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using std::max;
using std::min;

//------
// STRUCTURE for the initial states of the problem
//------
struct STATE 
{
	STATE (int s) : side(s) {};
	STATE () {};
	double  temperature;
	double  chempot;
	double 	baryon;
	double  velx;
	// for 3Dim euler
	/*
	double  vely;
	double  velz;
	*/
	double  pressure;
	double  energy;
	double  gamma;
	double  soundvel;
	int side;
};	

//-------
// include EOS (IDEAL GAS, QCD (BOCHUM IMPLEMENTED)
//-------

#include "EOS.h"

//-------
// for usefull energy and pressure, multiply with fermi^3 (Mev^3/fm^3)
//-------

double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));


//-------
// functions which reads the input.ini anf prints it out in the terminal
//-------

void readin(STATE * L, STATE * R, double *STOPTIME, double *DELTAT, int *GRIDPOINTS, double *GRIDLENGTH, double *SHEARCONST, string *SCHEME, string *SLIMITER, string *WHICHEOS, string *VISCOUS); 

void printout(STATE L, STATE R, double STOPTIME, double DELTAT, int GRIDPOINTS, double GRIDLENGTH, double SHEARCONST, string SCHEME, string SLIMITER, string WHICHEOS, string VISCOUS);


//-------
//include NUMERICAL SCHEMES 
//-------


#include "NUMSCHEME.h"

//-------
//include DISSIPATION-METHOD (Akamatsu et. al 2013)
//-------

#include "dissipation.h"



int main()
{
	//-------
	// initialize inistial States
	//-------
	STATE iLEFT(0), iRIGHT(1);


	//-------
	//initialize equation of state
	//-------
	EOS * eofstate;

	//-------
	//initialize numerical scheme 
	//-------
	NUMSCHEME * nscheme;

	//-------
	//some variables needed for the specification of the programm parameters like
	//computing the dissipation, bulk voscosity, which eos, which numerical scheme etc.
	//-------

	string SCHEME, SLIMITER, WHICHEOS, VISCOUS;
	int yesno=0, bulkyes=0;

	//-------
	//initialize system dependend variables like grid length finial time of the system, # of gridoints
	//-------

	double STOPTIME, DELTAT, GRIDLENGTH;
	int GRIDPOINTS;
	STOPTIME=DELTAT=GRIDLENGTH=0;
	GRIDPOINTS=0;
	double SHEARCONST=0, *BULKCONST;
	double **GVAR;
	double *X , *DX;


	//-------
	//VARIABLES FOR DISSIPATION
	//-------

	double **GVARINITIAL, *VINITIAL,*VRIEMANN, *TAUETA, *PXX, *BULK;



	//-------
	//initialize time for measurement of program duration
	//-------

	//time_t now, later;
	double now, later;
	//int mod;



	//-------
	//check if OPENMP is available --> check number of cores
	//-------
#ifdef _OPENMP
	printf("OPENMP:\n");
	printf("Number of usable CPUs: %d \n", omp_get_num_procs());
#else	
	{
		printf("OpenMP is not supported. \n");
	}
#endif

	//-------
	//choose equation of state, the need of dissipation, bulk viscosity 
	//-------

	//time(&now);
	now=getCPUTime( );

	//-------
	//READ IN initial parameter from "input.ini"
	//-------

	readin(&iLEFT, &iRIGHT, &STOPTIME, &DELTAT, &GRIDPOINTS, &GRIDLENGTH, &SHEARCONST, &SCHEME, &SLIMITER, &WHICHEOS, &VISCOUS); 
	if(WHICHEOS== "QCD")
		eofstate = new QCD;
	else if(WHICHEOS== "GAS")
		eofstate = new GAS;
	else {
		cout << "WRONG INPUT!!!" << endl;
		return 0;
	}

	if(VISCOUS== "Y") 
		yesno=1;
	else if(VISCOUS== "BULK")
		bulkyes=yesno=1;
	else if(VISCOUS== "N")
		yesno=0;
	else {
		cout << "WRONG INPUT!!!" << endl;
		return 0;
	}
	iLEFT.pressure=eofstate->get_pressure(iLEFT.temperature);
	iLEFT.energy=eofstate->get_energy(iLEFT.temperature);
	//gamma factor 1-Dim	
	iLEFT.gamma=1/sqrt(1-(pow(iLEFT.velx,2)));
	//fermi3 zum testen!!!
	iLEFT.baryon=eofstate->get_baryondensity(iLEFT.temperature, iLEFT.chempot)*fermi3;
	iLEFT.soundvel=eofstate->get_soundvel(iLEFT.temperature);
	iRIGHT.pressure=eofstate->get_pressure(iRIGHT.temperature);
	iRIGHT.energy=eofstate->get_energy(iRIGHT.temperature);
	//gamma factor, only for 1Dim
	iRIGHT.gamma=1/sqrt(1-(pow(iRIGHT.velx,2)));
	//fermi3 zum testen!!!
	iRIGHT.baryon=eofstate->get_baryondensity(iRIGHT.temperature, iRIGHT.chempot)*fermi3;
	iRIGHT.soundvel=eofstate->get_soundvel(iRIGHT.temperature);
	printout(iLEFT, iRIGHT, STOPTIME, DELTAT, GRIDPOINTS, GRIDLENGTH, SHEARCONST, SCHEME, SLIMITER, WHICHEOS, VISCOUS);



	X = new double [GRIDPOINTS+8];
	DX = new double [GRIDPOINTS+8];
	GVAR = new double * [GRIDPOINTS+8];
	for(int I=0; I<GRIDPOINTS+8;I++) GVAR[I] = new double [3];

	//------
	//initialize grid
	//Size of Grid: [START-4:STOP+4]
	//------

	double ZERO=4;
	double STOP=GRIDPOINTS+4;

	//------
	//CONSERVED VARIABLES FOR DISSIPATION
	//------

	VINITIAL = new double  [GRIDPOINTS+8];
	VRIEMANN= new double  [GRIDPOINTS+8];
	GVARINITIAL= new double * [GRIDPOINTS+8];
	for(int I=0; I<GRIDPOINTS+8;I++) GVARINITIAL[I] = new double [3];

	//------
	//RELAXATION TIME
	//------	

	TAUETA= new double  [GRIDPOINTS+8];

	//------
	//VARIABLES FOR SHEAR (PXX) und BULK VISCOSITY (BULK)
	//------

	PXX= new double  [GRIDPOINTS+8];
	BULK= new double  [GRIDPOINTS+8];
	double ETA=SHEARCONST*pow(10,3)/(fermi3);		


	//------
	//CONST FOR BULK VISCOSITY (NOT CONSTANT BECAUSE OF THE DEPENDENCE ON SOUND VELOCITY) 
	//------

	BULKCONST= new double  [GRIDPOINTS+8];

	//INITIALIZE ARRAYS

	for(int I=0; I<GRIDPOINTS+8;I++)
	{
		TAUETA[I]=BULK[I]=BULKCONST[I]=PXX[I]=VINITIAL[I]=VRIEMANN[I]=0;
		for(int K=0; K<3; K++)
			GVARINITIAL[I][K]=0;
	}

	//-------
	//START time evaluation with the choosen numerical scheme (PPM, MUSCL); 
	// (EOS, GRID POINTS, GRID LENGTH, L initial STATE, R initial STATE, X, DELTA X)
	//-------

	if(SCHEME== "PPM")
		nscheme = new PPM (eofstate,GRIDPOINTS,GRIDLENGTH,iLEFT,iRIGHT, GVAR,X,DX);
	else {
		if(SCHEME== "MUSCL")
			nscheme = new MUSCL (eofstate,GRIDPOINTS,GRIDLENGTH,iLEFT,iRIGHT, GVAR,X,DX,SLIMITER);
		else {
			cout << "(SCHEME) WRONG INPUT!!!" << endl;
			return 0;
		}
	}


	float progress = 0.0;
	int FN=1;

	double FTIME=0;

	//-------
	//output for the contourplot for p & v
	//-------
	fstream c_output_p, c_output_v, c_output_n;

	c_output_p.open("DATA/c_output_p.txt", ios::out);
	c_output_v.open("DATA/c_output_v.txt", ios::out);
	c_output_n.open("DATA/c_output_n.txt", ios::out);


	//-------
	//output of initial state in file "output_0.txt"
	//-------

	fstream output;

	output.open("DATA/output_0.txt", ios::out);

	for(int I=ZERO; I<STOP; I++)
	{
		double OUTTEMP=eofstate->get_temperature(GVAR[I][2]);
		output << X[I] << "\t";
		output << GVAR[I][0] << "\t";
		output << GVAR[I][1] << "\t";
		output << GVAR[I][2]*fermi3 << "\t ";
		output << eofstate->get_energy(OUTTEMP)*fermi3 << "\t";
		output << OUTTEMP;
		output << endl;

		c_output_p << X[I] << "\t";
		c_output_p << FTIME << "\t";
		c_output_p << GVAR[I][2]*fermi3 << "\t";
		c_output_p << endl;

		c_output_v << X[I] << "\t";
		c_output_v << FTIME << "\t";
		c_output_v << GVAR[I][1] << "\t";
		c_output_v << endl;

		c_output_n << X[I] << "\t";
		c_output_n << FTIME << "\t";
		c_output_n << GVAR[I][0] << "\t";
		c_output_n << endl;
	}
	c_output_p << endl;
	c_output_v << endl;	
	c_output_n << endl;	
	output.close();



	int barWidth=70;

	//-------
	//start of time loop
	//-------

	while(	FTIME <= STOPTIME)
	{
		//-------
		//creation of progress bar
		//-------

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "█";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress*100.0) << " %\r"; 
		std::cout.flush();

		nscheme->boundary(GVAR);
		//-------
		//SAVE PRIMITIVE VARIABLES FOR DISSIPATION & GET RELAXATION TIME
		//-------
		for(int I=ZERO; I< STOP; I++)
		{
			double FERMI=0.197*pow(10,3);
			double T=eofstate->get_temperature(GVAR[I][2]);
			//TAUETA[I]=(SHEARCONST*10*FERMI)/T;
			TAUETA[I]=1/(T/FERMI);
			if(bulkyes==1)
			{
				double C=1/3-eofstate->get_soundvel(T);
				BULKCONST[I]=15*pow(C,2)*ETA;
			}
			//CHECK IF DELTAT IS SMALL ENOUGH, BECAUSE DT<TAU!!!

			if(DELTAT>TAUETA[I] || DELTAT==TAUETA[I]) DELTAT=TAUETA[I]*0.1;

			VINITIAL[I]=GVAR[I][1];
			for(int K=0; K<3;K++)
				GVARINITIAL[I][K]=GVAR[I][K];
		}



		nscheme->evaluate(GVAR,DELTAT);

		//-------
		//SAVE PRIMITIVE VARIABLES FOR DISSIPATION
		//-------
		for(int I=ZERO; I< STOP; I++)
			VRIEMANN[I]=GVAR[I][1];


		//-------
		//DISSIPATION 
		//-------
		if(yesno==1)
		{
			DISSIPATION dissipation(eofstate, TAUETA, ETA,VINITIAL,VRIEMANN, GRIDPOINTS, DX, DELTAT, BULKCONST);

			dissipation.evolve_dissip(PXX,BULK);

			dissipation.get_PXX(PXX);
			dissipation.get_bulk(BULK);

			dissipation.evolve_complete(GVARINITIAL,PXX,BULK);
			nscheme->boundary(GVARINITIAL);

			nscheme->evaluate(GVARINITIAL,DELTAT);
			dissipation.evolve_complete(GVARINITIAL,PXX,BULK);

			for(int I=ZERO; I< STOP; I++)
			{
				for(int K=0; K<3;K++)
					GVAR[I][K]=GVARINITIAL[I][K];
			}
		}

		//CFL Number is 0.4  according to details for FIG 8.
		FTIME=FTIME+nscheme->get_timestep();
		

		//-------
		//output
		//------

		fstream output;

		output.open("DATA/output_"+std::to_string(FN)+".txt", ios::out);

		//output.open("out.dat", ios::out);
		for(int I=ZERO; I<STOP; I++)
		{
			double OUTTEMP=eofstate->get_temperature(GVAR[I][2]);
			output << X[I] << "\t";
			output << GVAR[I][0] << "\t";
			output << GVAR[I][1] << "\t";
			output << GVAR[I][2]*fermi3 << "\t";
			output << eofstate->get_energy(OUTTEMP)*fermi3 << "\t";
			output << OUTTEMP;	
			output << endl;


			c_output_p << X[I] << "\t";
			c_output_p << FTIME << "\t";
			c_output_p << GVAR[I][2]*fermi3 << "\t";
			c_output_p << endl;


			c_output_v << X[I] << "\t";
			c_output_v << FTIME << "\t";
			c_output_v << GVAR[I][1] << "\t";
			c_output_v << endl;

			c_output_n << X[I] << "\t";
			c_output_n << FTIME << "\t";
			c_output_n << GVAR[I][0] << "\t";
			c_output_n << endl;

		}
		c_output_p << endl;
		c_output_v << endl;	
		c_output_n << endl;	
		output.close();

		FN++;
		progress=FTIME/STOPTIME;
	}
	c_output_p.close();
	c_output_v.close();
	c_output_n.close();
	cout << endl;




	//time(&later);
	later = getCPUTime( );
	//mod= ((int) difftime(later,now))% 60;
	
	cout << "Computation Time: " << (later - now) << " seconds!" << endl; 
	return 0;
}

//-------
//READ IN the initial parameter from 'input.ini'
//-------
void readin(STATE * L, STATE * R, double *STOPTIME, double *DELTAT, int *GRIDPOINTS, double *GRIDLENGTH, double *SHEARCONST, string *SCHEME, string *SLIMITER, string *WHICHEOS, string *VISCOUS) 
{

	std::string templine; 
	ifstream sourcedata;
	sourcedata.open("input.ini", ios_base::in);

if (!sourcedata)
       cout << "Could not open data!" << endl;
else 
{                                       
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	L->temperature=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	L->chempot=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	L->velx=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	/*		
			L->vely=atof(templine.c_str());
			sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
			getline(sourcedata,templine);
			L->velz=atof(templine.c_str());
			sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
			getline(sourcedata,templine);
	*/
	R->temperature=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	R->chempot=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	R->velx=atof(templine.c_str());
	/*
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	   R->vely=atof(templine.c_str());
	   sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	   getline(sourcedata,templine);
	   R->velz=atof(templine.c_str());
	*/
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*STOPTIME=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*DELTAT=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*GRIDPOINTS=(int) atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*GRIDLENGTH=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*SHEARCONST=atof(templine.c_str());
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*SCHEME=templine.c_str();
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*SLIMITER=templine.c_str();
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*WHICHEOS=templine.c_str();
	sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
	getline(sourcedata,templine);
	*VISCOUS=templine.c_str();

} 
sourcedata.close();
}

void printout(STATE L, STATE R, double STOPTIME, double DELTAT, int GRIDPOINTS, double GRIDLENGTH, double SHEARCONST, string SCHEME, string SLIMITER, string WHICHEOS, string VISCOUS)
{
	cout << "---------------------" << endl;
	cout << "LEFT STATE:" << endl;
	cout << "---------------------" << endl;

	cout << "INITIAL TEMPERATURE [MeV]:" << endl; 
	cout << L.temperature << endl;
	cout << "INITIAL CHEMICAL POTENTIAL [MeV]:" << endl;
	cout << L.chempot << endl;
	cout << "PRESHOCK PRESSURE [MeV/fm^3]:";
	cout << "\n";
	cout << L.pressure*fermi3;
	cout << "\n";
	cout << "PRESHOCK ENERGYDENSITY [MeV/fm^3]:";
	cout << "\n";
	cout << L.energy*fermi3;
	cout << "\n";
	cout << "PRESHOCK BARYONDENSITY [MeV]:";
	cout << "\n";
	cout << L.baryon;
	cout << "\n"; 
	cout << "INITIAL VELOCITY X (UNITS OF C):" << endl; 
	cout << L.velx<< endl; 
/*
	cout << "INITIAL VELOCITY Y (UNITS OF C):" << endl; 
	cout << L.vely<< endl; 
	cout << "INITIAL VELOCITY Z (UNITS OF C):" << endl; 
	cout << L.velz<< endl; 
*/
	cout << "INITIAL SOUNDVElOCITY X (UNITS OF C):" << endl; 
	cout << L.soundvel << endl;
	cout << "GAMMA FACTOR:" << endl;
	cout << L.gamma<< endl;
	cout << "\n";

	cout << "---------------------" << endl;
	cout << "RIGHT STATE:" << endl;
	cout << "---------------------" << endl;
	cout << "INITIAL TEMPERATURE [MeV]:" << endl; 
	cout << R.temperature << endl;
	cout << "INITIAL CHEMICAL POTENTIAL [MeV]:" << endl;
	cout << R.chempot << endl;
	cout << "PRESHOCK PRESSURE [MeV/fm^3]:";
	cout << "\n";
	cout << R.pressure*fermi3;
	cout << "\n";
	cout << "PRESHOCK ENERGYDENSITY [MeV/fm^3]:";
	cout << "\n";
	cout << R.energy*fermi3;
	cout << "\n";
	cout << "PRESHOCK BARYONDENSITY [MeV]:";
	cout << "\n";
	cout << R.baryon;
	cout << "\n"; 
	cout << "INITIAL VELOCITY X (UNITS OF C):" << endl; 
	cout << R.velx<< endl; 
/*
	cout << "INITIAL VELOCITY Y (UNITS OF C):" << endl; 
	cout << R.vely<< endl; 
	cout << "INITIAL VELOCITY Z (UNITS OF C):" << endl; 
	cout << R.velz<< endl; 
*/
	cout << "INITIAL SOUNDVELOCITY X (UNITS OF C):" << endl; 
	cout << R.soundvel << endl;
	cout << "GAMMA FACTOR:" << endl;
	cout << R.gamma<< endl;
	cout << "\n";

	cout << "---------------------" << endl;
	cout << "TIME:" << endl;
	cout << "---------------------" << endl;
	cout << "TIME FOR EVALUTATION [fm]: " << STOPTIME << endl;
	cout << "ESTIMATED TIMESTEP [fm]: " << DELTAT << endl;
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "GRID:" << endl;
	cout << "---------------------" << endl;
	cout << "NUMBER OF GRIDPOINTS: " << GRIDPOINTS << endl;
       cout << "LENGTH OF THE SYSTEM [fm]: " << GRIDLENGTH << endl;	
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "DISSIPATION:" << endl;
	cout << "---------------------" << endl;
	cout << "SHEARVISCOSITY ETA/S [1/fm^3]: " << SHEARCONST << endl;
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "NUMERICAL SCHEME:" << endl;
	cout << "---------------------" << endl;
	cout << SCHEME << endl;
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "SLOPE LIMITER FOR MUSCL SCHEME:" << endl;
	cout << "---------------------" << endl;
	cout << SLIMITER << endl;
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "EQUATION OF STATE:" << endl;
	cout << "---------------------" << endl;
	cout << WHICHEOS << endl;
	cout << "\n";
	cout << "---------------------" << endl;
	cout << "DISSIPATION:" << endl;
	cout << "---------------------" << endl;
	cout << VISCOUS << endl;
}



