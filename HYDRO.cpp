//-------------------------------------------------
//-------------------------------------------------
// Program for the numerical simulation of Quark-Gluon-Plasma
// in heavy-ion collisions, utilizing relativistic hydrodynamics
//
//  Building Blocks:
//      + Computations in 1+1 Dimensions
//      + MUSCL-HANCOCK Scheme
//      + Can take dissipative effects under consideration
//          (using the ISRAEL-STEWARD-THEORY)
//      + Variable initial conditions:
//          - Shocktube Problem
//          - Landau-Model
//          - Integration of small pertubation
//  Compile via:
//	g++ HYDRO.cpp -std=c++11 -Wall -o hydro -lm -lgsl -fopenmp 
//
//
//  Author: Sascha Fleer
//-------------------------------------------------
//-------------------------------------------------


#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>
#include<limits>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<time.h>
#include <stdio.h>

//-------------------------------------------------
// Include library for OPENMP, if supported
//-------------------------------------------------

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using std::max;
using std::min;


//-------------------------------------------------
// Create STRUCTURE for initial conditions of problem
//-------------------------------------------------
struct STATE 
{
	STATE (int s) : side(s) {};
	STATE () {};
	double  temperature;
	double  chempot;
	double 	baryon;
	double  velx;
	double  pressure;
	double  energy;
	double  gamma;
	double  soundvel;
	int side;
};	


//-------------------------------------------------
// To get usefull units of pressure and energy, 
// multiply with fermi^3
// --> New Unit: Mev^3/fm^3
//-------------------------------------------------

double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));


//-------------------------------------------------
// Header with equations of state
//  + Ideal Free Gluon Gas, taken from
// 	  E. Molnar, H. Niemi and D.H. Rischke, 
// 	  The European Physical Journal C 65 (2009) 615
// 	+ QCD, taken from
// 	  Z. Fodor et al., 
// 	  Journal of High Energy Physics 2010 (2010) 77  
//-------------------------------------------------

#include "EOS.h"


//-------------------------------------------------
// Include numerical schemes
//-------------------------------------------------

#include "NUMSCHEME.h"


//-------------------------------------------------
// Include Dissipation
//  + Used Algorithms
//	    - M. Takamoto and S.i. Inutsuka,
//		  Journal of computational physics 230 (2011) 7002
//		- Y. Akamatsu et al., 
//		  arXiv.org (2013) 34, 1302.1665 
//-------------------------------------------------

#include "DISSIPATION.h"


//-------------------------------------------------
// Read in 'input.ini'
//-------------------------------------------------

void readin(std::string * filename, STATE * L1,STATE * L2, STATE * R1, STATE * R2, STATE * PERT, double *STOPTIME, double *DELTAT, int *GRIDPOINTS, double *GRIDLENGTH, double *SHEARCONST, string *SCHEME, string *SLIMITER, string *WHICHEOS, string *VISCOUS); 

bool is_file_exist(const char *fileName);

int main(int argc, char** argv)
{
	//-------------------------------------------------
	// Initialize the initial states
	//-------------------------------------------------

	STATE iLEFT1(0), iRIGHT1(1), iLEFT2(0), iRIGHT2(1);

	//-------------------------------------------------
    	// Initialize pertubation state
	//-------------------------------------------------

	STATE PERT(0);



	//-------------------------------------------------
   	// Initialize equation of state and numerical scheme
	//-------------------------------------------------

	EOS * eofstate;

	NUMSCHEME * nscheme;



	//-------------------------------------------------
    	// Initialization of some systemspecific variables
	//-------------------------------------------------


    // Scheme, Slope-Limiter, which EoS, dissipation?
	string SCHEME, SLIMITER, WHICHEOS, VISCOUS;
	int yesno=0, bulkyes=0;

	// end-time, time-discretization, length of system
	double STOPTIME, DELTAT, GRIDLENGTH;
	STOPTIME=DELTAT=GRIDLENGTH=0;

    // number of gridpoints
	int GRIDPOINTS;
	GRIDPOINTS=0;

    // Matrix for primitive variables in each gridpoint
	double **GVAR;

    // Vector for x-coordinate and spatial discretization in each point
	double *X , *DX;



	//-------------------------------------------------
    	// Initialization of variables needed for 
    	// the computation of dissipation
	//-------------------------------------------------

	// eta/s 
	double SHEARCONST=0;

	// zeta	
	double *BULKCONST;


	double **GVARINITIAL; 
	double *VINITIAL,*VRIEMANN, *TAUETA, *PXX, *BULK;




	//-------------------------------------------------
    	// Variables for tracking computation-time
	//-------------------------------------------------

	time_t now, later;
	int mod;


	//-------------------------------------------------
    	// Test if OPENMP is supported
    	// Print number of usable CPUs 
	//-------------------------------------------------

#ifdef _OPENMP
	printf("OPENMP:\n");
	printf("Number of usable CPUs: %d \n", omp_get_num_procs());
#else	
	{
		printf("OpenMP is not supported. \n");
	}
#endif


	time(&now);


	//-------------------------------------------------
	// Read in parameters from 'input.ini'
	//-------------------------------------------------
    
    if (argc != 3 || argv[1] == argv[2])
    {
       throw invalid_argument("Not enough arguments! \n Please provide the configuration file and the folder to save the results! \n e.g.: \n ./hydro input.ini DATA");
    }



    std::string filename = argv[1]; 
    std::string foldername = argv[2]; 
    std::cout << "Used Configuration File: " << filename << std::endl;
    std::cout << "Used Folder: " << foldername << std::endl;

    if( false == is_file_exist(filename.c_str()))
        throw invalid_argument("Configuration File not existing!"); 

    if( false == is_file_exist(foldername.c_str()))
        throw invalid_argument("Folder for Results not existing!"); 

	readin(&filename, &iLEFT1,&iLEFT2, &iRIGHT1, &iRIGHT2, &PERT, &STOPTIME, &DELTAT, &GRIDPOINTS, &GRIDLENGTH, &SHEARCONST, &SCHEME, &SLIMITER, &WHICHEOS, &VISCOUS); 

	if(WHICHEOS== "QCD")
		eofstate = new QCD;
	else if(WHICHEOS== "GAS")
		eofstate = new GAS;
	else {
		cout << "WRONG INPUT (EOS) !!!" << endl;
		return 0;
	}

	if(VISCOUS== "Y") 
		yesno=1;
	else if(VISCOUS== "BULK")
		bulkyes=yesno=1;
	else if(VISCOUS== "N")
		yesno=0;
	else {
		cout << "WRONG INPUT (VISCOSITY) !!!" << endl;
		return 0;
	}


	//-------------------------------------------------
    	// Fill structures with help of EoS
	//-------------------------------------------------

	iLEFT1.pressure=eofstate->get_pressure(iLEFT1.temperature);
	iLEFT1.energy=eofstate->get_energy(iLEFT1.temperature);
	//Gamma-Factor 1-Dim	
	iLEFT1.gamma=1/sqrt(1-(pow(iLEFT1.velx,2)));
	iLEFT1.baryon=eofstate->get_baryondensity(iLEFT1.temperature, iLEFT1.chempot)*fermi3;
	iLEFT1.soundvel=eofstate->get_soundvel(iLEFT1.temperature);
	iLEFT2.pressure=eofstate->get_pressure(iLEFT2.temperature);
	iLEFT2.energy=eofstate->get_energy(iLEFT2.temperature);
	iLEFT2.gamma=1/sqrt(1-(pow(iLEFT2.velx,2)));
	iLEFT2.baryon=eofstate->get_baryondensity(iLEFT2.temperature, iLEFT2.chempot)*fermi3;
	iLEFT2.soundvel=eofstate->get_soundvel(iLEFT2.temperature);
	iRIGHT1.pressure=eofstate->get_pressure(iRIGHT1.temperature);
	iRIGHT1.energy=eofstate->get_energy(iRIGHT1.temperature);
	iRIGHT1.gamma=1/sqrt(1-(pow(iRIGHT1.velx,2)));
	iRIGHT1.baryon=eofstate->get_baryondensity(iRIGHT1.temperature, iRIGHT1.chempot)*fermi3;
	iRIGHT1.soundvel=eofstate->get_soundvel(iRIGHT1.temperature);
	iRIGHT2.pressure=eofstate->get_pressure(iRIGHT2.temperature);
	iRIGHT2.energy=eofstate->get_energy(iRIGHT2.temperature);
	iRIGHT2.gamma=1/sqrt(1-(pow(iRIGHT2.velx,2)));
	iRIGHT2.baryon=eofstate->get_baryondensity(iRIGHT2.temperature, iRIGHT2.chempot)*fermi3;
	iRIGHT2.soundvel=eofstate->get_soundvel(iRIGHT2.temperature);

	//-------------------------------------------------
	// Pertubation
    	// If no pertubation applied, the gridpoint has
    	// the same values as L2
	//-------------------------------------------------
	PERT.pressure=eofstate->get_pressure(PERT.temperature + iLEFT2.temperature);
	PERT.energy=eofstate->get_energy(PERT.temperature + iLEFT2.temperature);
	PERT.gamma=1/sqrt(1-(pow(PERT.velx + iLEFT2.velx,2)));
	PERT.baryon=eofstate->get_baryondensity(PERT.temperature + iLEFT2.temperature, PERT.chempot)*fermi3;
	PERT.soundvel=eofstate->get_soundvel(PERT.temperature + iLEFT2.temperature);



	//-------------------------------------------------
    	// Initialize the Grid
    	// Gridsize [START-4 : STOP+4]
	//-------------------------------------------------

	X = new double [GRIDPOINTS+8];
	DX = new double [GRIDPOINTS+8];
	GVAR = new double * [GRIDPOINTS+8];
	for(int I=0; I<GRIDPOINTS+8;I++) GVAR[I] = new double [3];

	double ZERO=4;
	double STOP=GRIDPOINTS+4;


	//-------------------------------------------------
    	// Save variables for dissipationalgorithm
    	//  + Startvelocity
    	//  + Velocity - Riemann-Problem
    	//  + Matrix to save all variables temporary
	//-------------------------------------------------

	VINITIAL = new double  [GRIDPOINTS+8];
	VRIEMANN= new double  [GRIDPOINTS+8];
	GVARINITIAL= new double * [GRIDPOINTS+8];
	for(int I=0; I<GRIDPOINTS+8;I++) GVARINITIAL[I] = new double [3];


	//-------------------------------------------------
    	// Relaxation time
	//-------------------------------------------------

	TAUETA= new double  [GRIDPOINTS+8];


	//-------------------------------------------------
    	// Varaibles for shear and bulk viscosity
	//-------------------------------------------------

	PXX= new double  [GRIDPOINTS+8];
	BULK= new double  [GRIDPOINTS+8];
	double ETA=SHEARCONST*pow(10,3)/(fermi3);		

	//-------------------------------------------------
	// Const. zeta for bulk viscosity
	//-------------------------------------------------

	BULKCONST= new double  [GRIDPOINTS+8];


	for(int I=0; I<GRIDPOINTS+8;I++)
	{
		TAUETA[I]=BULK[I]=BULKCONST[I]=PXX[I]=VINITIAL[I]=VRIEMANN[I]=0;
		for(int K=0; K<3; K++)
			GVARINITIAL[I][K]=0;
	}

	//-------------------------------------------------
    	// Variable for Progress bar
	//-------------------------------------------------

	int barWidth=70;

	//-------------------------------------------------
    	// Initialize numerical scheme with initial 
    	// parameters
	//	+ MUSCL --> MUSCL-HANCOCK-Method
	//-------------------------------------------------

    if(SCHEME== "MUSCL")
        nscheme = new MUSCL (eofstate,GRIDPOINTS,GRIDLENGTH,iLEFT1,iLEFT2,iRIGHT1,iRIGHT2, PERT, GVAR,X,DX,SLIMITER);
    else {
        cout << "(SCHEME) WRONG INPUT!!!" << endl;
        return 0;
    }
	


	//-------------------------------------------------
    	// START of main function
	//	+ FN --> Timestep 
	//	+ FTIME --> Complete Time [fm]
	//-------------------------------------------------

	float progress = 0.0;
	int FN=1;

	double FTIME=0;


	//-------------------------------------------------
    	// Open files for contour-plot-data of primitve variables
	//-------------------------------------------------

	fstream c_output_p, c_output_v, c_output_n;

	c_output_p.open(foldername + "/c_output_p.txt", ios::out);
	c_output_v.open(foldername + "/c_output_v.txt", ios::out);
	c_output_n.open(foldername + "/c_output_n.txt", ios::out);


	//-------------------------------------------------
    	// Open files to store data of every time step
    	// Files are saved in "./DATA/" folder, which has to
    	// be creatd manually!
	//-------------------------------------------------

	fstream output;

	output.open(foldername + "/output_0.txt", ios::out);

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


	//-------------------------------------------------
    	// START of main-loop
	//-------------------------------------------------

	while(	FTIME <= STOPTIME)
	{

		//-------------------------------------------------
        	// Create progress bar
		//-------------------------------------------------

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress*100.0) << " %\r"; 
		std::cout.flush();


		//-------------------------------------------------
        	// Initialize boundary conditions
		//-------------------------------------------------

		nscheme->boundary(GVAR);


		//-------------------------------------------------
        	// Compute relaxation time
		//	+ tau=1/T
		//-------------------------------------------------


		for(int I=ZERO; I< STOP; I++)
		{
			double FERMI=0.197*pow(10,3);
			double T=eofstate->get_temperature(GVAR[I][2]);
			TAUETA[I]=1/(T/FERMI);
			if(bulkyes==1)
			{
				double C=1/3-eofstate->get_soundvel(T);
				BULKCONST[I]=15*pow(C,2)*ETA;
			}

            // Check if DELTAT is small enough as DELTAT>TAU
            // If not, decrease it

			if(DELTAT>TAUETA[I] || DELTAT==TAUETA[I]) DELTAT=TAUETA[I]*0.1;

			VINITIAL[I]=GVAR[I][1];
			for(int K=0; K<3;K++)
				GVARINITIAL[I][K]=GVAR[I][K];
		}


		//-------------------------------------------------
        	// Temporal evolution of primitive variables GVAR
        	// under influence of the chosen numerical scheme
		//-------------------------------------------------

		nscheme->evaluate(GVAR,DELTAT);


		//-------------------------------------------------
        	// Save primitive variables for dissipation
        	// algorithm
		//-------------------------------------------------
		for(int I=ZERO; I< STOP; I++)
			VRIEMANN[I]=GVAR[I][1];


		//-------------------------------------------------
		//DISSIPATION
		//-------------------------------------------------
		
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
		//-------------------------------------------------
		// Compute nex time-step (CFL-Condition, 
        	// CFL-Number:0.4)
		//-------------------------------------------------

		FTIME=FTIME+nscheme->get_timestep();

		//-------------------------------------------------
		// Save data to "DATA/output_*.txt"
		//-------------------------------------------------

		fstream output;

		output.open(foldername + "/output_"+std::to_string(FN)+".txt", ios::out);

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


	//-------------------------------------------------
    	// Close files for contour plots
	//-------------------------------------------------

	c_output_p.close();
	c_output_v.close();
	c_output_n.close();
	cout << endl;


	//-------------------------------------------------
	// Print computation time
	//-------------------------------------------------

	time(&later);
	mod= ((int) difftime(later,now))% 60;

	cout << "Computation Time: " << (difftime(later,now)- (double) mod)/60.<< " minutes & " << mod  << " seconds!" << endl; 
	return 0;
}


//-------------------------------------------------
// Read in parameters from 'input.ini'
//-------------------------------------------------

void readin(std::string * filename, STATE * L1 ,STATE * L2, STATE * R1, STATE * R2, STATE * PERT, double *STOPTIME, double *DELTAT, int *GRIDPOINTS, double *GRIDLENGTH, double *SHEARCONST, string *SCHEME, string *SLIMITER, string *WHICHEOS, string *VISCOUS) 
{

	std::string templine; 
	ifstream sourcedata;
	sourcedata.open(*filename, ios_base::in);

	if (!sourcedata)
		cout << "Could not open data!" << endl;
	else 
	{                                       
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		L1->temperature=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		L1->velx=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		L2->temperature=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		L2->velx=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		PERT->temperature=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		PERT->velx=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		R1->temperature=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		R1->velx=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		R2->temperature=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		R2->velx=atof(templine.c_str());
		sourcedata.ignore(numeric_limits<streamsize>::max(), '\n');
		getline(sourcedata,templine);
		L1->chempot=L2->chempot=R1->chempot=R2->chempot=atof(templine.c_str());
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

bool is_file_exist(const char *fileName)
    {
        std::ifstream infile(fileName);
        return infile.good();
    }
