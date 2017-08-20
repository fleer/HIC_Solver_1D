class EOS 
{
public:

virtual double get_pressure(double TEMP) =0;
virtual double get_energy(double TEMP) =0;
virtual double get_temperature(double P) =0; 
virtual double get_baryondensity(double TEMP, double NB) =0;
virtual double get_susceptibility(double TEMP) =0;
virtual double get_deriv_susceptibility(double TEMP) =0;
virtual double get_soundvel(double TEMP) =0;
double recovery(double * cons);

private:
double fprimitive( double * cons, double p);
double dfprimitive( double * cons, double p);
};

double EOS::fprimitive( double * cons, double p)
{
	double T=get_temperature(p);
	double e=get_energy(T);
	double g=1/(1-(pow(cons[1],2)/pow(cons[2]+p,2)));
//	cout << "G: " << g << endl;
//cout << "F: " << ((e+p)*g)-cons[2]-p << endl;
	return ((e+p)*g)-cons[2]-p;
}
				
double EOS::dfprimitive(double * cons, double p)
{
	double T=get_temperature(p);
	double e=get_energy(T);
	double g=1/(1-(pow(cons[1],2)/pow(cons[2]+p,2)));
	double nB=cons[0]/sqrt(g);
	//cout << "NB " << nB << endl;
double c=1/get_soundvel(T);
double chi=get_susceptibility(T);
double dchi=get_deriv_susceptibility(T);

double A=(nB/chi)*(1+((T*dchi)/chi)-c);

return (c+1)*g-1+(((A*nB)-2*(e+p))*(g/(cons[2]+p))*(g-1));
}

double EOS::recovery( double * cons)
{
#define RECOVERYMAX 1000
#define RECOVERYTOL 10e-8
double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
	/*
	//Generate Flux for Testing
	double cons[3]; 

	cons[0]=primitive.gamma*primitive.baryon;
	cons[1]=(primitive.energy+primitive.pressure)*pow(primitive.gamma,2)*primitive.velx;
	cons[2]=(cons[1]/primitive.velx)-primitive.pressure;

	cout << "D:" << cons[0] << endl;
	cout << "m:" << cons[1] << endl;
	cout << "E:" << cons[2] << endl;
*/	
	double IP,SP=0;
	double dif=0;
	double abb=0;
	double abb2=0;
	//Guess pressure:  Need Good apprixmation
	
	//IP=fabs(fabs(cons[1])-cons[2]-cons[0]);
	IP=get_pressure(400);
//cout << "GUESS PRESSURE:"<< IP*fermi3 << endl;
//IP=eos->get_pressure(300);

	for(int w=0; w<= RECOVERYMAX; w++)
	{
		//cout << w << endl;
		dif=fprimitive(cons,IP )/ dfprimitive(cons, IP);
	//	cout << "f: " << fprimitive(cons, IP) << "\t" << dfprimitive(cons, IP) << endl;


	//	cout << "DIF:" << dif << endl;
		SP=fabs(IP-dif);
		//abb=SP-IP;
		abb=(SP-IP)/(0.5*(IP+SP));
	//cout << "COND:" << fabs(abb) << endl;
		if(fabs(abb)<= RECOVERYTOL) 
		{

			//cout << "PRIM P:" << IP*fermi3 << endl;
			return IP; 
		}
		else 
		{	
			//Additional Routine if TOLERANZ is too big for some Pressure Values on grid points	


	/*		if(abb2<fabs(abb)+RECOVERYTOL && abb2>fabs(abb)-RECOVERYTOL)
			{
				//cout << "OUT OF BOUNDS!" << endl;
				return IP;
			}
*/
			if(w>100 && fabs(fabs(abb2)-fabs(abb)) < 1) return IP;

				abb2=abb;
				
				IP=SP;
		//		cout << "ITERATION P:" << SP << endl;
//cout << w << "\t" << get_temperature(IP) << "\t" << abb << endl;
		}

	}	
cout << "ITERATION P:" << SP*fermi3 << endl;
	cout << "TO MANY ITERATIONS IN FINDING P!!" << endl;
	exit(0);
	return 0;
}
///-----------------
//Code for QCD-Eos of Borsnyi et. al (2010)
//----------------
//definition of parameter structure for integration in QCD
#define h0 0.1396
#define h1 -0.18
#define h2 0.035
#define f0 2.76
#define f1 6.79
#define f2 -5.29
#define g1 -0.47
#define g2 1.04


class QCD : public EOS 
{
public:

static double f(double x, void * params)
	{

		double f= (exp ((-h1/(x/200.))-(h2/pow((x/200.),2.)))*(h0+(f0*(tanh ((f1*(x/200.))+f2)+1.)/(1+(g1*(x/200.))+(g2*pow((x/200.),2.))))))/x;
		return f;
	}
double get_pressure(double TEMP);
double get_energy(double TEMP);
double get_temperature(double P); 
double get_baryondensity(double TEMP, double NB);
double get_susceptibility(double TEMP);
double get_deriv_susceptibility(double TEMP);
double get_soundvel(double TEMP);
private:

//---------
//Definition of Anomaly I
//---------

double anomaly(double x)
	{
		return (exp ((-h1/(x/200.))-(h2/pow((x/200.),2.)))*(h0+(f0*(tanh ((f1*(x/200.))+f2)+1.)/(1+(g1*(x/200.))+(g2*pow((x/200.),2.))))));
	}

};




//----------
//compute pressure from temperature, integration with gsl_integration_qags
//---------
double QCD::get_pressure(double TEMP)
{
	double error;
	double result=0;
	double alpha=1.0;
		gsl_function F;
	 void* params_ptr = &alpha;
		F.params=params_ptr;
		F.function = f;
	//	cout << TEMP << endl;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	

		gsl_integration_qags(&F ,0, TEMP,0,1e-12,10000,w,&result,&error);
				//*fermi3 so be cautious with the units... not sure if ok!!
		gsl_integration_workspace_free (w);
	return result*pow(TEMP,4);

}

//---------
//compute pressure for given p
//get p through simple Newton-Rhapson-Algorithm with shiftet root to p
//Be car eful with maximal T... maybe some modifictations would be useful
//P'(T)=(e+P)/T
//---------
double QCD::get_temperature(double P)
{
double T=1000; //MAXIMAL Range for Temperature
// Maybe T is maxmimal T in initial conditions
double IT=0;
double abb2=0;
double abb=0;
#define ITERATIONS 100
#define TOL 1e-10

for(int i=0; i<=ITERATIONS; i++)
{
	IT=T-((get_pressure(T)-P)/((get_energy(T)+get_pressure(T))/T));
	abb=IT-T;
	if(fabs(abb)<=TOL) return T;
	else
	{

			//Additional Routine if TOLERANZ is too big for some Temperatures Values on grid points	


			if(abb2<fabs(abb)+TOL && abb2>fabs(abb)-TOL)
			{
				//cout << "OUT OF BOUNDS!" << endl;
				return T;
			}

			else
			{
				abb2=abb;
				T=IT;
				//cout << T << endl;
				//cout << "ITERATION P:" << SP << endl;
			}	
	}
}	
cout << "TO MANY ITERATIONS IN FINDING T!!" << endl;
cout << "P: " << P << endl;
return 0;
}


//-------
//computation of energydensity
//-------
double QCD::get_energy(double TEMP)
{
return (anomaly(TEMP)*pow(TEMP,4))+(3*get_pressure(TEMP));
}

//-------
//computation of baryondensity through susceptibility
//-------
double QCD::get_baryondensity(double TEMP, double NB)
{

double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return (a*pow(TEMP,2)*(1+(tanh((TEMP-T0)/(deltaT)))))*NB*fermi3;
}

double QCD::get_susceptibility(double TEMP)
{
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return (a*pow(TEMP,2)*(1+(tanh((TEMP-T0)/(deltaT)))));
}
//------
//computation of dervative of susceptibility
//for primitive recovery in PPM algorithmus
//------
double QCD::get_deriv_susceptibility(double TEMP)
{
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return a*TEMP*(2*(1+(tanh((TEMP-T0)/(deltaT))))+(TEMP*(1/(pow(cosh((TEMP-T0)/(deltaT)),2)*deltaT))));
}
//------
//computation of squared shock velocity
//------
double QCD::get_soundvel(double TEMP)
{
	//cS IS SQUARED!!!
//-------------	
	//BAD APPROXIMATION
	//return get_pressure(TEMP)/get_energy(TEMP);
//-------------	
//acosh = 1/cosh??
	// Derivative of anmomaly I	
	double A=1+(g1*TEMP*200.)+(g2*pow(TEMP*200.,2));
	double B=(f0*(tanh((f1*TEMP*200.)+f2)+1.))/A;
	double C=((f0*f1*200.)/pow(cosh((f1*200.*TEMP)+f2),2))/A;
	double DI=exp(-(h1/(TEMP*200.))-(h2/pow(TEMP*200.,2)))*((4*pow(TEMP,3)*(h0+B))+(pow(TEMP,4)*((((h1/(pow(TEMP,2)*200.))+(2*h2/(pow(TEMP,3)*pow(200.,2))))*(h0+B))+(C-(B*((g1*200.)+(2*g2*TEMP*pow(200.,2)))/A)))));
	/*double DI=exp(-((h2*pow(200,2))/pow(TEMP,2)) - (h1*200)/TEMP)*pow(TEMP,4)*((f0*(f1*200)*pow(acosh(f2 + (f1*200)*TEMP),2))/(1 + (g1*200)*TEMP + (g2*pow(200,2))*pow(TEMP,2)) - (f0*((g1*200) + 2*(g2*pow(200,2))*TEMP)*(1 + tanh(f2 + (f1*200)*TEMP)))/pow(1 + (g1*200)*TEMP + (g2*pow(200,2))*pow(TEMP,2),2)) + 4*exp(-((h2*pow(200,2)/pow(TEMP,2)) - (h1*200)/TEMP)*pow(TEMP,3)*(h0 + (f0*(1 + tanh(f2 + (f1*200)*TEMP)))/(1 + (g1*200)*TEMP + (g2*pow(200,2))*pow(TEMP,2))) + exp(-(h2/pow(TEMP,2)) - (h1*200)/TEMP)*((2*h2)/pow(TEMP,3) + (h1*200)/pow(TEMP,2))*pow(TEMP,4)*(h0 + (f0*(1 + tanh(f2 + (f1*200)*TEMP)))/(1 + (g1*200)*TEMP + (g2*pow(200,2))*pow(TEMP,2))));*/
 	double DP=(get_energy(TEMP)+get_pressure(TEMP))/TEMP;
     double DDP=(DI+(3*DP));
         return DP/DDP;
}
//-----------------
//Code for GAS-EOS from MOLNAR et al. (2009)
//----------------



class GAS : public EOS 
{
public:

double get_pressure(double TEMP);
double get_energy(double TEMP);
double get_temperature(double P); 
double get_baryondensity(double TEMP, double NB);
double get_susceptibility(double TEMP);
double get_deriv_susceptibility(double TEMP);
double get_soundvel(double TEMP);
};

//----------
//compute pressure from temperature, integration with gsl_integration_qags
//---------
double GAS::get_pressure(double TEMP)
{
		return ((16.*pow(TEMP,4.))/pow(M_PI,2.));
}

//---------
//compute pressure for given p
//get p through simple Newton-Rhapson-Algorithm with shiftet root to p
//Be careful with maximal T... maybe some modifictations would be useful
//P'(T)=(e+P)/T
//---------
double GAS::get_temperature(double P)
{
return pow(((P*M_PI*M_PI)/16.),0.25);
}


//-------
//computation of energydensity
//-------
double GAS::get_energy(double TEMP)
{
return ((48.*pow(TEMP,4.))/pow(M_PI,2.));
}

//-------
//computation of baryondensity through susceptibility
//-------
double GAS::get_baryondensity(double TEMP, double NB)
{
double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
		double a=0.15; //From QCD-EOS

		return a*pow(TEMP,2)*NB*fermi3;
}
double GAS::get_susceptibility(double TEMP)
{
		double a=0.15; //From QCD-EOS

		return a*pow(TEMP,2);
}

//------
//computation of dervative of susceptibility
//for primitive recovery in PPM algorithmus
//------
double GAS::get_deriv_susceptibility(double TEMP)
{
	double a=0.15;

	return 2*a*TEMP;
}

//------
//computation of squared shock velocity
//------
double GAS::get_soundvel(double TEMP)
{
return 1./3.;
}
