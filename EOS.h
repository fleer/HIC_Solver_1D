//-------------------------------------------------
// Header mit Zustandsgleichungen 
// 	+ Ideales Gluongas aus
// 	  E. Molnar, H. Niemi and D.H. Rischke, 
// 	  The European Physical Journal C 65 (2009) 615
// 	+ QCD
// 	  Z. Fodor et al., 
// 	  Journal of High Energy Physics 2010 (2010) 77  
//-------------------------------------------------


//-------------------------------------------------
// Klasse EOS als Template zur einfachen Implementierung weiterer
// Zustandsgleichungen
//-------------------------------------------------

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

		//-------------------------------------------------
		// Funktionen fuer den Algorithmus zur Wiederherstellung
		// der primitiven Variablen
		//-------------------------------------------------

		double fprimitive( double * cons, double p);
		double dfprimitive( double * cons, double p);
};

double EOS::fprimitive( double * cons, double p)
{
	double T=get_temperature(p);
	double e=get_energy(T);
	double g=1/(1-(pow(cons[1],2)/pow(cons[2]+p,2)));
	return ((e+p)*g)-cons[2]-p;
}
				
double EOS::dfprimitive(double * cons, double p)
{
	double T=get_temperature(p);
	double e=get_energy(T);
	double g=1/(1-(pow(cons[1],2)/pow(cons[2]+p,2)));
	double nB=cons[0]/sqrt(g);
	double c=1/get_soundvel(T);
	double chi=get_susceptibility(T);
	double dchi=get_deriv_susceptibility(T);

	double A=(nB/chi)*(1+((T*dchi)/chi)-c);

	return (c+1)*g-1+(((A*nB)-2*(e+p))*(g/(cons[2]+p))*(g-1));
}

double EOS::recovery( double * cons)
{
	int RECOVERYMAX=1000;
	double  RECOVERYTOL=10e-8;
	double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
	double IP,SP=0;
	double dif=0;
	double abb=0;
	double abb2=0;

	//-------------------------------------------------
	//
	// Newton-Rhapson-Verfahren
	//
	// Raten eines Anfangsdrucks fuer das Newton-Rhapson-Verfahren
	// Hier: Druck fuer Temperatur T=400 MeV
	// Sehr unelegant, fuehrt aber im Gegensatz zur eleganteren
	// Variante
	// 		IP=fabs(fabs(cons[1])-cons[2]-cons[0]);
	// nie zu Problemen
	//-------------------------------------------------

	IP=get_pressure(400);

	for(int w=0; w<= RECOVERYMAX; w++)
	{
		dif=fprimitive(cons,IP )/ dfprimitive(cons, IP);
		SP=fabs(IP-dif);
		abb=(SP-IP)/(0.5*(IP+SP));
		if(fabs(abb)<= RECOVERYTOL) 
		{
			return IP; 
		}
		else 
		{	
			if(w>100 && fabs(fabs(abb2)-fabs(abb)) < 1) return IP;
			abb2=abb;
			IP=SP;
		}
	}	
	cout << "ITERATION P:" << SP*fermi3 << endl;
	cout << "TO MANY ITERATIONS IN FINDING P!!" << endl;
	exit(0);
	return 0;
}

///-----------------
// QCD-Zustandsgelichung aus 
// Z. Fodor et al., Journal of High Energy Physics 2010 (2010) 77
//----------------

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


//-------------------------------------------------
// Definition der Anomalie I
//-------------------------------------------------

double anomaly(double x)
	{
		return (exp ((-h1/(x/200.))-(h2/pow((x/200.),2.)))*(h0+(f0*(tanh ((f1*(x/200.))+f2)+1.)/(1+(g1*(x/200.))+(g2*pow((x/200.),2.))))));
	}

};


//-------------------------------------------------
// Berechnung des Drucks aus der Temperatur
// Integration mit gsl_integration_qags aus der mathgl-Bibliothek
//-------------------------------------------------

double QCD::get_pressure(double TEMP)
{
	double error;
	double result=0;
	double alpha=1.0;
	gsl_function F;
	void* params_ptr = &alpha;
	F.params=params_ptr;
	F.function = f;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags(&F ,0, TEMP,0,1e-12,10000,w,&result,&error);
	gsl_integration_workspace_free (w);
	return result*pow(TEMP,4);

}


//-------------------------------------------------
// Berechnung der Temperatur aus gegebenem Druck
// Anfangstemperatur fuer das Newton-Rhapson-Verfahren ist T=1 GeV
// Nicht sehr schoen, aber es funktioniert
//-------------------------------------------------

double QCD::get_temperature(double P)
{
	double T=1000; 
	double IT=0;
	double abb2=0;
	double abb=0;
	int ITERATIONS=100;
	double TOL=1e-10;

	for(int i=0; i<=ITERATIONS; i++)
	{
		IT=T-((get_pressure(T)-P)/((get_energy(T)+get_pressure(T))/T));
		abb=IT-T;
		if(fabs(abb)<=TOL) return T;
		else
		{

			// Zusaetzliche Routine um Konvergenzprobleme abzufangen
			if(abb2<fabs(abb)+TOL && abb2>fabs(abb)-TOL)
			{
				return T;
			}

			else
			{
				abb2=abb;
				T=IT;
			}	
		}
	}	
	cout << "TO MANY ITERATIONS IN FINDING T!!" << endl;
	cout << "P: " << P << endl;
	return 0;
}


//-------------------------------------------------
// Berechnung der Energiedichte
//-------------------------------------------------

double QCD::get_energy(double TEMP)
{
return (anomaly(TEMP)*pow(TEMP,4))+(3*get_pressure(TEMP));
}


//-------------------------------------------------
// Berechnung der Beyondichte mit Hilfe der Suszeptibilitaet
//-------------------------------------------------

double QCD::get_baryondensity(double TEMP, double NB)
{

double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return (a*pow(TEMP,2)*(1+(tanh((TEMP-T0)/(deltaT)))))*NB*fermi3;
}


//-------------------------------------------------
// Berechnung der Suszeptibilitaet
//-------------------------------------------------

double QCD::get_susceptibility(double TEMP)
{
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return (a*pow(TEMP,2)*(1+(tanh((TEMP-T0)/(deltaT)))));
}


//-------------------------------------------------
// Berechnung der Ableitung der Suszeptibilitaet
//-------------------------------------------------

double QCD::get_deriv_susceptibility(double TEMP)
{
		double a=0.15;
		double T0=167;
		double deltaT=60;
		return a*TEMP*(2*(1+(tanh((TEMP-T0)/(deltaT))))+(TEMP*(1/(pow(cosh((TEMP-T0)/(deltaT)),2)*deltaT))));
}


//-------------------------------------------------
// Berechnung des quadrates der Schallgeschwindigkeit
//-------------------------------------------------

double QCD::get_soundvel(double TEMP)
{
	double A=1+(g1*TEMP*200.)+(g2*pow(TEMP*200.,2));
	double B=(f0*(tanh((f1*TEMP*200.)+f2)+1.))/A;
	double C=((f0*f1*200.)/pow(cosh((f1*200.*TEMP)+f2),2))/A;
	double DI=exp(-(h1/(TEMP*200.))-(h2/pow(TEMP*200.,2)))*((4*pow(TEMP,3)*(h0+B))+(pow(TEMP,4)*((((h1/(pow(TEMP,2)*200.))+(2*h2/(pow(TEMP,3)*pow(200.,2))))*(h0+B))+(C-(B*((g1*200.)+(2*g2*TEMP*pow(200.,2)))/A)))));
	double DP=(get_energy(TEMP)+get_pressure(TEMP))/TEMP;
	double DDP=(DI+(3*DP));
	return DP/DDP;
}


//-----------------
// Zustandsgelichung eines freien Gluon-Gases
// E. Molnar, H. Niemi and D.H. Rischke, 
// The European Physical Journal C 65 (2009) 615
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


//-------------------------------------------------
// Berechnung des Drucks aus der Temperatur
//-------------------------------------------------

double GAS::get_pressure(double TEMP)
{
		return ((16.*pow(TEMP,4.))/pow(M_PI,2.));
}


//-------------------------------------------------
// Berechnung der Temperatur aus gegebenem Druck
//-------------------------------------------------

double GAS::get_temperature(double P)
{
	return pow(((P*M_PI*M_PI)/16.),0.25);
}


//-------------------------------------------------
// Berechnung der Energiedichte
//-------------------------------------------------

double GAS::get_energy(double TEMP)
{
	return ((48.*pow(TEMP,4.))/pow(M_PI,2.));
}


//-------------------------------------------------
// Berechnung der Beyondichte mit Hilfe der Suszeptibilitaet
//-------------------------------------------------

double GAS::get_baryondensity(double TEMP, double NB)
{
	double fermi3= 1/(pow(0.197,3.)*pow(10.,9.));
	
	// Wert fuer auebernommen von QCD-Zustandsgleichung
	double a=0.15; 

	return a*pow(TEMP,2)*NB*fermi3;
}


//-------------------------------------------------
// Berechnung der Ableitung der Suszeptibilitaet
//-------------------------------------------------

double GAS::get_susceptibility(double TEMP)
{
	// Wert fuer auebernommen von QCD-Zustandsgleichung
	double a=0.15;

	return a*pow(TEMP,2);
}


//-------------------------------------------------
// Berechnung der Ableitung der Suszeptibilitaet
//-------------------------------------------------

double GAS::get_deriv_susceptibility(double TEMP)
{
	double a=0.15;

	return 2*a*TEMP;
}


//-------------------------------------------------
// Berechnung des quadrates der Schallgeschwindigkeit
//-------------------------------------------------

double GAS::get_soundvel(double TEMP)
{
return 1./3.;
}
