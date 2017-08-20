//-------
// include APPROXIMATE RIEMANN SOLVER (Akamatsu et. al 2013) 
//-------

#include "riemann.h"
#include <stdlib.h>

class NUMSCHEME
{
public:
	virtual void evaluate(double **VAR, double TIMESTEP)=0;
	virtual void boundary(double **VAR)=0;
	virtual double get_timestep()=0;
};


//-------
// include 1-DIM PPM-METHOD for time evolution (adapted version of Marti & Mueller 1996) 
//-------

#include "ppm.h"


//-------
// include 1-DIM MUSCL-HANCOCK SCHEME for time evolution 
//-------

#include "MUSCL.h"
