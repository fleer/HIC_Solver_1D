#include <stdlib.h>

//-------------------------------------------------
// Include Class for Approximative Riemann solver,
//  based on
//  Y. Akamatsu et al., arXiv.org (2013) 34, 1302.1665
//-------------------------------------------------

#include "RIEMANN.h"


//-------------------------------------------------
// Include class with numerical schemes
//-------------------------------------------------

class NUMSCHEME
{
public:
	virtual void evaluate(double **VAR, double TIMESTEP)=0;
	virtual void boundary(double **VAR)=0;
	virtual double get_timestep()=0;
};


//-------------------------------------------------
// 1+1 Dimensional MUSCL-HANCOCK-Methode
//-------------------------------------------------

#include "MUSCL.h"
