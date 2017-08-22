#include <stdlib.h>

//-------------------------------------------------
// Einbinden der Header mit dem approximaitven Riemann-Solver aus
// Y. Akamatsu et al., arXiv.org (2013) 34, 1302.1665
//-------------------------------------------------

#include "RIEMANN.h"


//-------------------------------------------------
// Header mit Klasse als Template zum einfachen implementieren von
// numerischen Schmes
//-------------------------------------------------

class NUMSCHEME
{
public:
	virtual void evaluate(double **VAR, double TIMESTEP)=0;
	virtual void boundary(double **VAR)=0;
	virtual double get_timestep()=0;
};


//-------------------------------------------------
// Einbinden der Header fuer die 1+1 Dimensionale MUSCL-HANCOCK-Methode
//-------------------------------------------------

#include "MUSCL.h"


//-------------------------------------------------
// Einbinden der Header fuer die 1+1 Dimensionale PPM-Methode
// (Wird in der Masterarbeit nicht verwendet)
//-------------------------------------------------

#include "PPM.h"
