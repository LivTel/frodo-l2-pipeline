/************************************************************************

 File:				frodo_red_arcfit.h
 Last Modified Date:     	07/03/11

************************************************************************/

char FRA_BLURB_FILE [200];

int FRARCFIT_VAR_REF_FIBRE_INDEX_CC		= 0;				// this is the ARRAY INDEX of the fibre used as a reference in cross correlation analysis (e.g. fibre 1 = 0)

int FRARCFIT_VAR_POLYORDER_LO			= 2;
int FRARCFIT_VAR_POLYORDER_HI			= 10;

char FRARCFIT_VAR_ACCURACY_COEFFS [10]		= "%.10e";
char FRARCFIT_VAR_ACCURACY_CHISQ [10]		= "%.2f";

double FRARCFIT_VAR_CHISQUARED_MIN		= 0.1;
double FRARCFIT_VAR_CHISQUARED_MAX		= 5;

char FRARCFIT_OUTPUTF_WAVFITS_FILE [100]	= "frarcfit_wavfits.dat";

