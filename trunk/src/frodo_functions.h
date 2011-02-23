/************************************************************************

 File:				frodo_functions.h
 Last Modified Date:     	31/01/11

************************************************************************/

#include <stdbool.h>

int calc_least_sq_fit (int, int, double [], double [], double [], double *);
int calculate_cross_correlation_offset (double [], double [], int, int, int *, double *);
int check_file_exists (char []);
int find_centroid_parabolic (double [], int [], int, double [], bool);
int find_peaks(int, double [], int [], int *, int, int, int, int, bool);
int find_peaks_contiguous(int, int, double **, int **, int *, int, int, int, int, int, bool);
int find_time (char []);
int flip_array_dbl(double [], int);
int lsearch_int(int [], int, int);
int populate_img_parameters(char [], fitsfile *, int, int *, int *, long [], int *, char []);
int populate_env_variable(char [], char []);
int print_file(char *);


bool INDEXING_CORRECTION			= TRUE;

char HEADER_FILE [100];

int SUCCESS_FLAG				= 0;

char FILE_WRITE_ACCESS [2] 			= "w";			// r readonly; w overwrite; a+ append

char REF_ERROR_CODES_FILE [100];					// variable to hold location of error codes reference file
char ERROR_CODES_FILE [100] 			= "error_codes";	// location of file to write error codes to
char ADDITIONAL_KEYS_FILE [100] 		= "additional_keys";	// location of file to write additional keys to

char ERROR_CODES_INITIAL_FILE_WRITE_ACCESS [2]	= "w";			// r readonly; w overwrite; a+ append
char ERROR_CODES_FILE_WRITE_ACCESS [2]	 	= "a+";			// r readonly; w overwrite; a+ append

// FRFIND

char FRFS_BLURB_FILE [200];

char FRFIND_OUTPUTF_PEAKS_FILE [100] 		= "frfind_peaks.dat";

int IMG_READ_ACCURACY				= 82;

// FRCLEAN

char FRFSC_BLURB_FILE [200];

char FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE [100] 	= "frclean_cleaned.dat";

// FRTRACE

char FRT_BLURB_FILE [200];

char FRTRACE_OUTPUTF_TRACES_FILE [100] 		= "frtrace_traces.dat";

int FRTRACE_VAR_POLYORDER_LO			= 2;
int FRTRACE_VAR_POLYORDER_HI			= 10;

char FRTRACE_VAR_ACCURACY_COEFFS [10]		= "%.10e";
char FRTRACE_VAR_ACCURACY_CHISQ [10]		= "%.2f";

double FRTRACE_VAR_CHISQUARED_MIN		= 0.1;
double FRTRACE_VAR_CHISQUARED_MAX		= 10;

// FREXTRACT

char FRES_BLURB_FILE [200];

int INTERMEDIATE_IMG_ACCURACY[2]		= {-64, 82};			// Definition follows types from http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node20.html {create_img, write_img}

int FREXTRACT_VAR_FIBRE_PROFILE_WIDTH 		= 7;				// in pixels

// FRARCFIT

char FRA_BLURB_FILE [200];

char FRARCFIT_OUTPUTF_WAVFITS_FILE [100]	= "frarcfit_wavfits.dat";

int FRARCFIT_VAR_REF_FIBRE_INDEX_CC		= 0;				// this is the ARRAY INDEX of the fibre used as a reference in cross correlation analysis (e.g. fibre 1 = 0)

int FRARCFIT_VAR_POLYORDER_LO			= 2;
int FRARCFIT_VAR_POLYORDER_HI			= 10;

char FRARCFIT_VAR_ACCURACY_COEFFS [10]		= "%.10e";
char FRARCFIT_VAR_ACCURACY_CHISQ [10]		= "%.2f";

double FRARCFIT_VAR_CHISQUARED_MIN		= 0.1;
double FRARCFIT_VAR_CHISQUARED_MAX		= 5;

// FRCORRECT

char FRCT_BLURB_FILE [200];

// IFU PARAMETERS

int FIBRES 					= 144;







