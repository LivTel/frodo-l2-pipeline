

// Includes

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>

// *********************************************************************

// Entry point

int main(int argc, char *argv[]) {

	if (argc != 5) {
/*
		printf("\n");
		printf("### frprocess_spectra_TSV ###\n\n");
		printf(" _____ ____   ___  ____   ___  ____\n");
		printf("|  ___|  _ \\ / _ \\|  _ \\ / _ \\/ ___| _ __   ___  ___\n");
		printf("| |_  | |_) | | | | | | | | | \\___ \\| '_ \\ / _ \\/ __|\n");
		printf("|  _| |  _ <| |_| | |_| | |_| |___) | |_) |  __/ (__ \n");
		printf("|_|   |_| \\_\\\\___/|____/ \\___/|____/| .__/ \\___|\\___|\n");
		printf("				    |_|\n");
		printf("\n");
		printf("Author: Rob Barnsley (LJMU)\n");
		printf("Email: rmb@astro.livjm.ac.uk\n");
		printf("\n");
		printf("Usage: frprocess [input_file, output_file, HDU, is_wavelength_calibrated?]\n");
		printf("\n");
		printf("Synopsis:\n");
		printf("\n");
		printf("Accesses FRODOSpec _2.fits images, takes the an extension and outputs a wavelength/flux TSV file.\n");
		printf("\n");
		printf("Example:\n");
		printf("\n");
		printf("frprocess r_e..._2.fits output.fits 6 1\n");
		printf("\n");
*/
		return 1;

	}

	// ***********************************************************************
	// Set input file parameters and open it

	fitsfile *inf_fptr;

	int inf_status = 0;
	int inf_bitpix, inf_naxis;
	long inf_naxes[2] = {1,1};

	int HDU = atoi(argv[3]);

	if(!fits_open_file(&inf_fptr, argv[1], READONLY, &inf_status)) { 
		
		// move to correct HDU

		fits_movabs_hdu(inf_fptr, HDU, IMAGE_HDU, &inf_status);

		if(!fits_get_img_param(inf_fptr, 2, &inf_bitpix, &inf_naxis, inf_naxes, &inf_status)) {

			printf("\n");
			printf("Successfully opened input frame.");
			printf("\n");
			printf("\n");
			printf("(INPUT FRAME) Bits per pixel: %d\n", inf_bitpix);
			printf("(INPUT FRAME) Number of axes: %d\n", inf_naxis);
			printf("(INPUT FRAME) First axis dimension: %d\n", inf_naxes[0]);
			printf("(INPUT FRAME) Second axis dimension: %d\n", inf_naxes[1]);
			printf("\n");

			//if (inf_naxis != 2) { printf("The number of axes must be 2. Quitting.\n\n"); return 0; }

		} else { fits_report_error(stdout, inf_status); printf("\n"); return 0; }

	} else { fits_report_error(stdout, inf_status); printf("\n"); return 0; }

	// ***********************************************************************
	// Set parameters to be used when reading files eg. first pixel to be read etc..

	int anynul;
	long fpixel[2] = {1,1};
	long nelements = inf_naxes[0];

	// ***********************************************************************
	// Allocate memory and check

	double *inf_pixels;	
	inf_pixels = (double *) malloc(nelements * sizeof(double));

        if (!inf_pixels) {

              printf("Memory allocation error. Quitting. \n\n");
              return(0);

        } else { printf("Memory allocation OK. \n\n"); }

	printf("Executing program..");

	// ***********************************************************************
	// Find wavelength related .FITS keys from source file if flag set

	int ii;

	int is_wavelength_calibrated = atoi(argv[4]);

	double wavelengths[nelements];

	if (is_wavelength_calibrated == 1) {

		double CDELT1value;
		double CRVAL1value;
		double CRPIX1value;

		char CDELT1comment[81];
		char CRVAL1comment[81];
		char CRPIX1comment[81];

		fits_read_key_dbl(inf_fptr, "CDELT1", &CDELT1value, CDELT1comment, &inf_status); 
		fits_read_key_dbl(inf_fptr, "CRVAL1", &CRVAL1value, CRVAL1comment, &inf_status); 
		fits_read_key_dbl(inf_fptr, "CRPIX1", &CRPIX1value, CRPIX1comment, &inf_status); 

		printf("\n\nCDELT1\t%f\n", CDELT1value);
		printf("CRVAL1\t%f\n", CRVAL1value);
		printf("CRPIX1\t%f\n\n", CRPIX1value);

		for (ii=0; ii<nelements; ii++) {

			wavelengths[ii] = CRVAL1value + (CRPIX1value * CDELT1value) + (ii*CDELT1value);

		}
	
	} else {

		for (ii=0; ii<nelements; ii++) {

			wavelengths[ii] = ii;

		}
	
	}

	// ***********************************************************************
	// Read row and write the values into a 1D array at the correct place

	if(!fits_read_pix(inf_fptr, TDOUBLE, fpixel, nelements, NULL, inf_pixels, NULL, &inf_status)) {

	} else { fits_report_error(stdout, inf_status); printf("\n"); return 0; }


	// ***********************************************************************
	// Get some fits keys relating to frame

	char DATEvalue[200];
	char UTSTARTvalue[200];
	char OBJECTvalue[200];
	char GRATIDvalue[200];
	char RAvalue[200];
	char DECvalue[200];

	char DATEcomment[200];
	char UTSTARTcomment[200];
	char OBJECTcomment[200];
	char GRATIDcomment[200];
	char RAcomment[200];
	char DECcomment[200];

	fits_read_key_str(inf_fptr, "DATE", DATEvalue, DATEcomment, &inf_status); 
	fits_read_key_str(inf_fptr, "UTSTART", UTSTARTvalue, UTSTARTcomment, &inf_status); 
	fits_read_key_str(inf_fptr, "OBJECT", OBJECTvalue, OBJECTcomment, &inf_status); 
	fits_read_key_str(inf_fptr, "GRATID", GRATIDvalue, GRATIDcomment, &inf_status); 
	fits_read_key_str(inf_fptr, "RA", RAvalue, RAcomment, &inf_status); 
	fits_read_key_str(inf_fptr, "DEC", DECvalue, DECcomment, &inf_status);

	// Set output file

  	 FILE *outputfile;
  	 outputfile = fopen(argv[2], "w");

	// and write

	/*

	fprintf(outputfile, "####\n\n");

	fprintf(outputfile, "Date Start:\t%s\n", DATEvalue);
	fprintf(outputfile, "Time Start:\t%s\n", UTSTARTvalue);
	fprintf(outputfile, "Object:\t%s\n", OBJECTvalue);
	fprintf(outputfile, "Config:\t%s\n", GRATIDvalue);
	fprintf(outputfile, "RA:\t%s\n", RAvalue);
	fprintf(outputfile, "DEC:\t%s\n\n", DECvalue);

	*/

	// ***********************************************************************
	// Write to file

	for (ii=0; ii<nelements; ii++) {

		fprintf(outputfile, "%f\t%f\n", wavelengths[ii], inf_pixels[ii]);

	}

	// ***********************************************************************
	// Close files
	
	printf("Closing files tidily.\n\n");

	fits_close_file(inf_fptr, &inf_status);
	fclose(outputfile);

}
