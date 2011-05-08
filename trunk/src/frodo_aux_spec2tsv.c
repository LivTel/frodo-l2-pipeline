/************************************************************************

 File:				frodo_aux_spec2tsv.c
 Last Modified Date:     	08/05/11

************************************************************************/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "frodo_error_handling.h"
#include "frodo_functions.h"
#include "frodo_config.h"
#include "frodo_aux_spec2tsv.h"

int main(int argc, char *argv[]) {

	if (argc != 5) {

		if(populate_env_variable(FAS_BLURB_FILE, "L2_FAS_BLURB_FILE")) {

			printf("Failed to get L2_FAS_BLURB_FILE environment variable\n\n");

		} else {

			print_file(FAS_BLURB_FILE);

		}

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f			= strdup(argv[1]);
		int HDU				= strtol(argv[2], NULL, 0);
		int is_wavelength_calibrated	= strtol(argv[3], NULL, 0);
		char *output_f			= strdup(argv[4]);

		// ***********************************************************************
		// Open input file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2, input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes [2] = {1,1};

		if(!fits_open_file(&input_f_ptr, input_f, IMG_READ_ACCURACY, &input_f_status)) {

			if(!fits_movabs_hdu(input_f_ptr, HDU+1, IMAGE_HDU, &input_f_status)) {		// HDU+1 so the index is 0 not 1

				if(!fits_get_img_param(input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status)) {

					if (input_f_naxes[1] != 1) {

						printf("Size of y axis != 1. This is not a spectral frame\n\n");
						fits_report_error(stdout, input_f_status); 

						free(input_f); 
						free(output_f);

						if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

						return 1; 

					}

					// printf("%d\n", input_f_naxes[0]);	// DEBUG

				} else { 

					printf("Failed to get image parameters of input frame\n\n");
					fits_report_error(stdout, input_f_status); 

					free(input_f); 
					free(output_f);

					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

			} else {

				printf("Failed to move to HDU %d\n\n", HDU);
				fits_report_error(stdout, input_f_status);

				free(input_f);
				free(output_f);

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		} else { 

			printf("Failed to open input FITS file (ARG 1)\n\n");
			fits_report_error(stdout, input_f_status); 

			free(input_f);
			free(output_f);

			return 1; 

		}	

		// ***********************************************************************
		// Set the range limits

		int cut_x [2] = {1, input_f_naxes[0]};
		int cut_y [2] = {1, input_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from input fits file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from input fits file (ARG 1)

		double input_f_pixels [nxelements];

		// ***********************************************************************
		// Get input fits file (ARG 1) values and store in 1D array

		if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

		} else { 

			fits_report_error(stdout, input_f_status);

			free(input_f);
			free(output_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

		}

		// ***********************************************************************
		// Get wavelength calibration header keys if applicable

		int ii;

		double values [nxelements];

		if (is_wavelength_calibrated == TRUE) {

			double key_value_CDELT1;
			double key_value_CRVAL1;
			double key_value_CRPIX1;

			char key_comment_CDELT1 [81];
			char key_comment_CRVAL1 [81];
			char key_comment_CRPIX1 [81];

			if(fits_read_key_dbl(input_f_ptr, "CDELT1", &key_value_CDELT1, key_comment_CDELT1, &input_f_status)) {

				printf("Failed to get CDELT1 key\n\n");
				fits_report_error(stdout, input_f_status); 

				free(input_f);
				free(output_f);

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 
		
			} 

			if(fits_read_key_dbl(input_f_ptr, "CRVAL1", &key_value_CRVAL1, key_comment_CRVAL1, &input_f_status)) {

				printf("Failed to get CRVAL1 key\n\n");
				fits_report_error(stdout, input_f_status);

				free(input_f);
				free(output_f); 

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

			if(fits_read_key_dbl(input_f_ptr, "CRPIX1", &key_value_CRPIX1, key_comment_CRPIX1, &input_f_status)) {

				printf("Failed to get CRPIX1 key\n\n");
				fits_report_error(stdout, input_f_status);

				free(input_f);
				free(output_f);

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

			for (ii=0; ii<nxelements; ii++) {

				values[ii] = key_value_CRVAL1 + (key_value_CRPIX1 * key_value_CDELT1) + (ii*key_value_CDELT1);

			}
	
		} else {

			for (ii=0; ii<nxelements; ii++) {

				values[ii] = ii;

			}
	
		}

		// ***********************************************************************
		// Create output file and print results

		FILE *outputfile;
		outputfile = fopen(output_f, FILE_WRITE_ACCESS);

		for (ii=0; ii<nxelements; ii++) {

			fprintf(outputfile, "%f\t%f\n", values[ii], input_f_pixels[ii]);

		}

		// ***********************************************************************
		// Clean up heap memory

		free(input_f);
		free(output_f); 

		// ***********************************************************************
		// Close input file (ARG 1) and output file (ARG 4)

		if(fclose(outputfile)) {

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			fits_report_error (stdout, input_f_status); 

			return 1; 

	    	}

		return 0;

	}

}

