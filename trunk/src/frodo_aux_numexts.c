/************************************************************************

 File:				frodo_aux_numexts.c
 Last Modified Date:     	05/05/11

************************************************************************/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include "frodo_error_handling.h"
#include "frodo_functions.h"
#include "frodo_config.h"
#include "frodo_aux_numexts.h"

int main(int argc, char *argv[]) {

	if (argc != 2) {

		if(populate_env_variable(FAN_BLURB_FILE, "L2_FAN_BLURB_FILE")) {

			printf("Failed to get L2_FAN_BLURB_FILE environment variable\n\n");

		} else {

			print_file(FAN_BLURB_FILE);

		}

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f		= strdup(argv[1]);

		// ***********************************************************************
		// Open input file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2, input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes [2] = {1,1};

		int num_HDUs;

		if(!fits_open_file(&input_f_ptr, input_f, IMG_READ_ACCURACY, &input_f_status)) {

			if (!fits_get_num_hdus(input_f_ptr, &num_HDUs, &input_f_status)) {	

			} else {

				printf("Failed to get number of HDUs in input frame\n\n");
				fits_report_error(stdout, input_f_status);

				free(input_f); 

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 	

			}

		} else { 

			printf("Failed to open input FITS file (ARG 1)\n\n");
			fits_report_error(stdout, input_f_status); 

			free(input_f); 

			return 1; 

		}

		int ii;

		int filled_HDUs = num_HDUs;

		for (ii=1; ii<=num_HDUs; ii++) { 

			if(!fits_movabs_hdu(input_f_ptr, ii, IMAGE_HDU, &input_f_status)) {

				if(!fits_get_img_param(input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status)) {

					if (input_f_naxis == 1) { 	// this is blank

						filled_HDUs--;

					}

					// printf("%d\n", input_f_naxis);	// DEBUG

				} else { 

					printf("Failed to get image parameters of input frame\n\n");
					fits_report_error(stdout, input_f_status); 

					free(input_f); 

					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

			} else {

				printf("Failed to move to HDU %d\n\n", ii);
				fits_report_error(stdout, input_f_status); 

				free(input_f); 

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		}

		printf("%d", filled_HDUs);

		// ***********************************************************************
		// Clean up heap memory

		free(input_f); 

		// ***********************************************************************
		// Close input file (ARG 1)

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			fits_report_error (stdout, input_f_status); 

			return 1; 

	    	}

		return 0;

	}

}

