/************************************************************************

 File:				frodo_red_extract_simple.c
 Last Modified Date:     	07/02/11

************************************************************************/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "frodo_error_handling.h"
#include "frodo_functions.h"

// *********************************************************************

int main ( int argc, char *argv [] ) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if ( argc != 5 ) {

		if(populate_env_variable(HEADER_FILE, "L2_HEADER_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", 1, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		if(populate_env_variable(FRES_BLURB_FILE, "L2_FRES_BLURB_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", 2, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		print_file(HEADER_FILE);
		print_file(FRES_BLURB_FILE);

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -1, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else { 

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f			= strdup(argv[1]);
		int half_aperture_num_pix	= atoi(argv[2]);
		bool flip_disp_axis		= atoi(argv[3]);
		char *output_f			= strdup(argv[4]);

		// ***********************************************************************
		// Is there EXPLICIT cross-talk? i.e. Is the extraction aperture width >
		// the average fibre profile width?

		if (half_aperture_num_pix > (FREXTRACT_VAR_FIBRE_PROFILE_WIDTH - 1) / 2) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", 3, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}
	
		// ***********************************************************************
		// Open input file (ARG 1), get parameters and perform any data 
		// format checks

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2, input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes [2] = {1,1};

		if(!fits_open_file(&input_f_ptr, input_f, IMG_READ_ACCURACY, &input_f_status)) {

			if(!populate_img_parameters(input_f, input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status, "INPUT FRAME")) {

				if (input_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -2, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -3, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -4, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status); 

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
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from input fits file (ARG 1)

		double input_f_pixels [nxelements];

		// ***********************************************************************
		// Open [FRTRACE_OUTPUTF_TRACES_FILE] file
	
		FILE *traces_file;
	
		if (!check_file_exists(FRTRACE_OUTPUTF_TRACES_FILE)) { 

			traces_file = fopen(FRTRACE_OUTPUTF_TRACES_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -5, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}

		// ***********************************************************************
		// Find some [FRTRACE_OUTPUTF_TRACES_FILE] file details

		char input_string [500];

		bool find_polynomialorder_comment = FALSE;

		int polynomial_order;	

		char search_string_1 [20] = "# Polynomial Order:\0";	// this is the comment to be found from the [FRTRACE_OUTPUTF_TRACES_FILE] file

		while(!feof(traces_file)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, traces_file);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &polynomial_order);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_polynomialorder_comment = 1;
				break;


			} 

		}

		if (find_polynomialorder_comment != 1) {	// error check - didn't find the comment in the [FRTRACE_OUTPUTF_TRACES_FILE] file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -6, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}

		// ***********************************************************************
		// Rewind and extract coefficients from [FRTRACE_OUTPUTF_TRACES_FILE] file 

		rewind(traces_file);

		int token_index;	// this variable will hold which token we're dealing with
		int coeff_index;	// this variable will hold which coefficient we're dealing with
		int this_fibre;	
		double this_coeff;
		double this_chisquared;
	
		char *token;

		double coeffs [FIBRES][polynomial_order+1];
		memset(coeffs, 0, sizeof(double)*FIBRES*(polynomial_order+1));

		while(!feof(traces_file)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, traces_file);

			token_index = 0;
			coeff_index = 0;

			if (atoi(&input_string[0]) > 0) { 		// check the line begins with a positive number

				// ***********************************************************************
				// String tokenisation loop: 
				//
				// 1. init calls strtok() loading the function with input_string
				// 2. terminate when token is null
				// 3. we keep assigning tokens of input_string to token until termination by calling strtok with a NULL first argument
				// 
				// n.b. searching for tab or newline separators ('\t' and '\n')

				for (token=strtok(input_string, "\t\n"); token !=NULL; token = strtok(NULL, "\t\n")) {	

					if (token_index == 0 ) {							// fibre token

						this_fibre = atoi(token);

					} else if ((token_index >= 1) && (token_index <= polynomial_order+1)) { 	// coeff token

						this_coeff = atof(token);
						// printf("%d\t%d\t%e\n", this_fibre, coeff_index, this_coeff);		// DEBUG
						coeffs[this_fibre-1][coeff_index] = this_coeff;
						coeff_index++;

					} else if (token_index == polynomial_order+2) {					// chisquared token

						this_chisquared = atof(token);

					}

					token_index++;

				}

			}

		}	

		// ***********************************************************************
		// Standard aperture extraction of flux using coefficients from
		// [FRTRACE_OUTPUTF_TRACES_FILE] file

		int ii, jj;

		double x;
		int x_int;
		int y_int;

		double input_frame_values [FIBRES][nyelements];
		memset(input_frame_values, 0, sizeof(double)*nyelements*FIBRES);

		// double total_flux = 0;	// DEBUG
		// int count = 0;		// DEBUG

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(input_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

				y_int = fpixel[1];

				for (ii=0; ii<FIBRES; ii++) {

					x = 0;

					// ***********************************************************************
					// Find tracing centroid from polynomial coefficients
	
					for (jj=0; jj<polynomial_order+1; jj++) {
					
						x = (coeffs[ii][jj]*(pow(y_int,jj))) + x;
	
					}

					x_int = rint(x) - INDEXING_CORRECTION;

					// ***********************************************************************
					// Does [x_int] violate the img boundaries?

					if ((x_int > nxelements) || (x_int <= 0)) {

						write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -7, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
						fits_report_error(stdout, input_f_status); 

						return 1;

					}

					// ***********************************************************************
					// Extract flux within pixel aperture

					for (jj=x_int-half_aperture_num_pix; jj<=x_int+half_aperture_num_pix; jj++) {

						input_frame_values[ii][y_int-1] += input_f_pixels[jj];
						// total_flux += input_f_pixels[jj];	// DEBUG

					}

					//printf("%f\t%f\t%f\t%f\n", x, input_f_pixels[x_int], input_f_pixels[x_int-1], input_f_pixels[x_int+1]);		// DEBUG

					// if ((input_f_pixels[x_int] > input_f_pixels[x_int-1]) && (input_f_pixels[x_int] > input_f_pixels[x_int+1])) count++;	// DEBUG

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -8, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				return 1; 

			}

		}

		// printf("%f\n", total_flux); // DEBUG
		// printf("%d\n", count); // DEBUG

		// ***********************************************************************
		// Set output frame parameters

		fitsfile *output_f_ptr;
	
		int output_f_status = 0;
		long output_f_naxes [2] = {nyelements,FIBRES};
	
		long output_f_fpixel = 1;

		// ***********************************************************************
		// Create [output_frame_values] array to hold the output data in the 
		// correct format

		int kk;

		double output_frame_values [FIBRES*nyelements];
		memset(output_frame_values, 0, sizeof(double)*FIBRES*nyelements);

		double intermediate_frame_values [nyelements];

		for (ii=0; ii<FIBRES; ii++) {

			memset(intermediate_frame_values, 0, sizeof(double)*nyelements);
	
			jj = ii * nyelements;
	
			for (kk=0; kk<nyelements; kk++) {
	
				intermediate_frame_values[kk] = input_frame_values[ii][kk];
	
			}

			if (flip_disp_axis == TRUE) {

				flip_array_dbl(intermediate_frame_values, nyelements);

			}	

			for (kk=0; kk<nyelements; kk++) {

				output_frame_values[jj] = intermediate_frame_values[kk];
				jj++;

			}
		
		}

		// ***********************************************************************
		// Create and write [output_frame_values] to output file (ARG 4)
	
		if (!fits_create_file(&output_f_ptr, output_f, &output_f_status)) {
	
			if (!fits_create_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, output_f_naxes, &output_f_status)) {

				if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, FIBRES * nyelements, output_frame_values, &output_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -9, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, output_f_status); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -10, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -11, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Close input file (ARG 1) and output file (ARG 4)

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -12, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, input_f_status); 

			return 1; 

	    	}

		if(fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -13, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Free arrays on heap

		free(input_f);
		free(output_f);

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", SUCCESS_FLAG, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}
