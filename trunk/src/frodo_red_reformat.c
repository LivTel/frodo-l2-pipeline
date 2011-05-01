/************************************************************************

 File:				frodo_red_reformat.c
 Last Modified Date:     	01/05/11

************************************************************************/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_sort_double.h>
#include "frodo_error_handling.h"
#include "frodo_functions.h"
#include "frodo_config.h"
#include "frodo_red_reformat.h"
#include "frodo_red_rebin.h"

// *********************************************************************/

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 5) {

		if(populate_env_variable(FRRF_BLURB_FILE, "L2_FRRF_BLURB_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", 2, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		print_file(FRRF_BLURB_FILE);

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -1, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f		= strdup(argv[1]);
		char *headers_f		= strdup(argv[2]);
		char *operation		= strdup(argv[3]);
		char *output_f		= strdup(argv[4]);

		// ***********************************************************************
		// Check operation choice is recognised

		if (strcmp(operation, "L1_IMAGE") && strcmp(operation, "RSS_NONSS") && strcmp(operation, "RSS_SS") && strcmp(operation, "CUBE_NONSS") && strcmp(operation, "CUBE_SS") && strcmp(operation, "SPEC_NONSS") && strcmp(operation, "SPEC_SS") && strcmp(operation, "COLCUBE_NONSS")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -25, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			return 1;

		}

		// ***********************************************************************
		// Check input file (ARG 1) exists and get parameters if so

		bool missing_input_file = FALSE;

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2;
		int input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes[2] = {1,1};

		if(check_file_exists(input_f)) {	// file doesn't exist	
	
			missing_input_file = TRUE;

		} else {

			if(!fits_open_file(&input_f_ptr, input_f, READONLY, &input_f_status)) {

				if(!populate_img_parameters(input_f, input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status, "INPUT FRAME")) {} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -6, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					return 1; 

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -7, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Open the headers file (ARG 2)

		fitsfile *headers_f_ptr;

		int headers_f_status = 0;

		if(!fits_open_file(&headers_f_ptr, headers_f, READONLY, &headers_f_status)) {} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -4, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, headers_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Create the output file (ARG 4) if it doesn't already exist, otherwise
		// open

		fitsfile *output_f_ptr;

		int output_f_status = 0;

		if(check_file_exists(output_f)) {	// file doesn't exist

			if (!fits_create_file(&output_f_ptr, output_f, &output_f_status)) {

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -2, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 
			
			}

		} else {

			if(!fits_open_file(&output_f_ptr, output_f, READWRITE, &output_f_status)) {

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -3, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Print operation specifics to screen

		int hdunum;

		if(!fits_get_num_hdus(output_f_ptr, &hdunum, &output_f_status)) {} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -5, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			return 1; 

		}

		printf("\nOperation:\t%s\n", operation);
		printf("Current HDU:\t%d\n", hdunum+1);

		// SPECIFIC DATA HANDLING
		// ***********************************************************************

		int ii, jj;

		// ***********************************************************************
		// Write blank extension to output file (ARG 4) if input file (ARG 1) is
		// missing

		if (missing_input_file == TRUE) {

			// CREATE BLANK EXTENSION IN OUTPUT FILE 			   (BLANK)
			// ***********************************************************************
			// 1.	Set output file (ARG 4) parameters
			long output_f_naxes[1] = {1};
			int output_f_bitpix = 16, output_f_naxis = 1;

			// 2.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -7, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

			RETURN_FLAG = 2;
			
		} else if (!strcmp(operation, "L1_IMAGE") || !strcmp(operation, "RSS_NONSS") || !strcmp(operation, "RSS_SS")) { 

			// COPY DATA FROM INPUT FILE TO OUTPUT FILE 			    (COPY)
			// ***********************************************************************
			// 1.	Set output file (ARG 4) parameters

			long output_f_naxes[2] = {input_f_naxes[0], input_f_naxes[1]};
			int output_f_bitpix = input_f_bitpix, output_f_naxis = input_f_naxis;

			// 2.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -8, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

			// 3.	Copy values across from primary extension of input file (ARG 1)

			if(!fits_copy_data(input_f_ptr, output_f_ptr, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -9, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

		} else if (!strcmp(operation, "CUBE_NONSS") || !strcmp(operation, "CUBE_SS")) {

			// REFORMAT DATA FROM INPUT FILE TO OUTPUT FILE 		(DATACUBE)
			// ***********************************************************************
			// 1.	Set the range limits using input fits file (ARG 1)

			int cut_x [2] = {1, input_f_naxes[0]};
			int cut_y [2] = {1, input_f_naxes[1]};

			// 2.	Set parameters used when reading data from input fits file (ARG 1)

			long fpixel [2] = {cut_x[0], cut_y[0]};
			long nxelements = (cut_x[1] - cut_x[0]) + 1;
			long nyelements = (cut_y[1] - cut_y[0]) + 1;

			// 3.	Set output file (ARG 4) parameters

			long output_f_naxes[3] = {FRREFORMAT_VAR_IFU_DIM_X, FRREFORMAT_VAR_IFU_DIM_Y, nxelements};
			int output_f_bitpix = input_f_bitpix, output_f_naxis = 3;
			long output_f_fpixel = 1;

			double output_frame_values [FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y*nxelements];
	     		memset(output_frame_values, 0, sizeof(double)*FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y*nxelements);

			// 4.	Create array to store pixel values from input fits file (ARG 1)

			double input_f_pixels [nxelements];

			// 5.	Read row of input frame (ARG 1)

			int this_fibre_index;

			double input_frame_values [nyelements][nxelements];
			memset(input_frame_values, 0, sizeof(double)*nyelements*nxelements);

			for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

				this_fibre_index = fpixel[1] - 1;

				memset(input_f_pixels, 0, sizeof(double)*nxelements);

				if(!fits_read_pix(input_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

					for (ii=0; ii<nxelements; ii++) {

						input_frame_values[this_fibre_index][ii] = input_f_pixels[ii];

					}

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -26, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					return 1; 

				}

			}

			// 6.	Rearrange input frame values [input_frame_values] into datacube format

			int output_frame_values_index = 0;

			for (ii=0; ii<nxelements; ii++) {

				for (jj=0; jj<nyelements; jj++) {

					output_frame_values[output_frame_values_index] = input_frame_values[FRREFORMAT_VAR_IFU_FIBRE_FORMAT[jj]][ii];
					output_frame_values_index++;

				}

			}

			// 7.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -27, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

			// 8.	Write datacube values [output_frame_values] to output file (ARG 4)

			if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y*nxelements, output_frame_values, &output_f_status)) {} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -29, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

		} else if (!strcmp(operation, "SPEC_NONSS") || !strcmp(operation, "SPEC_SS")) {

			// REFORMAT DATA FROM INPUT FILE TO OUTPUT FILE 		(SPECTRUM)
			// ***********************************************************************
			// 1.	Set the range limits using input fits file (ARG 1)

			int cut_x [2] = {1, input_f_naxes[0]};
			int cut_y [2] = {1, input_f_naxes[1]};

			// 2.	Set parameters used when reading data from input fits file (ARG 1)

			long fpixel [2] = {cut_x[0], cut_y[0]};
			long nxelements = (cut_x[1] - cut_x[0]) + 1;
			long nyelements = (cut_y[1] - cut_y[0]) + 1;

			// 3.	Set output file (ARG 4) parameters

			long output_f_naxes[2] = {nxelements, 1};
			int output_f_bitpix = input_f_bitpix, output_f_naxis = 2;
			long output_f_fpixel = 1;

			double output_frame_values [nxelements];
	     		memset(output_frame_values, 0, sizeof(double)*nxelements);

			// 4.	Create array to store pixel values from input fits file (ARG 1)

			double input_f_pixels [nxelements];
			double smoothed_input_f_pixels [nxelements];

			// 5.	Read row of input frame (ARG 1)

			int this_fibre_index;

			double total_fluxes [nyelements];
			memset(total_fluxes, 0, sizeof(double)*nyelements);	

			double input_frame_values [nyelements][nxelements];
			memset(input_frame_values, 0, sizeof(double)*nyelements*nxelements);

			for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

				this_fibre_index = fpixel[1] - 1;

				memset(input_f_pixels, 0, sizeof(double)*nxelements);
				memset(smoothed_input_f_pixels, 0, sizeof(double)*nxelements);

				if(!fits_read_pix(input_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

					// 6.	Apply a median smoothing filter to eliminate CR from brightest fibre selection

					median_filter(input_f_pixels, smoothed_input_f_pixels, nxelements, FRREFORMAT_VAR_SPECTRUM_MEDIAN_HALF_FILTER_SIZE);

					for (ii=0; ii<nxelements; ii++) {

						input_frame_values[this_fibre_index][ii] = input_f_pixels[ii];

						// 7.	Sum total flux and store

						total_fluxes[this_fibre_index] += smoothed_input_f_pixels[ii];

					}

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -26, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					return 1; 

				}

			}

			// 8.	Find indexes of the brightest [FRREFORMAT_VAR_SPECTRUM_NUM_FIBRES] fibres

			size_t brightest_fibre_indexes[FRREFORMAT_VAR_SPECTRUM_NUM_FIBRES];

			gsl_sort_largest_index(brightest_fibre_indexes, FRREFORMAT_VAR_SPECTRUM_NUM_FIBRES, total_fluxes, 1, nyelements);

			// 9.	Add flux from each of the brightest rows

			for (ii=0; ii<nxelements; ii++) {

				for (jj=0; jj<FRREFORMAT_VAR_SPECTRUM_NUM_FIBRES; jj++) {

					output_frame_values[ii] += input_frame_values[brightest_fibre_indexes[jj]][ii];

				}

			}

			// 10.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -27, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

			// 11.	Write spectrum values [output_frame_values] to output file (ARG 4)

			if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, nxelements, output_frame_values, &output_f_status)) {} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -29, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

		} else if (!strcmp(operation, "COLCUBE_NONSS")) {

			// REFORMAT DATA FROM INPUT FILE TO OUTPUT FILE 	      (IFU BUNDLE)
			// ***********************************************************************
			// 1.	Set the range limits using input fits file (ARG 1)

			int cut_x [2] = {1, input_f_naxes[0]};
			int cut_y [2] = {1, input_f_naxes[1]};

			// 2.	Set parameters used when reading data from input fits file (ARG 1)

			long fpixel [2] = {cut_x[0], cut_y[0]};
			long nxelements = (cut_x[1] - cut_x[0]) + 1;
			long nyelements = (cut_y[1] - cut_y[0]) + 1;

			// 3.	Set output file (ARG 4) parameters

			long output_f_naxes[3] = {FRREFORMAT_VAR_IFU_DIM_X, FRREFORMAT_VAR_IFU_DIM_Y};
			int output_f_bitpix = input_f_bitpix, output_f_naxis = 2;
			long output_f_fpixel = 1;

			double output_frame_values [FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y];
	     		memset(output_frame_values, 0, sizeof(double)*FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y);

			// 4.	Create array to store pixel values from input fits file (ARG 1)

			double input_f_pixels [nxelements];

			// 5.	Read row of input frame (ARG 1)

			int this_fibre_index;

			double total_fluxes [nyelements];
			memset(total_fluxes, 0, sizeof(double)*nyelements);

			for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

				this_fibre_index = fpixel[1] - 1;

				memset(input_f_pixels, 0, sizeof(double)*nxelements);

				if(!fits_read_pix(input_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

					for (ii=0; ii<nxelements; ii++) {

						// 6.	Sum total flux and store

						total_fluxes[this_fibre_index] += input_f_pixels[ii];

					}

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -26, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					return 1; 

				}

			}

			// 7.	Rearrange input frame values [input_frame_values] into IFU bundle format

			int output_frame_values_index = 0;

			for (jj=0; jj<nyelements; jj++) {

				output_frame_values[output_frame_values_index] = total_fluxes[FRREFORMAT_VAR_IFU_FIBRE_FORMAT[jj]];
				output_frame_values_index++;

			}

			// 8.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -27, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}

			// 9.	Write datacube values [output_frame_values] to output file (ARG 4)

			if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, FRREFORMAT_VAR_IFU_DIM_X*FRREFORMAT_VAR_IFU_DIM_Y, output_frame_values, &output_f_status)){} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -29, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				return 1; 

			}


		}

		// ***********************************************************************
		// Add EXTNAME keyword

		if (!fits_write_key(output_f_ptr, TSTRING, "EXTNAME", operation, NULL, &output_f_status)) {} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -32, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Copy SCALE keywords from input file (ARG 1) to output file (ARG 4) if
		// they exist

		if (missing_input_file == FALSE) {

			if(populate_env_variable(FITS_KEYS_TO_OMIT, "L2_FITS_KEYS_TO_OMIT_FILE")) {	// this file contains a list of keywords to be omitted from each extension

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", 3, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

			} else {

				char card [FLEN_CARD];
				bool found_key;

				int nkeys;

				if(!fits_get_hdrspace(input_f_ptr, &nkeys, NULL, &input_f_status)) {} else {

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -10, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, headers_f_status); 

					return 1; 

				} 

				for (ii=0; ii<nkeys; ii++) {

					found_key = FALSE;

					if(!fits_read_record(input_f_ptr, ii, card, &input_f_status)) {} else {

						write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -11, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
						fits_report_error(stdout, input_f_status); 

						return 1; 

					}

					if (fits_get_keyclass(card) == TYP_SCAL_KEY) {				// copy only scale keywords (BZERO/BSCALE)

						check_key_to_omit(FITS_KEYS_TO_OMIT, card, operation, &found_key);

						if (found_key == TRUE) {

							// printf("%d\n", TRUE);	// DEBUG

						} else {

							// printf("%d\t%s\n", fits_get_keyclass(card), card);
							fits_write_record(output_f_ptr, card, &output_f_status);

						}

					}

				}

				// ***********************************************************************
				// Copy headers from headers file (ARG 2) to output file (ARG 4)

				if(!fits_get_hdrspace(headers_f_ptr, &nkeys, NULL, &headers_f_status)) {} else {

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -12, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, headers_f_status); 

					return 1; 

				} 

				for (ii=0; ii<nkeys; ii++) {

				// printf("%d\t%d\n", ii, nkeys);	// DEBUG
		
					found_key = FALSE;

					if(!fits_read_record(headers_f_ptr, ii, card, &headers_f_status)) {} else {

						write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -13, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
						fits_report_error(stdout, headers_f_status); 

						return 1; 

					}

					if (fits_get_keyclass(card) >= TYP_WCS_KEY) {

						check_key_to_omit(FITS_KEYS_TO_OMIT, card, operation, &found_key);

						if (found_key == TRUE) {

							// printf("%d\n", TRUE);	// DEBUG

						} else {

							// printf("%d\t%s\n", fits_get_keyclass(card), card);

							fits_write_record(output_f_ptr, card, &output_f_status);

						}

					} // else { printf("%d\t%s\n", fits_get_keyclass(card), card); }	// DEBUG

				// printf("%d\n", &nkeys);	// DEBUG
	
				}

			}

		}

		// SPECIFIC HEADERS HANDLING
		// ***********************************************************************

		if (!strcmp(operation, "L1_IMAGE")) {	

		} else if (!strcmp(operation, "RSS_NONSS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "RSS_CALIBRATION", 2, &output_f_status);

		} else if (!strcmp(operation, "CUBE_NONSS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "CUBE_CALIBRATION", 2, &output_f_status);

		} else if (!strcmp(operation, "RSS_SS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "RSS_CALIBRATION", 2, &output_f_status);

		} else if (!strcmp(operation, "CUBE_SS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "CUBE_CALIBRATION", 2, &output_f_status);

		} else if (!strcmp(operation, "SPEC_NONSS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "SPEC_CALIBRATION", 2, &output_f_status);	

		} else if (!strcmp(operation, "SPEC_SS")) {

			write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "SPEC_CALIBRATION", 2, &output_f_status);

		} else if (!strcmp(operation, "COLCUBE_NONSS")) {

		}

		// Write L2 status keywords
		// ***********************************************************************

		write_error_codes_file_to_header(ERROR_CODES_FILE, output_f_ptr, &output_f_status);

		// ***********************************************************************
		// Write checksums

	  	fits_write_chksum(output_f_ptr, &output_f_status);

		// ***********************************************************************
		// Close input file (ARG 1), headers file (ARG 3) and output file (ARG 4)

		if (missing_input_file == FALSE) {

			if (fits_close_file(input_f_ptr, &input_f_status)) { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -14, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		}

		if (fits_close_file(headers_f_ptr, &headers_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -15, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, headers_f_status); 

			return 1; 

		    }

		if (fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -16, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Clean up heap memory

		// ***********************************************************************
		// Return success

	}

}
