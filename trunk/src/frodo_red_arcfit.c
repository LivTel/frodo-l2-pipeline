/************************************************************************

 File:				frodo_red_arcfit.c
 Last Modified Date:     	05/05/11

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
#include "frodo_red_arcfit.h"

#include <gsl/gsl_statistics_int.h>

// *********************************************************************/

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 19) {

		if(populate_env_variable(FRA_BLURB_FILE, "L2_FRA_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FRA_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -1, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {
			
		// ***********************************************************************
		// Redefine routine input parameters
	
		char *ext_arc_f			= strdup(argv[1]);
		char *ext_target_f		= strdup(argv[2]);
		char *ext_cont_f		= strdup(argv[3]);
		int max_cc_delay		= strtol(argv[4], NULL, 0);
		int min_dist			= strtol(argv[5], NULL, 0);
		int half_aperture_num_pix	= strtol(argv[6], NULL, 0);
		int derivative_tol		= strtol(argv[7], NULL, 0);
		int derivative_tol_ref_px	= strtol(argv[8], NULL, 0);
		int pix_tolerance		= strtol(argv[9], NULL, 0);
		int min_contiguous_lines	= strtol(argv[10], NULL, 0);
		char *arc_line_list_filename	= strdup(argv[11]);
		int max_pix_diff		= strtol(argv[12], NULL, 0);
		int min_matched_lines		= strtol(argv[13], NULL, 0);
		double max_av_wavelength_diff	= strtol(argv[14], NULL, 0);
		int fit_order			= strtol(argv[15], NULL, 0);
		char *cc_ext_arc_f		= strdup(argv[16]);
		char *cc_ext_target_f		= strdup(argv[17]);
		char *cc_ext_cont_f		= strdup(argv[18]);

		// ***********************************************************************
		// Open ext arc file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *ext_arc_f_ptr;

		int ext_arc_f_maxdim = 2, ext_arc_f_status = 0, ext_arc_f_bitpix, ext_arc_f_naxis;
		long ext_arc_f_naxes [2] = {1,1};

		if(!fits_open_file(&ext_arc_f_ptr, ext_arc_f, READONLY, &ext_arc_f_status)) {

			if(!populate_img_parameters(ext_arc_f, ext_arc_f_ptr, ext_arc_f_maxdim, &ext_arc_f_bitpix, &ext_arc_f_naxis, ext_arc_f_naxes, &ext_arc_f_status, "ARC FRAME")) {

				if (ext_arc_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -2, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -3, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_arc_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -4, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, ext_arc_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			return 1; 

		}

		// ***********************************************************************
		// Open ext target file (ARG 2), get parameters and perform any data 
		// format checks 

		fitsfile *ext_target_f_ptr;

		int ext_target_f_maxdim = 2;
		int ext_target_f_status = 0, ext_target_f_bitpix, ext_target_f_naxis;
		long ext_target_f_naxes [2] = {1,1};

		if(!fits_open_file(&ext_target_f_ptr, ext_target_f, READONLY, &ext_target_f_status)) {

			if(!populate_img_parameters(ext_target_f, ext_target_f_ptr, ext_target_f_maxdim, &ext_target_f_bitpix, &ext_target_f_naxis, ext_target_f_naxes, &ext_target_f_status, "TARGET FRAME")) {

				if (ext_target_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -5, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
					if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -6, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_target_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -7, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, ext_target_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Open extracted continuum file (ARG 3), get parameters and perform any 
		// data format checks 

		fitsfile *ext_cont_f_ptr;

		int ext_cont_f_maxdim = 2;
		int ext_cont_f_status = 0, ext_cont_f_bitpix, ext_cont_f_naxis;
		long ext_cont_f_naxes [2] = {1,1};

		if(!fits_open_file(&ext_cont_f_ptr, ext_cont_f, READONLY, &ext_cont_f_status)) {

			if(!populate_img_parameters(ext_cont_f, ext_cont_f_ptr, ext_cont_f_maxdim, &ext_cont_f_bitpix, &ext_cont_f_naxis, ext_cont_f_naxes, &ext_cont_f_status, "CONTINUUM FRAME")) {

				if (ext_cont_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -8, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
					if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
					if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -9, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_cont_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -10, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, ext_cont_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Check consistency of arc/target/continuum fits files (ARGS 1, 2 and 3)
	
		printf("\nConsistency check");
		printf("\n-----------------\n");

		printf("\nBits per pixel:\t\t");

		if (ext_arc_f_bitpix != ext_target_f_bitpix && ext_target_f_bitpix != ext_cont_f_bitpix) { 	// if a = b and b = c then a must = c

			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -11, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		} else { 

			printf("OK\n"); 

		} 

		printf("Number of axes:\t\t");

		if (ext_arc_f_naxis != ext_target_f_naxis && ext_target_f_naxis != ext_cont_f_naxis) {	
	
			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -12, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		} else { 

			printf("OK\n"); 

		} 
	
		printf("First axis dimension:\t");

		if (ext_arc_f_naxes[0] != ext_target_f_naxes[0] && ext_target_f_naxes[0] != ext_cont_f_naxes[0]) {	
	
			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -13, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		} else { 

			printf("OK\n"); 

		} 
	
		printf("Second axis dimension:\t"); 

		if (ext_arc_f_naxes[1] != ext_target_f_naxes[1] && ext_target_f_naxes[1] != ext_cont_f_naxes[1]) {
		
			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -14, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		} else { 

			printf("OK\n"); 

		} 

		// ***********************************************************************
		// Set the range limits using target fits file (ARG 2) n.b. this should
		// be an arbitrary choice if all files have identical parameters

		int cut_x [2] = {1, ext_target_f_naxes[0]};
		int cut_y [2] = {1, ext_target_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from arc/target/continuum fits 
		// files (ARGS 1, 2 and 3)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from arc/target/continuum fits
		// files (ARGS 1, 2 and 3)

		double ext_arc_f_pixels [nxelements];
		double ext_target_f_pixels [nxelements];
		double ext_cont_f_pixels [nxelements];

		// ***********************************************************************
		// Get arc fits file (ARG 1) values and store in 2D array

		int ii;

		double ext_arc_frame_values [nyelements][nxelements];
		memset(ext_arc_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(ext_arc_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(ext_arc_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, ext_arc_f_pixels, NULL, &ext_arc_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					ext_arc_frame_values[fpixel[1]-1][ii] = ext_arc_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -15, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_arc_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Get target fits file (ARG 2) values and store in 2D array

		double ext_target_frame_values [nyelements][nxelements];
		memset(ext_target_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(ext_target_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(ext_target_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, ext_target_f_pixels, NULL, &ext_target_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					ext_target_frame_values[fpixel[1]-1][ii] = ext_target_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -16, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_target_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Get continuum fits file (ARG 3) values and store in 2D array

		double ext_cont_frame_values [nyelements][nxelements];
		memset(ext_cont_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(ext_cont_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(ext_cont_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, ext_cont_f_pixels, NULL, &ext_cont_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					ext_cont_frame_values[fpixel[1]-1][ii] = ext_cont_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -17, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ext_cont_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

				return 1; 

			}

		}

		// CROSS CORRELATION AND DATA OFFSET/TRIM APPLICATION
		// ***********************************************************************
		// 1.	Cross correlate each row to find the best offsets using 
		//	[FRARCFIT_VAR_REF_FIBRE_INDEX_CC] as the reference row

		double reference_row_values [nxelements];
		memset(reference_row_values, 0, sizeof(double)*nxelements);

		for (ii=0; ii<nxelements; ii++) {

			reference_row_values[ii] = ext_arc_frame_values[FRARCFIT_VAR_REF_FIBRE_INDEX_CC][ii];
	
		}

		printf("\nOffset calculations");
		printf("\n-------------------\n");
		printf("\nFibre\tOffset\tR\n\n");

		int jj;

		double offset_row_values [nxelements];
		int best_offsets [nyelements];

		int this_best_offset;
		double this_best_r;

		for (jj=0; jj<nyelements; jj++) {

			memset(offset_row_values, 0, sizeof(double)*nxelements);

			for (ii=0; ii<nxelements; ii++) {

				offset_row_values[ii] = ext_arc_frame_values[jj][ii];
	
			}

			calculate_cross_correlation_offset(reference_row_values, offset_row_values, nxelements, max_cc_delay, &this_best_offset, &this_best_r);

			printf("%d\t%d\t%f\n", jj+1, this_best_offset, this_best_r);

			best_offsets[jj] = this_best_offset;

		}

		// 2.	Apply best offsets to arc/target/continuum 2D arrays and store

		double offset_ext_arc_frame_values [nyelements][nxelements];
		memset(offset_ext_arc_frame_values, 0, sizeof(double)*nxelements*nyelements);

		double offset_ext_target_frame_values [nyelements][nxelements];
		memset(offset_ext_target_frame_values, 0, sizeof(double)*nxelements*nyelements);

		double offset_ext_cont_frame_values [nyelements][nxelements];
		memset(offset_ext_cont_frame_values, 0, sizeof(double)*nxelements*nyelements);	

		int this_offset;

		for (jj=0; jj<nyelements; jj++) {

			this_offset = best_offsets[jj];

			for (ii=0; ii<nxelements; ii++) {

				if ((ii+this_offset < 0) || ii+this_offset >= nxelements) {	// make sure we don't insert into an out-of-bounds array element

					continue;

				} else {

					offset_ext_arc_frame_values[jj][ii] = ext_arc_frame_values[jj][ii+this_offset];
					offset_ext_target_frame_values[jj][ii] = ext_target_frame_values[jj][ii+this_offset];
					offset_ext_cont_frame_values[jj][ii] = ext_cont_frame_values[jj][ii+this_offset];

				}

			}

		}

		// 3.	Trim arc/target/continuum 2D arrays between [lag_offset] and
		//	[lead_offset] and store 

		int lead_offset = abs(gsl_stats_int_min(best_offsets, 1, nyelements));
		int lag_offset = abs(gsl_stats_int_max(best_offsets, 1, nyelements));

		int trim_nxelements = nxelements - abs(lead_offset) - abs(lag_offset);

		double **trim_offset_ext_arc_frame_values;							// Dynamically sized multidimensional arrays require ptr-to-ptr allocation
		trim_offset_ext_arc_frame_values = malloc(nyelements * sizeof(double *));			// and the corresponding initialisation of both columns (as pointers)

		for(jj=0; jj<nyelements; jj++) {

			trim_offset_ext_arc_frame_values[jj] = malloc(trim_nxelements * sizeof(double));	// and rows

		}

		double trim_offset_ext_target_frame_values [nyelements][trim_nxelements];
		memset(trim_offset_ext_target_frame_values, 0, sizeof(double)*trim_nxelements*nyelements);

		double trim_offset_ext_cont_frame_values [nyelements][trim_nxelements];
		memset(trim_offset_ext_cont_frame_values, 0, sizeof(double)*trim_nxelements*nyelements);	

		for (jj=0; jj<nyelements; jj++) {

			for (ii=0; ii<trim_nxelements; ii++) {

				trim_offset_ext_arc_frame_values[jj][ii] = offset_ext_arc_frame_values[jj][ii+lead_offset];
				trim_offset_ext_target_frame_values[jj][ii] = offset_ext_target_frame_values[jj][ii+lead_offset];
				trim_offset_ext_cont_frame_values[jj][ii] = offset_ext_cont_frame_values[jj][ii+lead_offset];

			}

		}

		// ARC LINE MATCHING
		// ***********************************************************************
		// 1.	Find contiguous peaks in arc frame

		int num_peaks = 0;	

		int **peaks;							// Dynamically sized multidimensional array to hold peak locations for all iterations
		peaks = malloc(nyelements * sizeof(int *));

		for(jj=0; jj<nyelements; jj++) {

			peaks[jj] = malloc(trim_nxelements * sizeof(int));

		}

		find_peaks_contiguous(trim_nxelements, nyelements, trim_offset_ext_arc_frame_values, peaks, &num_peaks, min_dist, half_aperture_num_pix, derivative_tol, derivative_tol_ref_px, pix_tolerance, INDEXING_CORRECTION);

		printf("\nArc line matching");
		printf("\n-----------------\n\n");
		printf("Candidate lines found:\t\t\t\t\t%d\n", num_peaks);

		if (num_peaks < min_contiguous_lines) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -18, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		}

		// 2.	Find parabolic centroids of contiguous peaks in arc 
		//	frame

		double row_values [trim_nxelements];

		double peak_centroids [nyelements][num_peaks];
		memset(peak_centroids, 0, sizeof(double)*nyelements*num_peaks);	

		double this_row_peak_centroids [num_peaks];

		int this_row_peaks [num_peaks];

		for (jj=0; jj<nyelements; jj++) {

			memset(row_values, 0, sizeof(double)*trim_nxelements);	

			for (ii=0; ii<trim_nxelements; ii++) {

				row_values[ii] = trim_offset_ext_arc_frame_values[jj][ii];

			}

			memset(this_row_peak_centroids, 0, sizeof(double)*num_peaks);	

			for (ii=0; ii<num_peaks; ii++) {

				this_row_peaks[ii] = peaks[jj][ii];

			}

			find_centroid_parabolic(row_values, this_row_peaks, num_peaks, this_row_peak_centroids, INDEXING_CORRECTION);

			for (ii=0; ii<num_peaks; ii++) {

				peak_centroids[jj][ii] = this_row_peak_centroids[ii];

			}

		}

		// 3.	Open reference arc line list file (ARG 11) and count
		//	number of reference lines
	
		FILE *arc_line_list_f;
	
		if (!check_file_exists(arc_line_list_filename)) { 

			arc_line_list_f = fopen(arc_line_list_filename , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -19, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1;

		}	
	
		char input_string[300];

		int arc_line_list_num_ref_lines = 0;
	
		while(!feof(arc_line_list_f)) {

			memset(input_string, '\0', sizeof(char)*300);

			fgets(input_string, 300, arc_line_list_f);	

			if (strtol(&input_string[0], NULL, 0) > 0) {		// check the line begins with a positive number (usable)
	
				arc_line_list_num_ref_lines++;				

			}
	
		}

		printf("Number of reference lines:\t\t\t\t%d\n", arc_line_list_num_ref_lines);

		// 4.	Read arc line list file (ARG 11) for reference
		//	wavelengths and their corresponding pixel positions

		rewind(arc_line_list_f);

		double arc_peak_centroids [arc_line_list_num_ref_lines];
		memset(arc_peak_centroids, 0, sizeof(double)*arc_line_list_num_ref_lines);

		double arc_peak_wavelengths [arc_line_list_num_ref_lines];
		memset(arc_peak_wavelengths, 0, sizeof(double)*arc_line_list_num_ref_lines);

		double pixel_channel, wavelength;

		int arc_line_list_line_index = 0;

		while(!feof(arc_line_list_f)) {

			memset(input_string, '\0', sizeof(char)*300);

			fgets(input_string, 300, arc_line_list_f);	

			if (strtol(&input_string[0], NULL, 0) > 0) {		// check the line begins with a positive number (usable)
	
				sscanf(input_string, "%lf\t%lf\t", &pixel_channel, &wavelength); 

				arc_peak_centroids[arc_line_list_line_index] = pixel_channel;
				arc_peak_wavelengths[arc_line_list_line_index] = wavelength;	

				arc_line_list_line_index++;	

			}
	
		}

		// for (ii=0; ii<arc_line_list_line_index; ii++) printf("\n%f\t%f", arc_peak_centroids[ii], arc_peak_wavelengths[ii]);	// DEBUG

		// 5.	Match identified peaks to lines in reference file

		double this_diff, least_diff, least_diff_index;
		bool matched_line;
		int matched_line_count = 0;

		int matched_line_indexes [num_peaks];

		for (ii=0; ii<num_peaks; ii++) {

			matched_line_indexes[ii] = -1;			// -1 indicates no match

		}

		double matched_line_diffs [num_peaks];
		memset(matched_line_diffs, 0, sizeof(double)*num_peaks);

		double average_pixel_channels [num_peaks];
		memset(average_pixel_channels, 0, sizeof(double)*num_peaks);	

		for (ii=0; ii<num_peaks; ii++) {

			for (jj=0; jj<nyelements; jj++) {

				average_pixel_channels[ii] += peak_centroids[jj][ii];

			}

			average_pixel_channels[ii] /= (double) nyelements;				// calculate the average pixel channel using all rows	

		}

		int duplicate_index;
		
		for (ii=0; ii<num_peaks; ii++) {							// for each identified candidate peak

			least_diff = 0;
			matched_line = FALSE;

			for (jj=0; jj<arc_line_list_num_ref_lines; jj++) {				// then take each line in the reference arc line list file

				this_diff = fabs(average_pixel_channels[ii] - arc_peak_centroids[jj]);	// calculate the difference between the identified and reference arc line

				if (this_diff <= max_pix_diff) {					// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
			
					if (matched_line == FALSE || this_diff < least_diff) {		// is this the first line found or a line with a smaller difference between ref arc line and matched?

						matched_line = TRUE;

						least_diff = this_diff;
						least_diff_index = jj;

					}
				
				}

			}

			if (matched_line == TRUE) {

				if (lsearch_int(matched_line_indexes, least_diff_index, num_peaks) != -1) {			// has this line already been allocated to another peak?

					duplicate_index = lsearch_int(matched_line_indexes, least_diff_index, num_peaks);

					if(least_diff < matched_line_diffs[duplicate_index]) {					// does it have a smaller difference? - comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

						matched_line_indexes[duplicate_index] = -1;					// reset the duplicate values, no increment required if duplicate
						matched_line_diffs[duplicate_index] = 0;

						matched_line_indexes[ii] = least_diff_index;
						matched_line_diffs[ii] = least_diff;
						
					}

				} else {

					matched_line_indexes[ii] = least_diff_index;
					matched_line_diffs[ii] = least_diff;
					matched_line_count++;

				}

			} else {

				matched_line_indexes[ii] = -1;						// not matched

			}

		}

		// 6.	Did we match enough lines?

		if (matched_line_count < min_matched_lines) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -20, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1;

		}

		printf("Number of matched lines:\t\t\t\t%d\n", matched_line_count);

		printf("\nIndex\tList wavelength\tList centroid\tAv. channel\tChannel difference\n");		// Channel difference is the difference between the reference arc list pixel location and the identified pixel location
		printf("\t(Å)\t\t(px)\t\t(px)\t\t(px)\n\n");	

		for (ii=0; ii<num_peaks; ii++) {

			if (matched_line_indexes[ii] == -1) { 	// this line was unmatched

				// printf("%d\t%s\t\t%s\t\t%.2f\t\t%s\n", ii, "", "", average_pixel_channels[ii], ""); 	// DEBUG
				continue;

			} else {

				printf("%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n", ii, arc_peak_wavelengths[matched_line_indexes[ii]], arc_peak_centroids[matched_line_indexes[ii]], average_pixel_channels[ii], matched_line_diffs[ii]);

			}

		}


		// 7.	And do they sample the distribution well?

		// for (ii=0; ii<arc_line_list_line_index; ii++) printf("\n%f\t%f\t%d", arc_peak_wavelengths[ii], arc_peak_wavelengths[matched_line_indexes[ii]], matched_line_indexes[ii]);	// DEBUG

		double list_total_dist = 0.0;
		double list_av_dist;

		for (ii=0; ii<arc_line_list_num_ref_lines-1; ii++) {

			for (jj=ii+1; jj<arc_line_list_num_ref_lines; jj++) {
				
				list_total_dist += arc_peak_wavelengths[jj] - arc_peak_wavelengths[ii];

				// printf("\n%f", arc_peak_wavelengths[jj] - arc_peak_wavelengths[ii]);		// DEBUG

			}

		}

		list_av_dist = list_total_dist / powf(arc_line_list_line_index, 2);

		double matched_line_wavelengths [matched_line_count];
		memset(matched_line_wavelengths, 0, sizeof(double)*matched_line_count);
		
		int this_matched_line_index = 0;

		for (ii=0; ii<num_peaks; ii++) {
		
			if (matched_line_indexes[ii] == -1) { 	// this line was unmatched

				continue;

			} else {

				matched_line_wavelengths[this_matched_line_index] = arc_peak_wavelengths[matched_line_indexes[ii]];
				this_matched_line_index++;

			}

		}

		double sample_total_dist = 0.0;
		double sample_av_dist;

		for (ii=0; ii<matched_line_count-1; ii++) {

			for (jj=ii+1; jj<matched_line_count; jj++) {
				
				sample_total_dist += matched_line_wavelengths[jj] - matched_line_wavelengths[ii];

				// printf("\n%f", matched_line_wavelengths[jj] - matched_line_wavelengths[ii]);		// DEBUG

			}

		}

		sample_av_dist = sample_total_dist / powf(matched_line_count, 2);

		double sample_list_diff = abs(list_av_dist-sample_av_dist);

		printf("\nAverage distance between lines in reference arc list:\t%.1f", list_av_dist);
		printf("\nAverage distance between lines in matched line list:\t%.1f\n", sample_av_dist);

		if (sample_list_diff > max_av_wavelength_diff) {	// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

			RETURN_FLAG = 2;

		}

		// FIND THE DISPERSION SOLUTION AND WRITE TO 
		// [FRARCFIT_OUTPUTF_WAVFITS_FILE] OUTPUT FILE
		// ***********************************************************************

		// 1.	Perform a few checks to ensure the input fitting parameters 
		// 	are sensible

		if ((fit_order < FRARCFIT_VAR_POLYORDER_LO) || (fit_order > FRARCFIT_VAR_POLYORDER_HI)) {	// Check [fit_order] is within config limits

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -21, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1; 

		}

		// 2.	Create [FRARCFIT_OUTPUTF_WAVFITS_FILE] output file and print a few 
		// 	parameters

		FILE *outputfile;
		outputfile = fopen(FRARCFIT_OUTPUTF_WAVFITS_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -22, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1;

		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", FRARCFIT_OUTPUTF_WAVFITS_FILE);
	        fprintf(outputfile, "# List of fibre pixel channel to wavelength dispersion solution coefficients and corresponding chi-squareds found using the frarcfit program.\n\n");
                fprintf(outputfile, "# Run Datetime:\t\t%s\n\n", timestr);
	        fprintf(outputfile, "# Target Filename:\t%s\n", ext_target_f);
	        fprintf(outputfile, "# Arc Filename:\t\t%s\n", ext_arc_f);
	        fprintf(outputfile, "# Continuum Filename:\t%s\n\n", ext_cont_f);
                fprintf(outputfile, "# Polynomial Order:\t%d\n\n", fit_order);

		// 3.	Populate arrays and perform fits

		double coord_x [matched_line_count];
		double coord_y [matched_line_count];

		double coeffs [fit_order+1];

		double this_chi_squared, chi_squared_min = 0.0, chi_squared_max = 0.0, chi_squared = 0.0;

		int line_count;

		for (jj=0; jj<nyelements; jj++) {

			memset(coord_x, 0, sizeof(double)*matched_line_count);
			memset(coord_y, 0, sizeof(double)*matched_line_count);
			memset(coeffs, 0, sizeof(double)*fit_order+1);

			line_count = 0;

			for (ii=0; ii<num_peaks; ii++) {
		
				if (matched_line_indexes[ii] == -1) { 				// this line was unmatched

					continue;

				} else {

					coord_x[line_count] = peak_centroids[jj][ii];
					coord_y[line_count] = arc_peak_wavelengths[matched_line_indexes[ii]];
					line_count++;

				}

			}

			if (calc_least_sq_fit(fit_order, matched_line_count, coord_x, coord_y, coeffs, &this_chi_squared)) {	

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -23, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				for(ii=0; ii<nyelements; ii++) {

					free(trim_offset_ext_arc_frame_values[ii]);
					free(peaks[ii]);

				}

				free(peaks);
				free(trim_offset_ext_arc_frame_values);

				fclose(arc_line_list_f);
				fclose(outputfile);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

				return 1; 

			}

			// Output solutions to [FRARCFIT_OUTPUTF_WAVFITS_FILE] file 

			fprintf(outputfile, "%d\t", jj+1);

			for (ii=0; ii<=fit_order; ii++) {
		  
				fprintf(outputfile, FRARCFIT_VAR_ACCURACY_COEFFS, coeffs[ii]);
				fprintf(outputfile, "\t");

			}

			fprintf(outputfile, FRARCFIT_VAR_ACCURACY_CHISQ, this_chi_squared);
			fprintf(outputfile, "\n");

			if ((jj==0) || (this_chi_squared < chi_squared_min)) { 

				chi_squared_min = this_chi_squared;

			} else if ((jj==0) || (this_chi_squared > chi_squared_max)) {

				chi_squared_max = this_chi_squared;

			}

			// printf("%d\t%f\n", ii, this_chi_squared);	// DEBUG

			chi_squared += this_chi_squared;

		}

		fprintf(outputfile, "%d", EOF);

		printf("\nFitting results");
		printf("\n---------------\n");
		printf("\nMin χ2:\t\t\t%.2f\n", chi_squared_min);
		printf("Max χ2:\t\t\t%.2f\n", chi_squared_max);
		printf("Average χ2:\t\t%.2f\n", chi_squared/nyelements);

		// 4.	Perform a few checks to ensure the chi squareds are sensible 

		if ((chi_squared_min < FRARCFIT_VAR_CHISQUARED_MIN) || (chi_squared_max > FRARCFIT_VAR_CHISQUARED_MAX)) {

			RETURN_FLAG = 3;

		}

		// ***********************************************************************
		// Create [cc_ext_arc/target/continuum_frame_values] array to hold the 
		// output data in the correct format

		int kk;

		double cc_ext_arc_frame_values [nyelements*trim_nxelements];
		memset(cc_ext_arc_frame_values, 0, sizeof(double)*nyelements*trim_nxelements);

		double cc_ext_target_frame_values [nyelements*trim_nxelements];
		memset(cc_ext_target_frame_values, 0, sizeof(double)*nyelements*trim_nxelements);

		double cc_ext_cont_frame_values [nyelements*trim_nxelements];
		memset(cc_ext_cont_frame_values, 0, sizeof(double)*nyelements*trim_nxelements);

		for (jj=0; jj<nyelements; jj++) {
	
			ii = jj * trim_nxelements;
	
			for (kk=0; kk<trim_nxelements; kk++) {
	
				cc_ext_arc_frame_values[ii] = trim_offset_ext_arc_frame_values[jj][kk];
				cc_ext_target_frame_values[ii] = trim_offset_ext_target_frame_values[jj][kk];
				cc_ext_cont_frame_values[ii] = trim_offset_ext_cont_frame_values[jj][kk];
				ii++;
	
			}
	
		}

		// ***********************************************************************
		// Set output frame (ARGS 16,17 and 18) parameters

		fitsfile *cc_ext_arc_f_ptr, *cc_ext_target_f_ptr, *cc_ext_cont_f_ptr;

		int cc_ext_arc_f_status = 0, cc_ext_target_f_status = 0, cc_ext_cont_f_status = 0;

		long cc_ext_arc_f_naxes [2] = {trim_nxelements, nyelements};
		long cc_ext_target_f_naxes [2] = {trim_nxelements, nyelements};
		long cc_ext_cont_f_naxes [2] = {trim_nxelements, nyelements};
	
		long cc_ext_arc_f_fpixel = 1, cc_ext_target_f_fpixel = 1, cc_ext_cont_f_fpixel = 1;

		// ***********************************************************************
		// Create and write [cc_ext_arc_frame_values] to output file (ARG 16)
	
		if (!fits_create_file(&cc_ext_arc_f_ptr, cc_ext_arc_f, &cc_ext_arc_f_status)) {
	
			if (!fits_create_img(cc_ext_arc_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, cc_ext_arc_f_naxes, &cc_ext_arc_f_status)) {

				if (!fits_write_img(cc_ext_arc_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], cc_ext_arc_f_fpixel, nyelements * trim_nxelements, cc_ext_arc_frame_values, &cc_ext_arc_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -24, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, cc_ext_arc_f_status); 

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					for(ii=0; ii<nyelements; ii++) {

						free(trim_offset_ext_arc_frame_values[ii]);
						free(peaks[ii]);

					}

					free(peaks);
					free(trim_offset_ext_arc_frame_values);

					fclose(arc_line_list_f);
					fclose(outputfile);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
					if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
					if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
					if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -25, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_arc_f_status); 


				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				for(ii=0; ii<nyelements; ii++) {

					free(trim_offset_ext_arc_frame_values[ii]);
					free(peaks[ii]);

				}

				free(peaks);
				free(trim_offset_ext_arc_frame_values);

				fclose(arc_line_list_f);
				fclose(outputfile);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
				if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -26, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cc_ext_arc_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);
			fclose(outputfile);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 

			return 1;

		}

		// ***********************************************************************
		// Create and write [cc_ext_target_frame_values] to output file (ARG 17)
	
		if (!fits_create_file(&cc_ext_target_f_ptr, cc_ext_target_f, &cc_ext_target_f_status)) {
	
			if (!fits_create_img(cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, cc_ext_target_f_naxes, &cc_ext_target_f_status)) {

				if (!fits_write_img(cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], cc_ext_target_f_fpixel, nyelements * trim_nxelements, cc_ext_target_frame_values, &cc_ext_target_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -27, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, cc_ext_target_f_status); 

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					for(ii=0; ii<nyelements; ii++) {

						free(trim_offset_ext_arc_frame_values[ii]);
						free(peaks[ii]);

					}

					free(peaks);
					free(trim_offset_ext_arc_frame_values);

					fclose(arc_line_list_f);
					fclose(outputfile);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
					if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
					if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
					if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
					if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -28, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_target_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				for(ii=0; ii<nyelements; ii++) {

					free(trim_offset_ext_arc_frame_values[ii]);
					free(peaks[ii]);

				}

				free(peaks);
				free(trim_offset_ext_arc_frame_values);

				fclose(arc_line_list_f);
				fclose(outputfile);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
				if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
				if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -29, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cc_ext_target_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);
			fclose(outputfile);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Create and write [cc_ext_cont_frame_values] to output file (ARG 18)
	
		if (!fits_create_file(&cc_ext_cont_f_ptr, cc_ext_cont_f, &cc_ext_cont_f_status)) {
	
			if (!fits_create_img(cc_ext_cont_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, cc_ext_cont_f_naxes, &cc_ext_cont_f_status)) {

				if (!fits_write_img(cc_ext_cont_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], cc_ext_cont_f_fpixel, nyelements * trim_nxelements, cc_ext_cont_frame_values, &cc_ext_cont_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -30, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, cc_ext_cont_f_status); 

					free(ext_arc_f);
					free(ext_target_f);
					free(ext_cont_f);
					free(arc_line_list_filename);
					free(cc_ext_arc_f);
					free(cc_ext_target_f);
					free(cc_ext_cont_f);

					for(ii=0; ii<nyelements; ii++) {

						free(trim_offset_ext_arc_frame_values[ii]);
						free(peaks[ii]);

					}

					free(peaks);
					free(trim_offset_ext_arc_frame_values);

					fclose(arc_line_list_f);
					fclose(outputfile);

					if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
					if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
					if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
					if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
					if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
					if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -31, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_cont_f_status); 

				free(ext_arc_f);
				free(ext_target_f);
				free(ext_cont_f);
				free(arc_line_list_filename);
				free(cc_ext_arc_f);
				free(cc_ext_target_f);
				free(cc_ext_cont_f);

				for(ii=0; ii<nyelements; ii++) {

					free(trim_offset_ext_arc_frame_values[ii]);
					free(peaks[ii]);

				}

				free(peaks);
				free(trim_offset_ext_arc_frame_values);

				fclose(arc_line_list_f);
				fclose(outputfile);

				if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
				if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
				if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
				if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
				if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
				if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -32, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cc_ext_cont_f_status); 

			free(ext_arc_f);
			free(ext_target_f);
			free(ext_cont_f);
			free(arc_line_list_filename);
			free(cc_ext_arc_f);
			free(cc_ext_target_f);
			free(cc_ext_cont_f);

			for(ii=0; ii<nyelements; ii++) {

				free(trim_offset_ext_arc_frame_values[ii]);
				free(peaks[ii]);

			}

			free(peaks);
			free(trim_offset_ext_arc_frame_values);

			fclose(arc_line_list_f);
			fclose(outputfile);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Clean up heap memory

		free(ext_arc_f);
		free(ext_target_f);
		free(ext_cont_f);
		free(arc_line_list_filename);
		free(cc_ext_arc_f);
		free(cc_ext_target_f);
		free(cc_ext_cont_f);

		for(ii=0; ii<nyelements; ii++) {

			free(trim_offset_ext_arc_frame_values[ii]);
			free(peaks[ii]);

		}

		free(peaks);
		free(trim_offset_ext_arc_frame_values);

		// ***********************************************************************
		// Close input files (ARGS 1,2 and 3), output files (ARGS 16, 17 and 18),
		// arc list file (ARG 11) and [FRARCFIT_OUTPUTF_WAVFITS_FILE] file

		if (fclose(arc_line_list_f)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -33, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(outputfile);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

		}

		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -34, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

			if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

		}

		if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -35, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, ext_arc_f_status); 

			if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) fits_report_error (stdout, ext_target_f_status); 
			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

	    	}

		if(fits_close_file(ext_target_f_ptr, &ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -36, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, ext_target_f_status); 

			if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) fits_report_error (stdout, ext_cont_f_status); 
			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

	    	}

		if(fits_close_file(ext_cont_f_ptr, &ext_cont_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -37, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, ext_cont_f_status); 

			if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) fits_report_error (stdout, cc_ext_arc_f_status); 
			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

	    	}

		if(fits_close_file(cc_ext_arc_f_ptr, &cc_ext_arc_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -38, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cc_ext_arc_f_status); 

			if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) fits_report_error (stdout, cc_ext_target_f_status); 
			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

	    	}

		if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -39, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cc_ext_target_f_status); 

			if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) fits_report_error (stdout, cc_ext_cont_f_status);

			return 1; 

	    	}

		if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -40, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cc_ext_cont_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", RETURN_FLAG, "Status flag for L2 frarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

