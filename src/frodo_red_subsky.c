/************************************************************************

 File:				frodo_red_subsky.c
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
#include "frodo_red_subsky.h"

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

// *********************************************************************/

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 5) {

		if(populate_env_variable(FRS_BLURB_FILE, "L2_FRS_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FRS_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -1, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *reb_cor_cc_ext_target_f		= strdup(argv[1]);
		int clip_sigma				= strtol(argv[2], NULL, 0);
		int thresh_median_sigma			= strtol(argv[3], NULL, 0);
		char *ss_reb_cor_cc_ext_target_f	= strdup(argv[4]);

		// ***********************************************************************
		// Open cc extracted target file (ARG 1), get parameters and perform any  
		// data format checks 

		fitsfile *reb_cor_cc_ext_target_f_ptr;

		int reb_cor_cc_ext_target_f_maxdim = 2;
		int reb_cor_cc_ext_target_f_status = 0, reb_cor_cc_ext_target_f_bitpix, reb_cor_cc_ext_target_f_naxis;
		long reb_cor_cc_ext_target_f_naxes [2] = {1,1};

		if(!fits_open_file(&reb_cor_cc_ext_target_f_ptr, reb_cor_cc_ext_target_f, READONLY, &reb_cor_cc_ext_target_f_status)) {

			if(!populate_img_parameters(reb_cor_cc_ext_target_f, reb_cor_cc_ext_target_f_ptr, reb_cor_cc_ext_target_f_maxdim, &reb_cor_cc_ext_target_f_bitpix, &reb_cor_cc_ext_target_f_naxis, reb_cor_cc_ext_target_f_naxes, &reb_cor_cc_ext_target_f_status, "TARGET FRAME")) {

				if (reb_cor_cc_ext_target_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -2, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(reb_cor_cc_ext_target_f);
					free(ss_reb_cor_cc_ext_target_f);

					if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -3, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, reb_cor_cc_ext_target_f_status); 

				free(reb_cor_cc_ext_target_f);
				free(ss_reb_cor_cc_ext_target_f);

				if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -4, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, reb_cor_cc_ext_target_f_status); 

			free(reb_cor_cc_ext_target_f);
			free(ss_reb_cor_cc_ext_target_f);

			return 1; 

		}

		// ***********************************************************************
		// Set the range limits using target fits file (ARG 1)

		int cut_x [2] = {1, reb_cor_cc_ext_target_f_naxes[0]};
		int cut_y [2] = {1, reb_cor_cc_ext_target_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from target fits file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create array to store pixel values from target fits file (ARG 1)

		double reb_cor_cc_ext_target_f_pixels [nxelements];

		// CHECK IF SKY FIBRES EXIST
		// ***********************************************************************
		// 1.	Read row of target frame

		int this_fibre_index;
		double this_fibre_flux;

		double fibre_fluxes [nyelements];
		memset(fibre_fluxes, 0, sizeof(double)*nyelements);

		int ii;

		double reb_cor_cc_ext_target_frame_values [nyelements][nxelements];
		memset(reb_cor_cc_ext_target_frame_values, 0, sizeof(double)*nyelements*nxelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			this_fibre_index = fpixel[1] - 1;
			this_fibre_flux = 0.0;

			memset(reb_cor_cc_ext_target_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(reb_cor_cc_ext_target_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, reb_cor_cc_ext_target_f_pixels, NULL, &reb_cor_cc_ext_target_f_status)) {

				// 2.	Sum the total flux in each fibre and store to 2D array (for use later in sky subtraction)

				for (ii=0; ii<nxelements; ii++) {

					this_fibre_flux += reb_cor_cc_ext_target_f_pixels[ii];

					reb_cor_cc_ext_target_frame_values[this_fibre_index][ii] = reb_cor_cc_ext_target_f_pixels[ii];

				}

				fibre_fluxes[this_fibre_index] = this_fibre_flux;

				// printf("%f\n", this_fibre_flux);	// DEBUG
				
			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -5, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, reb_cor_cc_ext_target_f_status); 

				free(reb_cor_cc_ext_target_f);
				free(ss_reb_cor_cc_ext_target_f);

				if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

				return 1; 

			}

		}

		// 3.	Perform iterative sigma clip

		int sky_indexes [nyelements];	// array of boolean values, TRUE = sky only

		printf("\nSigma clipping");
		printf("\n--------------\n");

		double final_mean, final_sd;
		int num_sky_fibres;

		iterative_sigma_clip(fibre_fluxes, nyelements, clip_sigma, sky_indexes, &final_mean, &final_sd, &num_sky_fibres);

		/* // !!SUBSKY TESTS!! - comment iterative_sigma_clip + make half of [sky_indexes] array TRUE and set [num_sky_fibres]
		for (ii=0; ii<72; ii++) sky_indexes[ii] = TRUE; 
		for (ii=72; ii<143; ii++) sky_indexes[ii] = FALSE; 
		num_sky_fibres = 72; */
		
		// /* !!SUBSKY TESTS!! - remove checks

		// REMOVE AVERAGE SKY CONTRIBUTION FROM ALL FIBRES
		// ***********************************************************************
		// 1.	Simple checks to see if either all fibres have been identified
		// 	as sky or if we have no sky fibres given the selection criteria
		
		if (num_sky_fibres == nyelements) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -6, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(reb_cor_cc_ext_target_f);
			free(ss_reb_cor_cc_ext_target_f);

			if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

			return 1; 	

		} else if (num_sky_fibres == 0) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -7, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(reb_cor_cc_ext_target_f);
			free(ss_reb_cor_cc_ext_target_f);

			if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

			return 1; 	

		} 

		// */ // !!SUBSKY TESTS!! - remove checks

		// 2.	Duplicate [fibre_fluxes] array, sort into ascending order and
		//	calculate median

		double fibre_fluxes_sorted [nyelements];
		memcpy(fibre_fluxes_sorted, fibre_fluxes, sizeof(double)*nyelements);	

		gsl_sort(fibre_fluxes_sorted, 1, nyelements);
		double fibre_fluxes_median = gsl_stats_median_from_sorted_data(fibre_fluxes_sorted, 1, nyelements);

		int jj;

		double ss_reb_cor_cc_ext_target_frame_values [nyelements][nxelements];
		memset(ss_reb_cor_cc_ext_target_frame_values, 0, sizeof(double)*nyelements*nxelements);

		// 3.	Check to see if median of initial sample is within 
		//	[thresh_median_sigma] sigma of final mean. If this is not the
		//	case, it indicates that the IFU is over 50% filled with object
		//	fibres and sky subtraction mightn't be reliable or the source
		//	is extended.

		// printf("%f\t%f\t%f\n", fibre_fluxes_median, final_mean, thresh_median_sigma*final_sd);	// DEBUG
 		
		// /* !!SUBSKY TESTS!! - remove checks	
		if ((fibre_fluxes_median > final_mean - thresh_median_sigma*final_sd) && (fibre_fluxes_median < final_mean + thresh_median_sigma*final_sd)) { // Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
		// */ // !!SUBSKY TESTS!! - remove checks	
			int sky_fibre_index;

			double this_wavelength_sky_fibre_values_median;

			// 4.	Take median of sky fibre fluxes

			for (ii=0; ii<nxelements; ii++) {

				double *this_wavelength_sky_fibre_values;
				this_wavelength_sky_fibre_values = (double *) malloc((num_sky_fibres)*sizeof(double));

				sky_fibre_index = 0;
			
				for (jj=0; jj<nyelements; jj++) {

					if (sky_indexes[jj] == TRUE) {

						this_wavelength_sky_fibre_values[sky_fibre_index] = reb_cor_cc_ext_target_frame_values[jj][ii];
						sky_fibre_index++;

					} 
	
				}

				gsl_sort(this_wavelength_sky_fibre_values, 1, num_sky_fibres);
				this_wavelength_sky_fibre_values_median = gsl_stats_median_from_sorted_data(this_wavelength_sky_fibre_values, 1, num_sky_fibres);

				// 5.	Subtract the median from the target frame values

				for (jj=0; jj<nyelements; jj++) {

					ss_reb_cor_cc_ext_target_frame_values[jj][ii] = reb_cor_cc_ext_target_frame_values[jj][ii] - this_wavelength_sky_fibre_values_median;

				}	

				free(this_wavelength_sky_fibre_values);

			}
		
		// /* !!SUBSKY TESTS!! - remove checks
		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -8, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(reb_cor_cc_ext_target_f);
			free(ss_reb_cor_cc_ext_target_f);

			if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

			return 1; 

		}
		// */ // !!SUBSKY TESTS!! - remove checks

		// ***********************************************************************
		// Set ss_reb_cor_cc_ext_target frame parameters

		fitsfile *ss_reb_cor_cc_ext_target_f_ptr;
	
		int ss_reb_cor_cc_ext_target_f_status = 0;
		long ss_reb_cor_cc_ext_target_f_naxes [2] = {nxelements,nyelements};
	
		long ss_reb_cor_cc_ext_target_f_fpixel = 1;

		// ***********************************************************************
		// Create [ss_reb_cor_cc_ext_target_output_frame_values] array to hold the
		// output data in the correct format

		double ss_reb_cor_cc_ext_target_output_frame_values [nxelements*nyelements];
		memset(ss_reb_cor_cc_ext_target_output_frame_values, 0, sizeof(double)*nxelements*nyelements);

		int kk;

		for (ii=0; ii<nyelements; ii++) {
	
			jj = ii * nxelements;
	
			for (kk=0; kk<nxelements; kk++) {
	
				ss_reb_cor_cc_ext_target_output_frame_values[jj] = ss_reb_cor_cc_ext_target_frame_values[ii][kk];
				jj++;

			}
		
		}

		// ***********************************************************************
		// Create and write [cor_cc_ext_target_output_frame_values] to output file
		// (ARG 4)	
	
		if (!fits_create_file(&ss_reb_cor_cc_ext_target_f_ptr, ss_reb_cor_cc_ext_target_f, &ss_reb_cor_cc_ext_target_f_status)) {
	
			if (!fits_create_img(ss_reb_cor_cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, ss_reb_cor_cc_ext_target_f_naxes, &ss_reb_cor_cc_ext_target_f_status)) {

				if (!fits_write_img(ss_reb_cor_cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], ss_reb_cor_cc_ext_target_f_fpixel, nxelements * nyelements, ss_reb_cor_cc_ext_target_output_frame_values, &ss_reb_cor_cc_ext_target_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -9, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, ss_reb_cor_cc_ext_target_f_status); 

					free(reb_cor_cc_ext_target_f);
					free(ss_reb_cor_cc_ext_target_f);

					if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 
					if(fits_close_file(ss_reb_cor_cc_ext_target_f_ptr, &ss_reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, ss_reb_cor_cc_ext_target_f_status); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -10, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, ss_reb_cor_cc_ext_target_f_status); 

				free(reb_cor_cc_ext_target_f);
				free(ss_reb_cor_cc_ext_target_f);

				if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 
				if(fits_close_file(ss_reb_cor_cc_ext_target_f_ptr, &ss_reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, ss_reb_cor_cc_ext_target_f_status); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -11, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, ss_reb_cor_cc_ext_target_f_status); 

			free(reb_cor_cc_ext_target_f);
			free(ss_reb_cor_cc_ext_target_f);

			if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Clean up heap memory

		free(reb_cor_cc_ext_target_f);
		free(ss_reb_cor_cc_ext_target_f);

		// ***********************************************************************
		// Close input file (ARG 1) and output file (ARG 4)

		if(fits_close_file(reb_cor_cc_ext_target_f_ptr, &reb_cor_cc_ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -12, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, reb_cor_cc_ext_target_f_status); 

			if(fits_close_file(ss_reb_cor_cc_ext_target_f_ptr, &ss_reb_cor_cc_ext_target_f_status)) fits_report_error (stdout, ss_reb_cor_cc_ext_target_f_status); 

			return 1; 

	    	}

		if(fits_close_file(ss_reb_cor_cc_ext_target_f_ptr, &ss_reb_cor_cc_ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", -13, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, ss_reb_cor_cc_ext_target_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATSU", RETURN_FLAG, "Status flag for L2 frsubsky routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

