/************************************************************************

 File:				frodo_red_correct_throughput.c
 Last Modified Date:     	03/03/11

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
#include "frodo_red_correct_throughput.h"
#include "frodo_red_arcfit.h"

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

// *********************************************************************/

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 6) {

		if(populate_env_variable(HEADER_FILE, "L2_HEADER_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", 1, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		if(populate_env_variable(FRCT_BLURB_FILE, "L2_FRCT_BLURB_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", 2, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		print_file(HEADER_FILE);
		print_file(FRCT_BLURB_FILE);

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -1, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters
	
		char *cc_ext_target_f		= strdup(argv[1]);
		char *cc_ext_cont_f		= strdup(argv[2]);
		double start_wav		= strtod(argv[3], NULL);
		double end_wav			= strtod(argv[4], NULL);
		char *cor_cc_ext_target_f	= strdup(argv[5]);

		// ***********************************************************************
		// Open cc extracted target file (ARG 1), get parameters and perform any  
		// data format checks 

		fitsfile *cc_ext_target_f_ptr;

		int cc_ext_target_f_maxdim = 2;
		int cc_ext_target_f_status = 0, cc_ext_target_f_bitpix, cc_ext_target_f_naxis;
		long cc_ext_target_f_naxes [2] = {1,1};

		if(!fits_open_file(&cc_ext_target_f_ptr, cc_ext_target_f, READONLY, &cc_ext_target_f_status)) {

			if(!populate_img_parameters(cc_ext_target_f, cc_ext_target_f_ptr, cc_ext_target_f_maxdim, &cc_ext_target_f_bitpix, &cc_ext_target_f_naxis, cc_ext_target_f_naxes, &cc_ext_target_f_status, "TARGET FRAME")) {

				if (cc_ext_target_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -2, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -3, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_target_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -4, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cc_ext_target_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Open cc extracted continuum file (ARG 2), get parameters and perform  
		// any data format checks 

		fitsfile *cc_ext_cont_f_ptr;

		int cc_ext_cont_f_maxdim = 2;
		int cc_ext_cont_f_status = 0, cc_ext_cont_f_bitpix, cc_ext_cont_f_naxis;
		long cc_ext_cont_f_naxes [2] = {1,1};

		if(!fits_open_file(&cc_ext_cont_f_ptr, cc_ext_cont_f, READONLY, &cc_ext_cont_f_status)) {

			if(!populate_img_parameters(cc_ext_cont_f, cc_ext_cont_f_ptr, cc_ext_cont_f_maxdim, &cc_ext_cont_f_bitpix, &cc_ext_cont_f_naxis, cc_ext_cont_f_naxes, &cc_ext_cont_f_status, "CONTINUUM FRAME")) {

				if (cc_ext_cont_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -5, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -6, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_cont_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -7, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cc_ext_cont_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Check consistency of extracted target/continuum fits files 
		// (ARGS 1 and 2)
	
		printf("\nConsistency check");
		printf("\n-----------------\n");

		printf("\nBits per pixel:\t\t");
		if (cc_ext_target_f_bitpix != cc_ext_cont_f_bitpix) {
			printf("FAIL\n"); 
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -8, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			return 1; 
		} else { 
			printf("OK\n"); 
		} 

		printf("Number of axes:\t\t");
		if (cc_ext_target_f_naxis != cc_ext_cont_f_naxis) {		
			printf("FAIL\n"); 
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -9, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			return 1; 
		} else { 
			printf("OK\n"); 
		} 
	
		printf("First axis dimension:\t");
		if (cc_ext_target_f_naxes[0] != cc_ext_cont_f_naxes[0]) {		
			printf("FAIL\n"); 
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -10, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			return 1; 
		} else { 
			printf("OK\n"); 
		} 
	
		printf("Second axis dimension:\t"); 
		if (cc_ext_target_f_naxes[1] != cc_ext_cont_f_naxes[1]) {		
			printf("FAIL\n"); 
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -11, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			return 1; 
		} else { 
			printf("OK\n"); 
		} 

		// ***********************************************************************
		// Set the range limits using target fits file (ARG 1) n.b. this should
		// be an arbitrary choice if all files have identical parameters

		int cut_x [2] = {1, cc_ext_target_f_naxes[0]};
		int cut_y [2] = {1, cc_ext_target_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from target/continuum fits 
		// files (ARGS 1 and 2)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from target/continuum fits
		// files (ARGS 1, 2 and 3)

		double cc_ext_target_f_pixels [nxelements];
		double cc_ext_cont_f_pixels [nxelements];

		// ***********************************************************************
		// Get cc extracted continuum fits file (ARG 3) values and store in 2D 
		// array

		int ii;

		double cc_ext_cont_frame_values [nyelements][nxelements];
		memset(cc_ext_cont_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(cc_ext_cont_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(cc_ext_cont_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, cc_ext_cont_f_pixels, NULL, &cc_ext_cont_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					cc_ext_cont_frame_values[fpixel[1]-1][ii] = cc_ext_cont_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -12, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_cont_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Open [FRARCFIT_OUTPUTF_WAVFITS_FILE] dispersion solutions file

		FILE *dispersion_solutions_f;
	
		if (!check_file_exists(FRARCFIT_OUTPUTF_WAVFITS_FILE)) { 

			dispersion_solutions_f = fopen(FRARCFIT_OUTPUTF_WAVFITS_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -13, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}	

		// ***********************************************************************
		// Find some [FRARCFIT_OUTPUTF_WAVFITS_FILE] file details

		char input_string [500];

		bool find_polynomialorder_comment = FALSE;

		int polynomial_order;	

		char search_string_1 [20] = "# Polynomial Order:\0";	// this is the comment to be found from the [FRARCFIT_OUTPUTF_WAVFITS_FILE] file

		while(!feof(dispersion_solutions_f)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, dispersion_solutions_f);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &polynomial_order);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_polynomialorder_comment = 1;
				break;


			} 

		}

		if (find_polynomialorder_comment != 1) {	// error check - didn't find the comment in the [FRARCFIT_OUTPUTF_WAVFITS_FILE] file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -14, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}

		// ***********************************************************************
		// Rewind and extract coefficients from [FRARCFIT_OUTPUTF_WAVFITS_FILE]
		// file 

		rewind(dispersion_solutions_f);

		int token_index;	// this variable will hold which token we're dealing with
		int coeff_index;	// this variable will hold which coefficient we're dealing with
		int this_fibre;	
		double this_coeff;
		double this_chisquared;
	
		char *token;

		double coeffs [FIBRES][polynomial_order+1];
		memset(coeffs, 0, sizeof(double)*FIBRES*(polynomial_order+1));

		while(!feof(dispersion_solutions_f)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, dispersion_solutions_f);

			token_index = 0;
			coeff_index = 0;

			if (strtol(&input_string[0], NULL, 0) > 0) { 		// check the line begins with a positive number

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

						this_fibre = strtol(token, NULL, 0);

					} else if ((token_index >= 1) && (token_index <= polynomial_order+1)) { 	// coeff token

						this_coeff = strtod(token, NULL);
						// printf("%d\t%d\t%e\n", this_fibre, coeff_index, this_coeff);		// DEBUG
						coeffs[this_fibre-1][coeff_index] = this_coeff;
						coeff_index++;

					} else if (token_index == polynomial_order+2) {					// chisquared token

						this_chisquared = strtod(token, NULL);

					}

					token_index++;

				}

			}

		}	

		// ***********************************************************************
		// Find wavelength extremities from [FRARCFIT_OUTPUTF_WAVFITS_FILE] file
		// and ensure the input constraints [start_wav] (ARG 3) and [end_wav]
		// (ARG 4) don't lie outside these boundaries

		double this_fibre_smallest_wav, this_fibre_largest_wav, smallest_wav, largest_wav;

		int jj;

		for (ii=0; ii<FIBRES; ii++) {

			this_fibre_smallest_wav = 0.0; this_fibre_largest_wav = 0.0;

			for (jj=0; jj<=polynomial_order; jj++) {
		    
				this_fibre_smallest_wav += coeffs[ii][jj]*pow(0+INDEXING_CORRECTION, jj);
				this_fibre_largest_wav += coeffs[ii][jj]*pow((cut_x[1]-1)+INDEXING_CORRECTION, jj);

			}

			if (ii==0) {

				smallest_wav = this_fibre_smallest_wav;
				largest_wav = this_fibre_largest_wav;

			} else if (this_fibre_smallest_wav > smallest_wav) {    // note the sign, we want to find the maximum value for the least wavelength. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				smallest_wav = this_fibre_smallest_wav;

			} else if (this_fibre_largest_wav < largest_wav) {      // logic is as above, except for least value of maximum wavelength. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				largest_wav = this_fibre_largest_wav;
			  
			}
		
		}

		// ***********************************************************************
	        // Need to find pixel indexes for starting/ending wavelength ranges for
		// each fibre

	        double this_element_wav;

	        int first_element_index_array [nyelements];
		memset(first_element_index_array, 0, sizeof(int)*nyelements);

	        int last_element_index_array [nyelements];
		memset(last_element_index_array, 0, sizeof(int)*nyelements);

	        int kk;

		for (jj=0; jj<nyelements; jj++) {

			for (ii=0; ii<nxelements; ii++) {

				this_element_wav = 0.0;
	
				for (kk=0; kk<=polynomial_order; kk++) {
		      
		        		this_element_wav += coeffs[jj][kk]*pow(ii,kk);

		      		}

		      		if (this_element_wav > start_wav) {	// the current index, ii, represents the first pixel with a wavelength >= start_wav. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

		       			break;

		     		}

		   	}

			first_element_index_array[jj] = ii;

			// printf("%d\t%d\t%f\n", jj, ii, this_element_wav);	// DEBUG

	        }

		for (jj=0; jj<nyelements; jj++) {

			for (ii=nxelements; ii>=0; ii--) {

				this_element_wav = 0.0;
	
				for (kk=0; kk<=polynomial_order; kk++) {
		      
		        		this_element_wav += coeffs[jj][kk]*pow(ii,kk);

		      		}

		      		if (this_element_wav < end_wav) {	// the current index, ii, represents the last pixel with a wavelength <= end_wav. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

		       			break;

		     		}

		   	}

			last_element_index_array[jj] = ii;

			// printf("%d\t%d\t%f\n", jj, ii, this_element_wav);	// DEBUG

	        }
	
		printf("\nWavelength boundaries");
		printf("\n---------------------\n");

		printf("\nInherent minimum wavelength:\t%.2f Å", smallest_wav);
		printf("\nInherent maximum wavelength:\t%.2f Å\n", largest_wav);

		if (start_wav < smallest_wav) {		// Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
		  
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -15, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		} else if (end_wav > largest_wav) {	// Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -16, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		}

		// CALCULATE AND APPLY THROUGHPUT COEFFICIENTS TO TARGET FRAME (ARG 1)
		// ***********************************************************************
		// 1.	Calculate total flux through each fibre between the wavelength
		// 	ranges [start_wav] and [end_wav]

		double disp_val;
	
		double continuum_spectra_totals [nyelements];
		memset(continuum_spectra_totals, 0, sizeof(double)*nyelements);
	
		for(jj=0; jj<nyelements; jj++) {

			// printf("\n%d\t%d", first_element_index_array[jj], last_element_index_array[jj]);		// DEBUG

			disp_val = 0.0;

			for(ii=first_element_index_array[jj]; ii<=last_element_index_array[jj]; ii++) {

				disp_val += cc_ext_cont_frame_values[jj][ii];

			}

			continuum_spectra_totals[jj] = disp_val;
	
		}

		// 2.	Duplicate [continuum_spectra_totals] array and sort into ascending order

		double continuum_spectra_totals_sorted [nyelements];
		memcpy(continuum_spectra_totals_sorted, continuum_spectra_totals, sizeof(double)*nyelements);	

		gsl_sort(continuum_spectra_totals_sorted, 1, nyelements);
		double continuum_spectra_median = gsl_stats_median_from_sorted_data(continuum_spectra_totals_sorted, 1, nyelements);

		// printf("%f\n", continuum_spectra_median);	// DEBUG

		// 3.	Calculate throughput coefficients

		printf("\nFibre fluxes and throughput coefficients");
		printf("\n----------------------------------------\n");
		printf("\nMinimum Flux:\t%.2e\n", continuum_spectra_totals_sorted[0]);
		printf("Maximum Flux:\t%.2e\n", continuum_spectra_totals_sorted[nyelements-1]);
		printf("Median Flux:\t%.2e\n", continuum_spectra_median);
		printf("\nFibre\tTotal Flux\tCoefficient\n\n");

		double throughput_coefficients [nyelements];
		memset(throughput_coefficients, 0, sizeof(double)*nyelements);

		for(jj=0; jj<nyelements; jj++) {

			throughput_coefficients[jj] = continuum_spectra_median / continuum_spectra_totals[jj];

			printf("%d\t%.2e\t%.2f\n", jj+1, continuum_spectra_totals[jj], throughput_coefficients[jj]);

		}

		// 4.	Get cc extracted target fits file (ARG 1) values, apply throughput
		//	corrections and store to 2D array

		double cor_cc_ext_target_frame_values [nyelements][nxelements];
		memset(cor_cc_ext_target_frame_values, 0, sizeof(double)*nxelements*nyelements);

		int this_fibre_index;

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			this_fibre_index = fpixel[1] - 1;

			memset(cc_ext_target_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(cc_ext_target_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, cc_ext_target_f_pixels, NULL, &cc_ext_target_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					cor_cc_ext_target_frame_values[this_fibre_index][ii] = throughput_coefficients[this_fibre_index] * cc_ext_target_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -17, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cc_ext_target_f_status); 

				return 1; 

			}

		}

		// 5.	Set cor_cc_ext_target frame parameters

		fitsfile *cor_cc_ext_target_f_ptr;
	
		int cor_cc_ext_target_f_status = 0;
		long cor_cc_ext_target_f_naxes [2] = {nxelements,FIBRES};
	
		long cor_cc_ext_target_f_fpixel = 1;

		// 6.	Create [cor_cc_ext_target_output_frame_values] array to hold the
		//	output data in the correct format

		double cor_cc_ext_target_output_frame_values [nxelements*nyelements];
		memset(cor_cc_ext_target_output_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (ii=0; ii<nyelements; ii++) {
	
			jj = ii * nxelements;
	
			for (kk=0; kk<nxelements; kk++) {
	
				cor_cc_ext_target_output_frame_values[jj] = cor_cc_ext_target_frame_values[ii][kk];
				jj++;

			}
		
		}

		// 7.	Create and write [cor_cc_ext_target_output_frame_values] to output
		//	file (ARG 5)
	
		if (!fits_create_file(&cor_cc_ext_target_f_ptr, cor_cc_ext_target_f, &cor_cc_ext_target_f_status)) {
	
			if (!fits_create_img(cor_cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, cor_cc_ext_target_f_naxes, &cor_cc_ext_target_f_status)) {

				if (!fits_write_img(cor_cc_ext_target_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], cor_cc_ext_target_f_fpixel, nxelements * nyelements, cor_cc_ext_target_output_frame_values, &cor_cc_ext_target_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -18, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, cor_cc_ext_target_f_status); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -19, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cor_cc_ext_target_f_status); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -20, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cor_cc_ext_target_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Close input files (ARGS 1 and 2), output file (ARG 5) and 
		// [FRARCFIT_OUTPUTF_WAVFITS_FILE] file

		if(fits_close_file(cc_ext_target_f_ptr, &cc_ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -21, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cc_ext_target_f_status); 

			return 1; 

	    	}

		if(fits_close_file(cc_ext_cont_f_ptr, &cc_ext_cont_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -22, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cc_ext_cont_f_status); 

			return 1; 

	    	}

		if(fits_close_file(cor_cc_ext_target_f_ptr, &cor_cc_ext_target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -23, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cor_cc_ext_target_f_status); 

			return 1; 

	    	}

		if (fclose(dispersion_solutions_f)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -24, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		}

		// ***********************************************************************
		// Clean up heap memory

		free(cc_ext_target_f);
		free(cc_ext_cont_f);
		free(cor_cc_ext_target_f);

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", SUCCESS_FLAG, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}
