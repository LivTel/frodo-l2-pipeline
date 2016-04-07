#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gsl/gsl_statistics_double.h>

#include "frodo_error_handling.h"
#include "frodo_functions.h"
#include "frodo_config.h"
#include "frodo_red_extract.h"
#include "frodo_red_trace.h"

int main ( int argc, char *argv [] ) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if ( argc != 8 ) {

		if(populate_env_variable(FRE_BLURB_FILE, "L2_FRE_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FRE_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -1, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else { 

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f			= strdup(argv[1]);
		int half_aperture_num_pix	= strtol(argv[2], NULL, 0);
		int flip_disp_axis		= strtol(argv[3], NULL, 0);
		int optimal			= strtol(argv[4], NULL, 0);
		int cr_reject_sigma		= strtol(argv[5], NULL, 0);
		char *weights_f			= strdup(argv[6]);
		char *output_f			= strdup(argv[7]);

		// ***********************************************************************
		// Is there EXPLICIT cross-talk? i.e. Is the extraction aperture width >
		// the average fibre profile width?

		if (half_aperture_num_pix > (FREXTRACT_VAR_FIBRE_PROFILE_WIDTH - 1) / 2) {

			RETURN_FLAG = 2;

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

					free(input_f);
					free(weights_f);
					free(output_f);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -3, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				free(input_f);
				free(weights_f);
				free(output_f);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -4, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status); 

			free(input_f);
			free(weights_f);
			free(output_f);

			return 1; 

		}

		// ***********************************************************************
		// Get READNOIS and GAIN keywords (required for variance estimation)

		int input_f_status_read_key = 0;
		float readnoise;
		if (fits_read_key(input_f_ptr, TFLOAT, "READNOIS", &readnoise, NULL, &input_f_status_read_key)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -17, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status_read_key); 

			free(input_f);
			free(weights_f);
			free(output_f);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}
		float gain;
		if (fits_read_key(input_f_ptr, TFLOAT, "GAIN", &gain, NULL, &input_f_status_read_key)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -17, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status_read_key); 

			free(input_f);
			free(weights_f);
			free(output_f);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

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

			free(input_f);
			free(weights_f);
			free(output_f);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

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
				find_polynomialorder_comment = TRUE;
				break;


			} 

		}

		if (find_polynomialorder_comment == FALSE) {	// error check - didn't find the comment in the [FRTRACE_OUTPUTF_TRACES_FILE] file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -6, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(weights_f);
			free(output_f);
			fclose(traces_file);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

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

		double input_frame_values [FIBRES][nyelements];
		if (optimal == TRUE) {

			// ***********************************************************************
			// Open weights.dat file (only if optimal extraction flag set)
	
			FILE *weights_file;
	
			if (!check_file_exists(weights_f)) { 

				weights_file = fopen(weights_f , "r");

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -15, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

				free(input_f);
				free(weights_f);
				free(output_f);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1;

			}

			// PARSE WEIGHTS FOR OPTIMAL EXTRACTION
			// ***********************************************************************
			// Find some weights.dat file details

			find_polynomialorder_comment = FALSE;

			int polynomial_order_2;	

			char search_string_2 [20] = "# Polynomial Order:\0";	// this is the comment to be found from the weights.dat file

			while(!feof(weights_file)) {

				memset(input_string, '\0', sizeof(char)*500);
	
				fgets(input_string, 500, weights_file);	

				if (strncmp(input_string, search_string_2, strlen(search_string_2)) == 0) { 

					sscanf(input_string, "%*[^\t]%d", &polynomial_order_2);		// read all data up to tab as string ([^\t]), but do not store (*)
					find_polynomialorder_comment = TRUE;
					break;


				} 

			}

			if (find_polynomialorder_comment == FALSE) {	// error check - didn't find the comment in the weights.dat file

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -16, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

				free(input_f);
				free(weights_f);
				free(output_f);
				fclose(traces_file);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1;

			}

			// ***********************************************************************
			// Rewind and extract coefficients from weights.dat file 

			rewind(weights_file);

			double coeffs_2 [FIBRES][polynomial_order_2+1];
			memset(coeffs_2, 0, sizeof(double)*FIBRES*(polynomial_order_2+1));

			while(!feof(weights_file)) {

				memset(input_string, '\0', sizeof(char)*500);
	
				fgets(input_string, 500, weights_file);

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

						} else if ((token_index >= 1) && (token_index <= polynomial_order_2+1)) { 	// coeff token

							this_coeff = strtod(token, NULL);
							// printf("%d\t%d\t%e\n", this_fibre, coeff_index, this_coeff);		// DEBUG
							coeffs_2[this_fibre-1][coeff_index] = this_coeff;
							coeff_index++;

						} else if (token_index == polynomial_order_2+2) {				// chisquared token

							this_chisquared = strtod(token, NULL);

						}

						token_index++;

					}

				}
	
			}

			// ***********************************************************************
			// Optimal aperture extraction of flux using coefficients from
			// [FRTRACE_OUTPUTF_TRACES_FILE] file and weights from weights.dat file

			int ii, jj;
			int y_int;

			memset(input_frame_values, 0, sizeof(double)*nyelements*FIBRES);

			int n_reject = 0;
			for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

				memset(input_f_pixels, 0, sizeof(double)*nxelements);

				if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

					y_int = fpixel[1];

					for (ii=0; ii<FIBRES; ii++) {

						double x = 0.0;

						// ***********************************************************************
						// Find tracing centroid from polynomial coefficients
	
						for (jj=0; jj<polynomial_order+1; jj++) {
					
							x += (coeffs[ii][jj]*(pow(y_int,jj)));
	
						}

						x -= INDEXING_CORRECTION;

						// ***********************************************************************
						// Calculate sigma for this lambda from weights.dat 

						double sigma = 	0.0;
						double fwhm =	0.0;
						for (jj=0; jj<polynomial_order_2+1; jj++) {
					
							fwhm += (coeffs_2[ii][jj]*(pow(y_int,jj)));
	
						}

						sigma = fwhm/2.35;		// file values are FWHM

						// ***********************************************************************
						// 1.	Input check. Does [x] violate the img boundaries?

						if ((x > nxelements) || (x <= 0)) {

							write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -7, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
							fits_report_error(stdout, input_f_status); 

							free(input_f);
							free(weights_f);
							free(output_f);
							fclose(traces_file);
							if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

							return 1;

						}

						// ***********************************************************************
						// 2.	Construct model for spatial profile (P) as a function of x for 
						//	each wavelength.
						//
						//	This is done by averaging [smoothing_half_aperture_num_pix] pixels  
						//	in the dispersion direction to negate the influence of cosmic rays. 

						// this is sliding aperture, taking into account row boundaries

						int smoothing_half_aperture_num_pix = 3;

						int y_int_this_aperture_lo = y_int-smoothing_half_aperture_num_pix;
						int y_int_this_aperture_hi = y_int+smoothing_half_aperture_num_pix;
						int d;
						if (y_int_this_aperture_lo < cut_y[0]) {
							d = cut_y[0]-y_int_this_aperture_lo;
							y_int_this_aperture_lo+=d;
							y_int_this_aperture_hi+=d;
						} else if (y_int_this_aperture_hi > cut_y[1]) {
							d = y_int_this_aperture_hi-cut_y[1];
							y_int_this_aperture_lo-=d;
							y_int_this_aperture_hi-=d;
						}

						int y_aper_idx = 0;
						int nxelements_this_aperture = 1+(half_aperture_num_pix*2);
						double input_f_pixels_this_aperture[nxelements_this_aperture];
						for (y_aper_idx=y_int_this_aperture_lo; y_aper_idx<=y_int_this_aperture_hi; y_aper_idx++) {
							// find the centroid of the aperture
							double x_this_aperture = 0.0;
							for (jj=0; jj<polynomial_order+1; jj++) {
								x_this_aperture += (coeffs[ii][jj]*(pow(y_int,jj)));
							}

							int x_lo = round(x_this_aperture-half_aperture_num_pix);
							int x_hi = round(x_this_aperture+half_aperture_num_pix);	

							// extract values for this aperture
							memset(input_f_pixels_this_aperture, 0, sizeof(double)*nxelements_this_aperture);
							long fpixel_this_aperture[2] = {x_lo, y_aper_idx};
							if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel_this_aperture, nxelements_this_aperture, NULL, input_f_pixels_this_aperture, NULL, &input_f_status)) {
								//int aa=0; for (aa=0; aa<nxelements_this_aperture; aa++) printf("%f\t", input_f_pixels_this_aperture[aa]); printf("\n");
							} else {

								write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -8, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
								fits_report_error(stdout, input_f_status); 

								free(input_f);
								free(weights_f);
								free(output_f);
								fclose(traces_file);
								if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

								return 1; 

							}

						}

						printf("\n");
	
						double x_low, x_high;

						x_low = x-half_aperture_num_pix-0.5;		// -0.5 as trace assumes centroid of pixel is integer (e.g. 100 instead of 100.5)
						x_high = x+half_aperture_num_pix+1-0.5;		

						long x_low_floor, x_high_floor;
	
						x_low_floor = floor(x_low);
						x_high_floor = floor(x_high);
			

						double total_flux = 0;

						for (jj=x_low_floor; jj<=x_high_floor; jj++) {
							double partial_fraction_of_bin;
							if (jj == x_low_floor) {					// outside pixel where partial flux needs to be taken into account
								partial_fraction_of_bin = (x_low_floor + 1) - x_low;
							} else if (jj == x_high_floor) {				// outside pixel where partial flux needs to be taken into account
								partial_fraction_of_bin = x_high - x_high_floor;
							} else {
								partial_fraction_of_bin = 1;
							}
							total_flux += partial_fraction_of_bin * input_f_pixels[jj];
						}


						// ***********************************************************************
						// 3.	Cosmic ray filter (see Horne, 1986)
						//	
						//	The procedure for this is as follows. 
						//	a) cycle through each pixel in the aperture
						//	b) for each pixel, calculate how much of the flux would be
						//	   expected by taking the product of the weight (w) and the total
						//	   flux in the aperture.
						//	c) compare the squared difference between this product and that of 
						//	   the observed flux, and calculate how much this difference is a 
						//	   multiple of the pixel variance. 
						//
						//	The output is a boolean array, designating each pixel in the
						//	aperture as either rejected (true) or not (false).

						int n_px = x_high_floor-x_low_floor+1;
						bool jj_reject[n_px];
						int kk;
						for (kk=0; kk<n_px; kk++) { jj_reject[kk] = false; }	// init all elements to 0
						if (cr_reject_sigma != 0) {
							int n_iter = 0;
							bool is_iterating = true;
							while (is_iterating) {				// keep iterating until no more pixels have been rejected
								n_iter++;
								int idx = 0;
								for (jj=x_low_floor; jj<=x_high_floor; jj++) {
									// skip this pixel if it's been rejected in a previous iteration
									if (jj_reject[idx] == true) {
										idx++;
										continue;
									}
								
									// establish the actual/expected flux and error
									// n.b. [area] is required to keep track of the true fraction of flux expected when pixels have been removed
									double partial_fraction_of_bin, w, area=1;
									if (jj == x_low_floor) {
										partial_fraction_of_bin = (x_low_floor + 1) - x_low;
										w = get_gaussian_integral_between_two_limits(x, sigma, x_low-0.5, (x_low_floor + 1)-0.5);
									} else if (jj == x_high_floor) {
										partial_fraction_of_bin = x_high - x_high_floor;
										w = get_gaussian_integral_between_two_limits(x, sigma, x_high_floor-0.5, x_high-0.5);
									} else {
										partial_fraction_of_bin = 1;
										w = get_gaussian_integral_between_two_limits(x, sigma, jj-0.5, jj+1-0.5);
									}
									double this_obs_value = (partial_fraction_of_bin*input_f_pixels[jj]);				// obs
									double this_exp_value = ((w/area)*total_flux);							// exp (model), 
									double diff = this_obs_value - this_exp_value;							// obs - exp 
									double diff_sq = pow(diff, 2);
									double this_pixel_value = partial_fraction_of_bin*input_f_pixels[jj];
									double var = pow(readnoise, 2) + (fabs(this_pixel_value*w)/gain);				// eq 12 of Horne (1986)
									double obs_sigma = diff_sq/var;									// observed multiple of variance

									// flag this pixel if required, decrement total_flux and area accordingly
									if (obs_sigma > cr_reject_sigma) {
										area -= w;			// remove this pixel's weight from the normed [area]
										total_flux -= this_obs_value;	// remove this pixel's value from the [total_flux]
										jj_reject[idx] = true;		// flag this pixel as a CR
										n_reject++;
										break;
									}
									idx++;

								}
		
								// have we reached final pixel without setting the rejection flag in this iteration?
								if (idx == n_px) {
									is_iterating = false;
								}
							} 
						}
				
						// ***********************************************************************
						// 3.	Populate [input_frame_values] array with extracted and weighted 
						//	flux, accounting for CRs

						int idx = 0;
						double numerator = 0, denominator = 0;
						for (jj=x_low_floor; jj<=x_high_floor; jj++) {
							double partial_fraction_of_bin, w;
							if (jj == x_low_floor) {
								partial_fraction_of_bin = (x_low_floor + 1) - x_low;
								w = get_gaussian_integral_between_two_limits(x, sigma, x_low-0.5, (x_low_floor + 1)-0.5);
							} else if (jj == x_high_floor) {
								partial_fraction_of_bin = x_high - x_high_floor;
								w = get_gaussian_integral_between_two_limits(x, sigma, x_high_floor-0.5, x_high-0.5);
							} else {
								partial_fraction_of_bin = 1;
								w = get_gaussian_integral_between_two_limits(x, sigma, jj-0.5, jj+1-0.5);
							}	

							double this_pixel_value = partial_fraction_of_bin*input_f_pixels[jj];
							double var = pow(readnoise, 2) + (fabs(this_pixel_value*w)/gain);		// eq 12 of Horne (1986)

							if (jj_reject[idx] == false) {
								//numerator += this_pixel_value;	// DEBUG (no weight)
								//denominator = 1;			// DEBUG (no weight)
								numerator += ((this_pixel_value*w)/var);
								denominator += (pow(w, 2)/var);
							}
							
							idx++;
						}

						input_frame_values[ii][y_int-1] = numerator/denominator;

						/*if (ii==11 && y_int==896) {	// DEBUG
							for (jj=x_low_floor; jj<=x_high_floor; jj++) { printf("%f\t", input_f_pixels[jj]);}
							for (kk=0; kk<n_px; kk++) { printf("%d\t", jj_reject[kk]);}
							printf("\n");
						}*/
	

					}

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -8, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					free(input_f);
					free(weights_f);
					free(output_f);
					fclose(traces_file);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

			}	

			printf("\nCosmic ray rejection");
			printf("\n--------------------\n");
			printf("\nNumber of cosmics rejected:\t%d\n", n_reject);
		
		} else {

			// ***********************************************************************
			// Standard aperture extraction of flux using coefficients from
			// [FRTRACE_OUTPUTF_TRACES_FILE] file

			int ii, jj;

			int y_int;

			memset(input_frame_values, 0, sizeof(double)*nyelements*FIBRES);

			//double total_flux = 0;	// DEBUG
			// int count = 0;		// DEBUG

			for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

				memset(input_f_pixels, 0, sizeof(double)*nxelements);

				if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

					y_int = fpixel[1];

					for (ii=0; ii<FIBRES; ii++) {

						double x = 0.0;

						// ***********************************************************************
						// Find tracing centroid from polynomial coefficients
	
						for (jj=0; jj<polynomial_order+1; jj++) {
					
							x += (coeffs[ii][jj]*(pow(y_int,jj)));
	
						}

						x -= INDEXING_CORRECTION;

						// ***********************************************************************
						// Does [x] violate the img boundaries?

						if ((x > nxelements) || (x <= 0)) {

							write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -7, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
							fits_report_error(stdout, input_f_status); 

							free(input_f);
							free(weights_f);
							free(output_f);
							fclose(traces_file);
							if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

							return 1;

						}

						// ***********************************************************************
						// Extract flux within aperture

						double x_low, x_high;

						x_low = x-half_aperture_num_pix-0.5;
						x_high = x+half_aperture_num_pix+0.5;			

						int x_low_floor, x_high_floor;
	
						x_low_floor = floor(x-half_aperture_num_pix-0.5);
						x_high_floor = floor(x+half_aperture_num_pix+0.5);
			
						for (jj=x_low_floor; jj<=x_high_floor; jj++) {

							if (jj == x_low_floor) {			// outside pixel where partial flux needs to be taken into account

								double partial_fraction_of_bin = (x_low_floor + 1) - x_low;
								input_frame_values[ii][y_int-1] += partial_fraction_of_bin * input_f_pixels[jj];
								// total_flux += input_f_pixels[jj];	// DEBUG

							} else if (jj == x_high_floor) {		// outside pixel where partial flux needs to be taken into account

								double partial_fraction_of_bin = x_high - x_high_floor;
								input_frame_values[ii][y_int-1] += partial_fraction_of_bin * input_f_pixels[jj];
								// total_flux += input_f_pixels[jj];	// DEBUG

							} else {

								input_frame_values[ii][y_int-1] += input_f_pixels[jj];
								// total_flux += input_f_pixels[jj];	// DEBUG
				
							}

						}

					}

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -8, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					free(input_f);
					free(weights_f);
					free(output_f);
					fclose(traces_file);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

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

                int ii, jj;

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

					free(input_f);
					free(weights_f);
					free(output_f);
					fclose(traces_file);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
					if(fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status);

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -10, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(weights_f);
				free(output_f);
				fclose(traces_file);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if(fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status);

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -11, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			free(input_f);
			free(weights_f);
			free(output_f);
			fclose(traces_file);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Free arrays on heap

		free(input_f);
		free(weights_f);
		free(output_f);

		// ***********************************************************************
		// Close [FRTRACE_OUTPUTF_TRACES_FILE] file, input file (ARG 1) and output 
		// file (ARG 4)

		if (fclose(traces_file)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -12, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		}

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -13, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, input_f_status); 

			return 1; 

	    	}

		if(fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", -14, "Status flag for L2 frextract routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATEX", RETURN_FLAG, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

