/************************************************************************

 File:				frodo_red_findpeaks_simple.c
 Last Modified Date:     	01/05/11

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
#include "frodo_red_findpeaks_simple.h"

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 8) {

		if(populate_env_variable(FRFS_BLURB_FILE, "L2_FRFS_BLURB_FILE")) {

			RETURN_FLAG = 1;

		}

		print_file(FRFS_BLURB_FILE);

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -1, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *cont_f			= strdup(argv[1]);
		int min_dist			= strtol(argv[2], NULL, 0);
		int half_aperture_num_pix	= strtol(argv[3], NULL, 0);
		int der_tol_min			= strtol(argv[4], NULL, 0);
		int der_tol_max			= strtol(argv[5], NULL, 0);
		int der_tol_ref_px		= strtol(argv[6], NULL, 0);
		int min_rows			= strtol(argv[7], NULL, 0);

		// ***********************************************************************
		// Open cont file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *cont_f_ptr;

		int cont_f_maxdim = 2, cont_f_status = 0, cont_f_bitpix, cont_f_naxis;
		long cont_f_naxes [2] = {1,1};

		if(!fits_open_file(&cont_f_ptr, cont_f, IMG_READ_ACCURACY, &cont_f_status)) {

			if(!populate_img_parameters(cont_f, cont_f_ptr, cont_f_maxdim, &cont_f_bitpix, &cont_f_naxis, cont_f_naxes, &cont_f_status, "CONTINUUM FRAME")) {

				if (cont_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -2, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -3, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cont_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -4, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cont_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Set the range limits

		int cut_x [2] = {1, cont_f_naxes[0]};
		int cut_y [2] = {1, cont_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from continuum fits file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from continuum fits file (ARG 1)

		double cont_f_pixels [nxelements];

		// ***********************************************************************
		// Get continuum fits file (ARG 1) values and store in 2D array

		int ii;

		double cont_frame_values [nyelements][nxelements];
		memset(cont_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(cont_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(cont_f_ptr, TDOUBLE, fpixel, nxelements, NULL, cont_f_pixels, NULL, &cont_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					cont_frame_values[fpixel[1]-1][ii] = cont_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -5, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cont_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Find best derivative tolerance between limits

		int jj;

		int peaks [nxelements];

		int this_num_peaks;

		int this_derivative_tol, this_derivative_tol_rows;
		int best_derivative_tol, best_derivative_tol_rows;

		double row_values [nxelements];

		//printf("Best Tolerance\tNumber of rows\n\n");		// DEBUG

		for (this_derivative_tol=der_tol_min; this_derivative_tol<=der_tol_max; this_derivative_tol++) { 

			this_derivative_tol_rows = 0;

			for (jj=0; jj<nyelements; jj++) {		// cycle all rows in frame

				memset(row_values, 0, sizeof(double)*nxelements);
				memset(peaks, 0, sizeof(int)*nxelements);
	
				for (ii=0; ii<nxelements; ii++) {	// populate row_values array for this row

						row_values[ii] = cont_frame_values[jj][ii];
	
				}

				find_peaks(nxelements, row_values, peaks, &this_num_peaks, min_dist, half_aperture_num_pix, this_derivative_tol, der_tol_ref_px, INDEXING_CORRECTION);	

				if (this_num_peaks == FIBRES) {				// if we find the correct number of peaks..

					this_derivative_tol_rows++;			// ..increment counter

				}

			}

			if ((this_derivative_tol==der_tol_min) || (this_derivative_tol_rows > best_derivative_tol_rows)) {

				best_derivative_tol = this_derivative_tol;
				best_derivative_tol_rows = this_derivative_tol_rows;

				//printf("%d\t\t%d\n", best_derivative_tol, best_derivative_tol_rows);		// DEBUG

			}

		}

		printf("\nPeak finding results");
		printf("\n--------------------\n");
		printf("\nFound %d rows (%d required) with %d peaks using a tolerance of %d.\n", best_derivative_tol_rows, min_rows, FIBRES, best_derivative_tol);

		if (best_derivative_tol_rows < min_rows) {	// break if there are too few rows found

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -6, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			return 1;

		} 

		// ***********************************************************************
		// Create [FRFIND_OUTPUTF_PEAKS_FILE] output file and print a few 
		// parameters

		FILE *outputfile;
		outputfile = fopen(FRFIND_OUTPUTF_PEAKS_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -7, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			return 1;


		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", FRFIND_OUTPUTF_PEAKS_FILE);
		fprintf(outputfile, "# Lists the coordinates of the peaks found using the frfind routine.\n\n");
		fprintf(outputfile, "# Run filename:\t%s\n", cont_f);
		fprintf(outputfile, "# Run datetime:\t%s\n\n", timestr);
		fprintf(outputfile, "# Rows found:\t%d\n", best_derivative_tol_rows);

		// ***********************************************************************
		// Execute code with best derivative tolerance value to find centroids 
		// and store to [FRFIND_OUTPUTF_PEAKS_FILE] output file

		double peak_centroids [FIBRES];

		int count = 1;

		for (jj=0; jj<nyelements; jj++) {		// cycle all rows in frame

			memset(row_values, 0, sizeof(double)*nxelements);
			memset(peak_centroids, 0, sizeof(double)*FIBRES);

			for (ii=0; ii<nxelements; ii++) {	// populate row_values array for this row

				row_values[ii] = cont_frame_values[jj][ii];
	
			}

			find_peaks(nxelements, row_values, peaks, &this_num_peaks, min_dist, half_aperture_num_pix, best_derivative_tol, der_tol_ref_px, INDEXING_CORRECTION);

			if (find_centroid_parabolic(row_values, peaks, this_num_peaks, peak_centroids, INDEXING_CORRECTION)) {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -8, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

				return 1; 		

			}

			if (this_num_peaks == FIBRES) {					// if we find the correct number of peaks..	

				fprintf(outputfile, "\n# ROW:\t%d\n\n", count);		// ..output the row information to file

				for (ii=0; ii<FIBRES; ii++) {

					fprintf(outputfile, "%d\t", ii+1);
					fprintf(outputfile, "%.2f\t", peak_centroids[ii]);		// x
					fprintf(outputfile, "%d\t\n", jj+1);				// y

				}
				
			count++;

			}

		}

		fprintf(outputfile, "%d", EOF);

		// ***********************************************************************
		// Close [FRFIND_OUTPUTF_PEAKS_FILE] output file and continuum file 
		// (ARG 1)

		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -9, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			return 1; 

		}

		if(fits_close_file(cont_f_ptr, &cont_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -10, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cont_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Clean up heap memory

		free(cont_f);

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", RETURN_FLAG, "Status flag for L2 frfind routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

		return 0;

	}

}


