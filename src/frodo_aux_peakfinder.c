#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "frodo_error_handling.h"
#include "frodo_functions.h"
#include "frodo_config.h"
#include "frodo_aux_peakfinder.h"

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc !=9) {

		if(populate_env_variable(FAP_BLURB_FILE, "L2_FAP_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FAP_BLURB_FILE);

		}
	
		return 1;

	} else {

		char time_start [80];
		memset(time_start, '\0', sizeof(char)*80);

		find_time(time_start);

		// ***********************************************************************
		// Redefine routine input parameters

		char *cont_f			= strdup(argv[1]);
		int min_dist			= strtol(argv[2], NULL, 0);
		int half_aperture_num_pix	= strtol(argv[3], NULL, 0);
		int der_tol			= strtol(argv[4], NULL, 0);
		int der_tol_ref_px		= strtol(argv[5], NULL, 0);
		int min_rows			= strtol(argv[6], NULL, 0);
		int min_peaks			= strtol(argv[7], NULL, 0);
		int max_peaks			= strtol(argv[8], NULL, 0);

		// ***********************************************************************
		// Open cont file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *cont_f_ptr;

		int cont_f_maxdim = 2, cont_f_status = 0, cont_f_bitpix, cont_f_naxis;
		long cont_f_naxes [2] = {1,1};

		if(!fits_open_file(&cont_f_ptr, cont_f, IMG_READ_ACCURACY, &cont_f_status)) {

			if(!populate_img_parameters(cont_f, cont_f_ptr, cont_f_maxdim, &cont_f_bitpix, &cont_f_naxis, cont_f_naxes, &cont_f_status, "ARC FRAME")) {

				if (cont_f_naxis != 2) {	// any data format checks here
	
					free(cont_f);
					if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

					return 1;
	
				}

			} else { 

				fits_report_error(stdout, cont_f_status); 

				free(cont_f);
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 

			}

		} else { 

			fits_report_error(stdout, cont_f_status); 

			free(cont_f);

			return 1; 

		}

		// ***********************************************************************
		// Set the range limits

		int cut_x [2] = {1, cont_f_naxes[0]};
		int cut_y [2] = {1, cont_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from arc fits file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from arc fits file (ARG 1)

		double cont_f_pixels [nxelements];

		// ***********************************************************************
		// Get arc fits file (ARG 1) values and store in 2D array

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

				fits_report_error(stdout, cont_f_status); 

				free(cont_f);
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Find number of rows for this derivative tolerance

		int jj;

		int peaks [nxelements];

		int this_num_peaks;

		int this_derivative_tol = der_tol;

		double row_values [nxelements];

		int this_derivative_tol_rows = 0;

		for (jj=0; jj<nyelements; jj++) {		// cycle all rows in frame

			memset(row_values, 0, sizeof(double)*nxelements);
			memset(peaks, 0, sizeof(int)*nxelements);
	
			for (ii=0; ii<nxelements; ii++) {	// populate row_values array for this row

					row_values[ii] = cont_frame_values[jj][ii];
	
			}

			find_peaks(nxelements, row_values, peaks, &this_num_peaks, min_dist, half_aperture_num_pix, this_derivative_tol, der_tol_ref_px, INDEXING_CORRECTION);	

			if ((this_num_peaks <= max_peaks) && (this_num_peaks >= min_peaks)) { // satisfies row selection criteria?

				this_derivative_tol_rows++;				// ..increment counter

			}

		}

		printf("\nPeak finding results");
		printf("\n--------------------\n");
		printf("\nFound %d rows (%d required) with between %d and %d peaks using a tolerance of %d.\n", this_derivative_tol_rows, min_rows, min_peaks, max_peaks, this_derivative_tol);

		if (this_derivative_tol_rows < min_rows) {	// break if there are too few rows found

			free(cont_f);
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1;

		} 

		// ***********************************************************************
		// Create [FRPEAKFINDER_OUTPUTF_PEAKS_FILE] output file and print a few 
		// parameters

		FILE *outputfile;
		outputfile = fopen(FRPEAKFINDER_OUTPUTF_PEAKS_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			free(cont_f);
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1;


		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", FRPEAKFINDER_OUTPUTF_PEAKS_FILE);
		fprintf(outputfile, "# Lists the coordinates of the peaks found using the frfind routine.\n\n");
		fprintf(outputfile, "# Run filename:\t%s\n", cont_f);
		fprintf(outputfile, "# Run datetime:\t%s\n\n", timestr);
		fprintf(outputfile, "# Rows found:\t%d\n", this_derivative_tol_rows);

		// ***********************************************************************
		// Execute code to find centroids and store to 
		// [FRPEAKFINDER_OUTPUTF_PEAKS_FILE] output file

		double peak_centroids [max_peaks];

		int count = 1;

		for (jj=0; jj<nyelements; jj++) {		// cycle all rows in frame

			memset(row_values, 0, sizeof(double)*nxelements);
			memset(peak_centroids, 0, sizeof(double)*max_peaks);

			for (ii=0; ii<nxelements; ii++) {	// populate row_values array for this row

				row_values[ii] = cont_frame_values[jj][ii];
	
			}

			find_peaks(nxelements, row_values, peaks, &this_num_peaks, min_dist, half_aperture_num_pix, this_derivative_tol, der_tol_ref_px, INDEXING_CORRECTION);

			if (find_centroid_parabolic(row_values, peaks, this_num_peaks, peak_centroids, INDEXING_CORRECTION)) {

				free(cont_f);
				fclose(outputfile);
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 		

			}

			if ((this_num_peaks <= max_peaks) && (this_num_peaks >= min_peaks)) {	// satisfies row selection criteria	

				fprintf(outputfile, "\n# ROW:\t%d\n\n", count);		// ..output the row information to file

				for (ii=0; ii<this_num_peaks; ii++) {

					fprintf(outputfile, "%d\t", ii+1);
					fprintf(outputfile, "%.2f\t", peak_centroids[ii]);		// x
					fprintf(outputfile, "%d\t\n", jj+1);				// y

				}
				
			count++;

			}

		}

		fprintf(outputfile, "%d", EOF);

		// ***********************************************************************
		// Clean up heap memory

		free(cont_f);

		// ***********************************************************************
		// Close [FRPEAKFINDER_OUTPUTF_PEAKS_FILE] output file and arc file 
		// (ARG 1)

		if (fclose(outputfile)) {

			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1; 

		}

		if(fits_close_file(cont_f_ptr, &cont_f_status)) { 

			fits_report_error (stdout, cont_f_status); 

			return 1; 

	    	}

		return 0;

	}

}


