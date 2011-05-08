/************************************************************************

 File:				frodo_red_trace.c
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
#include "frodo_red_trace.h"
#include "frodo_red_findpeaks_simple_clean.h"

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 4) {

		if(populate_env_variable(FRT_BLURB_FILE, "L2_FRT_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FRT_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -1, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		int order			= strtol(argv[1], NULL, 0);
		int bins			= strtol(argv[2], NULL, 0);	
		int min_rows_per_bin		= strtol(argv[3], NULL, 0);

		// ***********************************************************************
		// Open [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file
	
		FILE *inputfile;
	
		if (!check_file_exists(FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE)) { 

			inputfile = fopen(FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -2, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}

		// ***********************************************************************
		// Find some [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file details

		char input_string [150];

		bool find_ymin_comment = FALSE;
		bool find_ymax_comment = FALSE;
		bool find_ybad_comment = FALSE;

		int y_min, y_max, y_bad;	

		char search_string_1 [15] = "# Min y coord:\0";			// this is the comment to be found from the [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file
		char search_string_2 [15] = "# Max y coord:\0";			// this is another comment to be found from the [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file
		char search_string_3 [12] = "# Bad rows:\0";			// this is another comment to be found from the [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file

		while(!feof(inputfile)) {

			if ((find_ymin_comment == TRUE) && (find_ymax_comment == TRUE) && (find_ybad_comment == TRUE)) {		// have we found what we need?

				break;

			}

			memset(input_string, '\0', sizeof(char)*150);
	
			fgets(input_string, 150, inputfile);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &y_min);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_ymin_comment = TRUE;


			} else if (strncmp(input_string, search_string_2, strlen(search_string_2)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &y_max);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_ymax_comment = TRUE;


			} else if (strncmp(input_string, search_string_3, strlen(search_string_3)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &y_bad);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_ybad_comment = TRUE;


			} 

		}

		if (find_ymin_comment == FALSE) {	// error check - didn't find the comment in the [FRCLEAN_OUTPUTF_PEAKSCLEANSED_FILE] input file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -3, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;

		} else if (find_ymax_comment == FALSE) {	// error check - didn't find the comment in the [FRCLEAN_OUTPUTF_PEAKSCLEANSED_FILE] input file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -4, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;

		} else if (find_ybad_comment == FALSE) {	// error check - didn't find the comment in the [FRCLEAN_OUTPUTF_PEAKSCLEANSED_FILE] input file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -5, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;

		}

		// ***********************************************************************
		// Perform a few checks to ensure the input tracing parameters 
		// are sensible

		if ((order < FRTRACE_VAR_POLYORDER_LO) || (order > FRTRACE_VAR_POLYORDER_HI)) {	// Check [order] is within config limits

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -6, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1; 

		}

		int range = y_max-y_min;

		if (bins > range) { 								// Check there aren't more bins than rows available

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -7, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1; 

		}

		// ***********************************************************************
		// Calculate the bin sizes

		int ii;

		double width_of_bin = (double) range/bins;
		double current_bin_position = (double) y_min;			// store the current bin position [current_bin_position] as a double to avoid rounding errors

		// printf("%f\n", width_of_bin);	// DEBUG
	
		int bin_limits [bins+1];
		memset(bin_limits, 0, sizeof(int)*(bins+1));
	
		for(ii=0; ii<bins+1; ii++) {
	
			bin_limits[ii] = lrint(current_bin_position);		// round [current_bin_position] and store as an integer to use as the current bin limit

			//printf("%f\t%ld\t%d\n", current_bin_position, lrint(current_bin_position), bin_limits[ii]);	// DEBUG

			current_bin_position += width_of_bin;
	
		}

		printf("\nBinning results");
		printf("\n--------------------\n");
		printf("\nMinimum y coordinate:\t%d", y_min);
		printf("\nMaximum y coordinate:\t%d", y_max);
		printf("\nRange:\t\t\t%d\n", range);
		printf("\nNumber of bins:\t\t%d\n", bins);
		printf("Bin width:\t\t%.2f\n", width_of_bin);

		// ***********************************************************************
		// Perform a further check to ensure the input tracing parameters 
		// are sensible

		if ((min_rows_per_bin > bins) || (min_rows_per_bin == 0)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -8, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1; 

		}

		// ***********************************************************************
		// Calculate average spaxel bin coordinates for each fibre

		rewind(inputfile);

		double x_coords_binned [FIBRES][bins];
		memset(x_coords_binned, 0, sizeof(double)*(FIBRES*bins));

		double y_coords_binned [FIBRES][bins];						// we can't just use midpoint of two bins as the majority of the rows may be grouped towards one bin
		memset(y_coords_binned, 0, sizeof(double)*(FIBRES*bins));

		int y_num_rows [FIBRES][bins];
		memset(y_num_rows, 0, sizeof(int)*(FIBRES*bins));	

		double coord_x;
		int fibre_number, coord_y;

		int this_bin = 0;

		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*150);
	
			fgets(input_string, 150, inputfile);	

			if (strtol(&input_string[0], NULL, 0) > 0) {					// else check the line begins with a positive number (usable)
	
				sscanf(input_string, "%d\t%lf\t%d\t", &fibre_number, &coord_x, &coord_y);
	
				if (coord_y > bin_limits[this_bin+1]) {
		
					this_bin++;

				} 

				y_num_rows [fibre_number-1][this_bin]++;			// this is the number of rows in each bin (n.b. for all fibres)
				x_coords_binned[fibre_number-1][this_bin] += coord_x;		// this is the cumulative x coordinate for each bin/fibre
				y_coords_binned[fibre_number-1][this_bin] += (double) coord_y;	// this is the cumulative y coordinate for each bin/fibre

			}

		}

		int jj;

		for (ii=0; ii<FIBRES; ii++) {

			for (jj=0; jj<bins; jj++) {

				if ((y_num_rows[ii][jj]) != 0) {	// are there any rows for this fibre/bin?

					x_coords_binned[ii][jj] /= (double) y_num_rows[ii][jj];		// must divide by [y_num_rows] to get average coordinate 
					y_coords_binned[ii][jj] /= (double) y_num_rows[ii][jj];		// must divide by [y_num_rows] to get average coordinate

				}

				// printf("%f\t%f\t%d\n", x_coords_binned[ii][jj], y_coords_binned[ii][jj], y_num_rows[ii][jj]);	// DEBUG		

			}

		}

		// ***********************************************************************
		// Create [FRTRACE_OUTPUTF_TRACES_FILE] output file and print a few 
		// parameters

		FILE *outputfile;
		outputfile = fopen(FRTRACE_OUTPUTF_TRACES_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -9, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;


		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", FRTRACE_OUTPUTF_TRACES_FILE);
		fprintf(outputfile, "# List of fibre tramline trace coefficients and corresponding chi-squareds found using the frtrace program.\n\n");
		fprintf(outputfile, "# Run datetime:\t\t%s\n\n", timestr);
		fprintf(outputfile, "# Polynomial Order:\t%d\n\n", order);

		// ***********************************************************************
		// Fit each fibre and store results to [FRTRACE_OUTPUTF_TRACES_FILE] file

		double coords_x [bins];
		double coords_y [bins];

		double coeffs [order+1];

		double this_chi_squared, chi_squared_min = 0.0, chi_squared_max = 0.0, chi_squared = 0.0;

		int bins_deficit;	// var to store how many bins have zero rows
		int bins_used;		// var to store how many bins we're using

		for (ii=0; ii<FIBRES; ii++) {

			bins_deficit = 0;
			bins_used = 0;

			memset(coords_y, 0, sizeof(double)*(bins));
			memset(coeffs, 0, sizeof(double)*(order+1));
			memset(coords_x, 0, sizeof(double)*(bins));	

			for (jj=0; jj<bins; jj++) {

				if (y_num_rows[ii][jj] < min_rows_per_bin) { 	// are there enough [min_rows_per_bin] rows for this fibre/bin?

					bins_deficit++;

				} else {

					coords_x[jj-bins_deficit] = x_coords_binned[ii][jj];
					coords_y[jj-bins_deficit] = y_coords_binned[ii][jj];

					bins_used++;
	
				}

			}

			if (calc_least_sq_fit(order, bins-bins_deficit, coords_y, coords_x, coeffs, &this_chi_squared)) {	// reversed [coord_y] and [coord_x] as want to find x = f(y) not y = f(x)

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -10, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

				fclose(inputfile);
				fclose(outputfile);

				return 1; 

			}

			fprintf(outputfile, "%d\t", ii+1);

			for (jj=0; jj<=order; jj++) {

				fprintf(outputfile, FRTRACE_VAR_ACCURACY_COEFFS, coeffs[jj]);
				fprintf(outputfile, "\t");

			}

			fprintf(outputfile, FRTRACE_VAR_ACCURACY_CHISQ, this_chi_squared);
			fprintf(outputfile, "\n");

			if ((ii==0) || (this_chi_squared < chi_squared_min)) { 		// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				chi_squared_min = this_chi_squared;

			} else if ((ii==0) || (this_chi_squared > chi_squared_max)) {	// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				chi_squared_max = this_chi_squared;

			}

			// printf("%d\t%f\n", ii, this_chi_squared);	// DEBUG

			chi_squared += this_chi_squared;

		}

		fprintf(outputfile, "%d", EOF);

		printf("\nFitting results");
		printf("\n--------------------\n");
		printf("\nNumber of bins used:\t%d\n", bins_used);
		printf("\nMin χ2:\t\t\t%.2f\n", chi_squared_min);
		printf("Max χ2:\t\t\t%.2f\n", chi_squared_max);
		printf("Average χ2:\t\t%.2f\n", chi_squared/FIBRES);

		// ***********************************************************************
		// Perform a few checks to ensure the chi squareds are sensible 

		if ((chi_squared_min < FRTRACE_VAR_CHISQUARED_MIN) || (chi_squared_max > FRTRACE_VAR_CHISQUARED_MAX)) {	// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

			RETURN_FLAG = 2;

		}

		// ***********************************************************************
		// Close [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] input file and
		// [FRTRACE_OUTPUTF_TRACES_FILE] output file

		if (fclose(inputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -11, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(outputfile);

			return 1; 

		}

		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", -12, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATTR", RETURN_FLAG, "Status flag for L2 frtrace routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}
