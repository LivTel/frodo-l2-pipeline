/************************************************************************

 File:				frodo_red_findpeaks_simple_clean.c
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
#include "frodo_red_findpeaks_simple_clean.h"
#include "frodo_red_findpeaks_simple.h"

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 3) {

		if(populate_env_variable(FRFSC_BLURB_FILE, "L2_FRFSC_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(FRFSC_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -1, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		double pixel_tolerance		= strtod(argv[1], NULL);
		int max_bad_rows		= strtol(argv[2], NULL, 0);
	
		// ***********************************************************************
		// Open [FRFIND_OUTPUTF_PEAKS] input file
	
		FILE *inputfile;
	
		if (!check_file_exists(FRFIND_OUTPUTF_PEAKS_FILE)) { 

			inputfile = fopen(FRFIND_OUTPUTF_PEAKS_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -2, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}

		// ***********************************************************************
		// Use the [FRFIND_OUTPUTF_PEAKS_FILE] to calculate the average x
		// coordinate for each fibre and find [y_min]/[y_max]

		char input_string [150];

		double array_coord_x_av [FIBRES];
		memset(array_coord_x_av, 0, sizeof(double) * FIBRES);

		int fibre_number;
		double coord_x;

		int coord_y;

		int y_min, y_max;	

		int rows;

		char search_string_1 [7] = "# ROW:\0";	// this is the comment to be found from the [FRFIND_OUTPUTF_PEAKS_FILE] input file
	
		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*150);

			fgets(input_string, 150, inputfile);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 		// get row number

				sscanf(input_string, "%*[^\t]%d", &rows);					// read all data up to tab as string ([^\t]), but do not store (*)		

			} else if (strtol(&input_string[0], NULL, 0) > 0) {					// else check the line begins with a positive number (usable)
	
				sscanf(input_string, "%d\t%lf\t%d\t", &fibre_number, &coord_x, &coord_y);
	
				array_coord_x_av [fibre_number-1] += coord_x;

				if ((rows == 1) || (coord_y < y_min)) {

					y_min = coord_y;

				} else if ((rows == 1) || (coord_y > y_min)) {

					y_max = coord_y;

				}
	
			}
	
		}

		int ii;

		for (ii=0; ii<FIBRES; ii++) {

			array_coord_x_av[ii] /= (double) rows;

			// printf("%d\t%f\n", ii+1, array_coord_x_av[ii]);	// DEBUG

		}

		// ***********************************************************************
		// Find rows where ANY peak x-coordinate lies more than [pixel_tolerance] 
		// from the average for that fibre and store in [bad_rows] array

		rewind(inputfile);

		int bad_rows [rows];
		memset(bad_rows, 0, sizeof(int)*rows);

		int bad_rows_number [rows];
		memset(bad_rows_number, 0, sizeof(int)*rows);

		int bad_row_count = 0;

		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*150);

			fgets(input_string, 150, inputfile);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 	// get row number

				sscanf(input_string, "%*[^\t]%d", &rows);	// read all data up to tab as string ([^\t]), but do not store (*)	

			} else if (strtol(&input_string[0], NULL, 0) > 0) {	// check the line begins with a positive number (usable)
	
				sscanf(input_string, "%d\t%lf\t%d\t", &fibre_number, &coord_x, &coord_y);

				if (fabs(coord_x-array_coord_x_av[fibre_number-1]) > pixel_tolerance)	{	// comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

					if (lsearch_int(bad_rows, coord_y, rows) == -1) {	// no entry for this [coord_y] already in array

						bad_rows[bad_row_count] = coord_y;		// so add it
						bad_rows_number[bad_row_count] = rows;		// and row number
						bad_row_count++;

					}

				}	
	
			}
	
		}	

		printf("\nRow cleaning results");
		printf("\n--------------------\n");
		printf("\nMinimum y coordinate:\t%d", y_min);
		printf("\nMaximum y coordinate:\t%d\n", y_max);

		printf("\nFound %d bad row(s):\n", bad_row_count);

		for (ii=0; ii<bad_row_count; ii++) {

			if (ii==0) {

				printf("\n");

			}

			printf("%d\n", bad_rows[ii]);

		}

		printf("\nUsable rows:\t\t%d\n", rows - bad_row_count);

		if (bad_row_count > max_bad_rows) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -3, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;


		} 

		// ***********************************************************************
		// Create [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] output file and print a few 
		// parameters

		FILE *outputfile;
		outputfile = fopen(FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -4, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(inputfile);

			return 1;


		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE);
		fprintf(outputfile, "# File cleaned using the frclean program.\n\n");
		fprintf(outputfile, "# Run datetime:\t%s\n", timestr);
		fprintf(outputfile, "\n# Min y coord:\t%d", y_min);
		fprintf(outputfile, "\n# Max y coord:\t%d\n", y_max);
		fprintf(outputfile, "\n# Bad rows:\t%d\n", bad_row_count);

		for (ii=0; ii<bad_row_count; ii++) {

			if (ii==0) {

				fprintf(outputfile, "\n");

			}

			fprintf(outputfile, "# %d\n", bad_rows[ii]);

		}

		fprintf(outputfile, "\n");
		fprintf(outputfile, "# Usable rows:\t%d\n\n", rows - bad_row_count);

		// ***********************************************************************
		// Rewrite [FRFIND_OUTPUTF_PEAKS_FILE]>[FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE]
		// erasing rows in the [bad_rows] array

		rewind(inputfile);

		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*150);

			fgets(input_string, 150, inputfile);	

			if (strtol(&input_string[0], NULL, 0) > 0) {	// check the line begins with a positive number (usable)
	
				sscanf(input_string, "%d\t%lf\t%d\t", &fibre_number, &coord_x, &coord_y);

				if (lsearch_int(bad_rows, coord_y, bad_row_count) == -1)	{		// not a bad row

					fprintf(outputfile, "%s", input_string);				// print to file as is else print nothing
					// if(fibre_number == 72) printf("%f\t%d\n", coord_x, coord_y);		// DEBUG

				} 
	
			} else {

				fprintf(outputfile, "%s", input_string);	// not a numeric entry, copy to file as is

			}

		}

		// ***********************************************************************
		// Close [FRFIND_OUTPUTF_PEAKS_FILE] input file and
		// [FRCLEAN_OUTPUTF_PEAKSCLEANED_FILE] output file

		if (fclose(inputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -5, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(outputfile);

			return 1; 

		}

		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -6, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1; 

		}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", RETURN_FLAG, "Status flag for L2 frclean routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}
	

}





