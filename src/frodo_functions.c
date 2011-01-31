/************************************************************************

 File:				frodo_functions.c
 Last Modified Date:     	27/01/11

************************************************************************/

#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/stat.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>

/************************************************************************

 Function:		calc_least_sq_fit
 Last Modified Date:    25/01/11
 Purpose:		finds the coefficients of the [order] fit, for
			[equations], whose values are defined in [array_x]
			and [array_y]. stores in [coeff] where coeff[0] 
			is the 0th order of fit with a corresponding 
			[chi_squared].
 Required By:		frodo_red_findpeaks_simple.c
			frodo_red_arcfit.c
 Additional Notes:


 Performs least-squares fitting to a general model y = Xc using the GSL
 library.

 			Array formats
		 	-------------

 where m is the order and n is the number of equations.

 y = [	y_0, y_1, y_2, ..., y_n;	]

 X = [	1, (x_0)^1, (x_0)^2, ..., (x_0)^m;
	1, (x_1)^1, (x_1)^2, ..., (x_1)^m;
	1, (x_2)^1, (x_2)^2, ..., (x_2)^m;
	1, (x_3)^1, (x_3)^2, ..., (x_3)^m;
	.	.	.	.	.
	.	.	.	.	.
	.	.	.	.	.
	1, (x_n)^1, (x_n)^2, ..., (x_n)^m; ]

 c = [	c_0, c_1, c_2, ..., c_n+1;	]

 n.b. matrices are are denoted in (rows x columns) format.

************************************************************************/

int calc_least_sq_fit(int order, int equations, double array_x [], double array_y [], double coeffs [], double *chi_squared) {

	if (order >= equations)	{	// then we don't have enough equations to solve for all parameters

		return 1;

	}

	int ii, jj;

	gsl_matrix *X, *cov;
	gsl_vector *y, *c;

	double chisq;

	X = gsl_matrix_alloc(equations, order+1);	// matrix of predictor variables
	y = gsl_vector_alloc(equations);		// vector of observed values

	c = gsl_vector_alloc(order+1);			// coefficient vector
	cov = gsl_matrix_alloc(order+1, order+1);	// covariance matrix

	for (ii=0; ii<equations; ii++) {					// for each equation
			
		for (jj=0; jj<order+1; jj++) {
			   
			gsl_matrix_set(X, ii, jj, powf(array_x[ii], jj)); 	// populate X matrix with the value of x^m for each order m.

		}
			   
		gsl_vector_set(y, ii, array_y[ii]);				// and populate y matrix

	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(equations, order+1);	// define a workspace
	gsl_multifit_linear(X, y, c, cov, &chisq, work);					// perform the fit
	gsl_multifit_linear_free(work);								// free the workspace

	#define C(i) (gsl_vector_get(c,(i)))

	for (ii=0; ii<order+1; ii++) {

		coeffs[ii] = C(ii);	// populate coeff array
		
		// printf("%e\t", coeffs[ii]);	if (ii==order) { printf("\n"); } 	// DEBUG

	}

	*chi_squared = chisq;		// write chisq
	     
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	     
	return 0;

}

/************************************************************************

 Function:		calculate_cross_correlation_offset
 Last Modified Date:    17/01/11
 Purpose:		calculates the best offset for two given
			datasets, x and y.
 Required By:		frodo_red_arcfit.c
 Additional Notes:	

 The cross correlation coefficient, r, is determined by the equation:

 r =

		sum_i[(x(i) - mean_x) * y(i-offset) - mean_y)]
 ------------------------------------------------------------------------
 sqrt_(sum_i((x(i) - mean_x)^2)) * sqrt_(sum_i((y(i-offset) - mean_y)^2)) 

 a +ve offset indicates the series lags the reference series
 a -ve offset indicates the series leads the references series		

************************************************************************/

int calculate_cross_correlation_offset(double x [], double y [], int n, int max_offset, int *this_best_offset, double *this_best_r) {

	double mean_x = gsl_stats_mean(x, 1, n);
	double mean_y = gsl_stats_mean(y, 1, n);

	int offset, best_offset;
	double r, best_r;

	for (offset=-max_offset; offset<=max_offset; offset++) {
	
		double sx, sy, numerator = 0, den_x_term = 0, den_y_term = 0, denominator;

		int ii, jj;
	
		for (ii=0; ii<n; ii++) {
		
			jj = ii + offset;

			if (jj < 0 || jj >= n) {	// element is outside of array boundaries

				continue;

			} else {

				sx = (x[ii] - mean_x);
				sy = (y[jj] - mean_y);

				numerator += sx * sy;
				den_x_term += powf(sx, 2);
				den_y_term += powf(sy, 2);
			
			}

			denominator = sqrt(den_x_term) * sqrt(den_y_term);
	
		}

		r = numerator / denominator;

		if ((offset == -max_offset) || (r > best_r)) {

			best_r = r;
			best_offset = offset;

		} 

		// printf("%d\t%f\t%f\t%f\n", offset, numerator/denominator, numerator, denominator);	// DEBUG

	}

	*this_best_offset = best_offset;
	*this_best_r = best_r;

	// printf("%d\t%f\n", best_offset, best_r);		// DEBUG

	return 0;

}

/************************************************************************

 Function:		check_file_exists
 Last Modified Date:    25/01/11
 Purpose:		checks a file to see if it exists
 Required By:		frodo_red_findpeaks_simple_clean.c
			frodo_red_trace.c
			frodo_red_extract_simple.c
			frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int check_file_exists(char filename []) {

	struct stat frfind_outputf_peaks_file;
	stat(filename, &frfind_outputf_peaks_file);
	
	if (S_ISREG(frfind_outputf_peaks_file.st_mode) == 1) {
		
		return 0;
	
	} else {
	
		return 1;
	
	}

}

/************************************************************************

 Function:		find_centroid_parabolic	
 Last Modified Date:    25/01/11
 Purpose:		finds the parabolic centroid of each peak in a
			dataset given the position of the peaks
 Required By:		frodo_red_findpeaks_simple.c
			frodo_red_arcfit.c
 Additional Notes:	None


 [INDEXING CORRECTION]:

 This variable defines how if we offset the arrays.

 If the correction is set to TRUE, [peaks_centroids] will store the x 
 value of the identified peak and not the ii value (x-1).

************************************************************************/

int find_centroid_parabolic(double row_values [], int peaks [], int num_peaks, double peak_centroids [], bool INDEXING_CORRECTION) {
	     
	int ii, jj;

	int order = 2;		// this is a quadratic procedure
	int equations = 3;	// using 3 pixels (i.e. fully defined)

	double array_x [equations];
	memset(array_x, 0, sizeof(double)*equations);

	double array_y [equations];
	memset(array_y, 0, sizeof(double)*equations);

	double coeffs [order+1];  
	memset(coeffs, 0, sizeof(double)*order+1);
 
	double chi_squared;

	for (ii=0; ii<num_peaks; ii++) {		// for each peak

		for (jj=0; jj<=order; jj++) {		// populate array_x and array_y with the appropriate variables
			
			array_x[jj] = peaks[ii] + (jj-1);			
			array_y[jj] = row_values[peaks[ii]+(jj-1)-INDEXING_CORRECTION];	// see INDEXING_CORRECTION

		}

		if (calc_least_sq_fit(order, equations, array_x, array_y, coeffs, &chi_squared)) {

			return 1;

		} else {

			peak_centroids[ii] = -coeffs[1]/(2*coeffs[2]);	// store centroid (derivative of y = ax^2 + bx + c is -b/2a)

		}

	}

	return 0;

}

/************************************************************************

 Function:		find_peaks
 Last Modified Date:    17/01/11	
 Purpose:		finds the peaks in a dataset
 Required By:		frodo_red_findpeaks_simple.c
			frodo_red_arcfit.c > frodo_functions.c 	
						(find_peaks_contiguous)
 Additional Notes:	


 [INDEXING CORRECTION]:

 This variable defines how if we offset the arrays.

 If the correction is set to TRUE, [peaks] will store the x value of the
 identified peak and not the ii value (x-1).

************************************************************************/

int find_peaks(int nxelements, double row_values [], int peaks [], int * num_peaks, int min_dist, int half_aperture_num_pix, int derivative_tol, int derivative_tol_ref_px, bool INDEXING_CORRECTION) {

	int ii, jj;
	int peak_count = 0;
	
	for (ii = half_aperture_num_pix; ii < nxelements - half_aperture_num_pix; ii++) {  		// this is the pixel we're considering (ii)
	
		for (jj = ii - half_aperture_num_pix; jj <= ii + half_aperture_num_pix; jj++) {		// this is the aperture we're considering

			if (jj == ii) { 

				continue; // this is the pixel in question, skip

			} else if (row_values[jj] > row_values[ii]) { 

				break;	  // this pixel is not a peak within the specified aperture

			} else if ((jj == ii + half_aperture_num_pix) && ((row_values[ii] - row_values[ii-derivative_tol_ref_px]) > derivative_tol && (row_values[ii] - row_values[ii+derivative_tol_ref_px]) > derivative_tol)) {	// reached last pixel in aperture => pixel ii is a peak within it. Check that it is sufficiently brighter than the derivative tolerance value / pixel

				if (peak_count == 0) {	// if this is the first peak, we can just add it without checking the minimum distance parameter (no other peaks exist)

					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;

				} else if (ii - peaks[peak_count-1] > min_dist) { // else we need to check the minimum distance criteria is satisfied
	
					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;
	
				}
	
			}
	
		}
	
	}
		
	*num_peaks = peak_count;

	return 0;

}

/************************************************************************

 Function:		find_peaks_contiguous
 Last Modified Date:    26/01/11	
 Purpose:		finds contiguous peaks in a dataset
 Required By:		frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int find_peaks_contiguous(int nxelements, int nyelements, double frame_values [nyelements][nxelements], int peaks [nyelements][nxelements], int * num_peaks, int min_dist, int half_aperture_num_pix, int derivative_tol, int derivative_tol_ref_px, int pix_tolerance, bool INDEXING_CORRECTION) {

	int ii, jj;

	double row_values [nxelements];

	int this_row_peaks [nxelements];					// array to store peak locations for a single iteration

	int this_row_num_peaks;							// variable to store number of peaks for a single iteration

	int all_rows_peaks [nyelements][nxelements];				// array to store peak locations for all iterations
	memset(all_rows_peaks, 0, sizeof(int)*nyelements*nxelements);	

	int all_rows_num_peaks [nyelements];					// array to store number of peaks for all iterations
	memset(all_rows_num_peaks, 0, sizeof(int)*nyelements*nxelements);

	for (jj=0; jj<nyelements; jj++) {

		memset(row_values, 0, sizeof(double)*nxelements);
		memset(this_row_peaks, 0, sizeof(int)*nxelements);	

		for (ii=0; ii<nxelements; ii++) {

			row_values[ii] = frame_values[jj][ii];

		}

		find_peaks(nxelements, row_values, this_row_peaks, &this_row_num_peaks, min_dist, half_aperture_num_pix, derivative_tol, derivative_tol_ref_px, INDEXING_CORRECTION);

		for (ii=0; ii<this_row_num_peaks; ii++) {

			all_rows_peaks[jj][ii] = this_row_peaks[ii];		// copy peaks locations across

		}

		all_rows_num_peaks[jj] = this_row_num_peaks;			// copy number of peaks across

	}

	int kk, index_to_insert = 0;

	bool insert;

	int considered_peaks [nyelements];	// array to hold locations of the considered peak

	double this_pix_difference = 0;
	double least_pix_difference;

	for (ii=0; ii<all_rows_num_peaks[0]; ii++) {									// for each peak on the first row

		memset(considered_peaks, 0, sizeof(int)*nyelements);

		considered_peaks[0] = all_rows_peaks[0][ii];

		for (jj=1; jj<nyelements; jj++) {									// cycle each futher row

			insert = FALSE;

			for (kk=0; kk<all_rows_num_peaks[jj]; kk++) {							// and each peak on the row

				this_pix_difference = abs(all_rows_peaks[jj][kk] - all_rows_peaks[0][ii]);		// calculate the pixel distance between this peak and the peak from the first row

				if (this_pix_difference <= pix_tolerance) {						// found a suitable candidate peak within the right tolerance

					if ((insert == FALSE) || (this_pix_difference < least_pix_difference)) {	// if it's the first cycle OR there's a closer peak

						insert = TRUE;
						considered_peaks[jj] = all_rows_peaks[jj][kk];
						least_pix_difference = this_pix_difference;				// this is always set on each iteration of jj as insert is set to FALSE

					}
					
				} 

			}

			if (insert == FALSE) {

				break;

			}

		}

		if (insert == TRUE) {	// considered_peaks can be inserted

			for (jj=0; jj<nyelements; jj++) {

				peaks[jj][index_to_insert] = considered_peaks[jj];

			}

			*num_peaks = *num_peaks + 1;
			index_to_insert++;

		}

	}

	return 0;

}

/************************************************************************

 Function:		find_time
 Last Modified Date:    25/01/11
 Purpose:		finds the time
 Required By:		frodo_red_findpeaks_simple.c
			frodo_red_findpeaks_simple_clean.c
			frodo_red_trace.c
			frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int find_time (char timestr []) {

	struct tm *time_now;
	time_t ltime;

	ltime = time(NULL);

	time_now = localtime(&ltime);

	strftime(timestr, 80, "%d-%m-%Y %H:%M:%S", time_now);

	return 0;

}

/************************************************************************

 Function:		flip_array_dbl
 Last Modified Date:    19/01/11	
 Purpose:		Flips an array
 Required By:		frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int flip_array_dbl(double array [], int size) {

	int ii;

	double flip_array [size];
	memset(flip_array, 0, sizeof(double)*size);

	for (ii=0; ii<size; ii++) {

		flip_array[size-1-ii] = array[ii];	// store reversed order of array elements to temporary array [flip_array]
	
	}

	memset(array, 0, sizeof(double)*size);

	for (ii=0; ii<size; ii++) {

		array[ii] = flip_array[ii];		// rewrite array [array] with the reversed array [flip_array]
	
	}

	return 0;

}

/************************************************************************

 Function:		lsearch_int
 Last Modified Date:    25/01/11	
 Purpose:		Performs a search for value [key] in int [array]
			of size [size].
 Required By:		frodo_red_findpeaks_simple_clean.c
			frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int lsearch_int(int array [], int key, int size) {

	int n;

	for (n=0; n<size; n++) {

		if (array[n] == key) { 

			return n;	// returns array index where value was found

		} 

	} 

	return -1;			// otherwise exited normally but didn't find value

}

/************************************************************************

 Function:		populate_env_variable
 Last Modified Date:    27/01/11	
 Purpose:		populates global variables with corresponding
			environment variables from the shell
 Required By:		all
 Additional Notes:	None

************************************************************************/

int populate_env_variable(char var_to_populate [], char env_var_name []) {

	if (!getenv(env_var_name)) {

		return 1;

	} else {

		strcpy(var_to_populate, getenv(env_var_name));

		return 0;

	}

}

/************************************************************************

 Function:		populate_img_parameters
 Last Modified Date:    11/01/11	
 Purpose:		populates variables with corresponding parameters
			and prints to stdout
 Required By:		frodo_red_findpeaks_simple.c
			frodo_red_extract_simple.c
			frodo_red_arcfit.c
 Additional Notes:	None

************************************************************************/

int populate_img_parameters(char f [], fitsfile *f_ptr, int maxdim, int *bitpix, int *naxis, long naxes [], int *status, char title_of_img []) {

	if(!fits_get_img_param(f_ptr, maxdim, bitpix, naxis, naxes, status)) {

		printf("\nFile parameters");
		printf("\n---------------\n");
		printf("\n");
		printf("Frame name:\t\t%s\n", title_of_img);
		printf("Relative path:\t\t%s\n", f);
		printf("Bits per pixel:\t\t%d\n", *bitpix);
		printf("Number of axes:\t\t%d\n", *naxis);
		printf("First axis dimension:\t%ld\n", naxes[0]);
		printf("Second axis dimension:\t%ld\n", naxes[1]);

		return 0;
		
	} else {

		return 1;

	}

}

/************************************************************************

 Function:		print_file
 Last Modified Date:    26/01/11
 Purpose:		prints content of file to screen
 Required By:		all
 Additional Notes:	None

************************************************************************/

int print_file(char TEXT_FILE [200]) {

	FILE *TEXT_FILE_fptr;
	TEXT_FILE_fptr = fopen(TEXT_FILE, "r");

	char input_string [300];

	if (TEXT_FILE_fptr) {
	
		while(!feof(TEXT_FILE_fptr)) {

			memset(input_string, '\0', 300);

			fgets (input_string, 300, TEXT_FILE_fptr);
	
			printf("%s", input_string);
	
		}

	} else {

		printf("\nWARNING:\tUnable to print file. File %s doesn't exist.\n", TEXT_FILE);
		return 1;

	}

	return 0;

}

