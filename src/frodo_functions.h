/************************************************************************

 File:				frodo_functions.h
 Last Modified Date:     	08/03/11

************************************************************************/

#include <stdbool.h>

int calc_least_sq_fit (int, int, double [], double [], double [], double *);
int calculate_cross_correlation_offset (double [], double [], int, int, int *, double *);
int check_file_exists (char []);
int find_centroid_parabolic (double [], int [], int, double [], bool);
int find_peaks(int, double [], int [], int *, int, int, int, int, bool);
int find_peaks_contiguous(int, int, double **, int **, int *, int, int, int, int, int, bool);
int find_time (char []);
int flip_array_dbl(double [], int);
int interpolate(char [], double [], double [], int, double, double, double, double []);
int iterative_sigma_clip(double [], int, int, int [], double *, double *, int *);
int lsearch_int(int [], int, int);
int populate_img_parameters(char [], fitsfile *, int, int *, int *, long [], int *, char []);
int populate_env_variable(char [], char []);
int print_file(char *);



