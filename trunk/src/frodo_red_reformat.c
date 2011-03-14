/************************************************************************

 File:				frodo_red_reformat.c
 Last Modified Date:     	14/03/11

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
#include "frodo_red_reformat.h"

// *********************************************************************/

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 2) {

		if(populate_env_variable(HEADER_FILE, "L2_HEADER_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", 1, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		if(populate_env_variable(FRRF_BLURB_FILE, "L2_FRRF_BLURB_FILE")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", 2, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		}

		print_file(HEADER_FILE);
		print_file(FRRF_BLURB_FILE);

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -1, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

	}

}
