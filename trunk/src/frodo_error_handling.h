/************************************************************************

 File:				frodo_error_handling.h
 Last Modified Date:     	07/03/11

************************************************************************/

int get_error_description(FILE, int, char *, char *, char *);
int write_key_to_file(char *, char *, char *, int, char *, char *);

char ERROR_CODES_INITIAL_FILE_WRITE_ACCESS [2]	= "w";			// r readonly; w overwrite; a+ append
char ERROR_CODES_FILE_WRITE_ACCESS [2]	 	= "a+";			// r readonly; w overwrite; a+ append

char REF_ERROR_CODES_FILE [100];					// variable to hold location of error codes reference file
char ERROR_CODES_FILE [100] 			= "error_codes";	// location of file to write error codes to
char ADDITIONAL_KEYS_FILE [100] 		= "additional_keys";	// location of file to write additional keys to

