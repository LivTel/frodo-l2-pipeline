/************************************************************************

 File:				frodo_config.h
 Last Modified Date:     	08/03/11

************************************************************************/

char HEADER_FILE [100];
int SUCCESS_FLAG				= 0;

bool INDEXING_CORRECTION			= TRUE;

int FIBRES 					= 144;

char FILE_WRITE_ACCESS [2] 			= "w";				// r readonly; w overwrite; a+ append

int IMG_READ_ACCURACY				= 82;
int INTERMEDIATE_IMG_ACCURACY[2]		= {-64, 82};			// Definition follows types from http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node20.html {create_img, write_img}





