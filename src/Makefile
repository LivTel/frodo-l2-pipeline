##################################
# Makefile 			 #
# Last Modified Date:   08/05/11 #
##################################

#############################
# compilers, archivers etc. #
#############################

CC		= /usr/bin/gcc
AR		= /usr/bin/ar

##################
# compiler flags #
##################

CFLAGS		= -Wall -g

#########################
# csh terminal commands #
#########################

CAT		= /bin/cat
GREP		= /bin/grep
MV		= /bin/mv
RM		= /bin/rm

################
# $(LIBS) flag #
################

#LIBS 		= /space/home/eng/frodo_pipeline/cfitsio3.13/lib/libcfitsio.a /usr/lib/libm.a# 					RMB/DEPLOY/LTDEVSRV
LIBS		= /usr/lib/libcfitsio.so /usr/lib/libm.so /usr/lib/libgsl.so /usr/lib/libgslcblas.so $(L2_LIB_DIR)/libfrodo.a#	RMB/DEVEL/LAPTOP

####################
# $(INCLUDES) flag #
####################

#INCLUDES	= -I/space/home/eng/frodo_pipeline/cfitsio3.13/include#								RMB/DEPLOY/LTDEVSRV
INCLUDES	= -I/usr/include/cfitsio -I/usr/include/gsl#									RMB/DEVEL/LAPTOP

################
# make default #
################

all: frodo_red_findpeaks_simple frodo_red_findpeaks_simple_clean frodo_red_trace frodo_red_extract_simple frodo_red_arcfit frodo_red_correct_throughput frodo_red_rebin frodo_red_subsky frodo_red_reformat modhead imcopy liststruc frodo_aux_spec2tsv frodo_aux_numexts

#################
# make routines #
#################

frodo_red_findpeaks_simple: frodo_red_findpeaks_simple.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frfind frodo_red_findpeaks_simple.c $(LIBS)
	$(MV) frfind $(L2_BIN_DIR)

frodo_red_findpeaks_simple_clean: frodo_red_findpeaks_simple_clean.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frclean frodo_red_findpeaks_simple_clean.c $(LIBS)
	$(MV) frclean $(L2_BIN_DIR)

frodo_red_trace: frodo_red_trace.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frtrace frodo_red_trace.c $(LIBS)
	$(MV) frtrace $(L2_BIN_DIR)

frodo_red_extract_simple: frodo_red_extract_simple.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frextract frodo_red_extract_simple.c $(LIBS)
	$(MV) frextract $(L2_BIN_DIR)

frodo_red_arcfit: frodo_red_arcfit.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frarcfit frodo_red_arcfit.c $(LIBS)
	$(MV) frarcfit $(L2_BIN_DIR)

frodo_red_correct_throughput: frodo_red_correct_throughput.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frcorrect frodo_red_correct_throughput.c $(LIBS)
	$(MV) frcorrect $(L2_BIN_DIR)		

frodo_red_rebin: frodo_red_rebin.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frrebin frodo_red_rebin.c $(LIBS)
	$(MV) frrebin $(L2_BIN_DIR)

frodo_red_subsky: frodo_red_subsky.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frsubsky frodo_red_subsky.c $(LIBS)
	$(MV) frsubsky $(L2_BIN_DIR)

frodo_red_reformat: frodo_red_reformat.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frreformat frodo_red_reformat.c $(LIBS)
	$(MV) frreformat $(L2_BIN_DIR)

modhead: modhead.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o modhead modhead.c $(LIBS)
	$(MV) modhead $(L2_BIN_DIR)	

imcopy: imcopy.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o imcopy imcopy.c $(LIBS)
	$(MV) imcopy $(L2_BIN_DIR)

liststruc: liststruc.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o liststruc liststruc.c $(LIBS)
	$(MV) liststruc $(L2_BIN_DIR)	

frodo_aux_spec2tsv: frodo_aux_spec2tsv.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frspec2tsv frodo_aux_spec2tsv.c $(LIBS)
	$(MV) frspec2tsv $(L2_BIN_DIR)	

frodo_aux_numexts: frodo_aux_numexts.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) -o frnumexts frodo_aux_numexts.c $(LIBS)
	$(MV) frnumexts $(L2_BIN_DIR)	

frodo_functions: frodo_functions.c
	$(CC) $(CFLAGS) $(INCLUDES) -c frodo_functions.c

frodo_error_handling: frodo_error_handling.c
	$(CC) $(CFLAGS) $(INCLUDES) -c frodo_error_handling.c

############
# libfrodo #
############

libfrodo: frodo_functions frodo_error_handling
	$(AR) rcs libfrodo.a frodo_functions.o frodo_error_handling.o
	$(MV) libfrodo.a $(L2_LIB_DIR)
	$(RM) *.o

##############
# make clean #
##############

clean:
	$(RM) $(L2_LIB_DIR)/*
	$(RM) $(L2_BIN_DIR)/*

