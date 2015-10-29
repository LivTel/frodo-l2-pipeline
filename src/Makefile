#############################
# compilers, archivers etc. #
#############################

CC		= /usr/bin/gcc
AR		= /usr/bin/ar

##################
# compiler flags #
##################

CFLAGS		= -Wall -g -std=c99

#########################
# csh terminal commands #
#########################

CAT		= /bin/cat
GREP		= /bin/grep
MV		= /bin/mv
RM		= /bin/rm -i

################
# $(LIBS) flag #
################

LIBS		= -L$(L2_LIB_DIR)#						rmb/laptop/tower/lt-qc
#LIBS		= -L/space/home/dev/src/cfitsio3310 -L$(L2_LIB_DIR)#		lt/ltdevsrv

####################
# $(INCLUDES) flag #
####################

INCLUDES	= -I/usr/include/cfitsio -I/usr/include/gsl#			rmb/laptop/tower/lt-qc
#INCLUDES	= -I/space/home/dev/src/cfitsio3310 -I/usr/include/gsl#		lt/ltdevsrv

################
# make default #
################

all: frodo_red_findpeaks_simple frodo_red_findpeaks_simple_clean frodo_red_trace frodo_red_extract_simple frodo_red_arcfit frodo_red_correct_throughput frodo_red_rebin frodo_red_subsky frodo_red_reformat modhead imcopy liststruc frodo_aux_spec2tsv frodo_aux_numexts frodo_aux_peakfinder FITS2jpeg

#################
# make routines #
#################

frodo_red_findpeaks_simple: frodo_red_findpeaks_simple.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_findpeaks_simple.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frfind
	$(MV) frfind $(L2_BIN_DIR)

frodo_red_findpeaks_simple_clean: frodo_red_findpeaks_simple_clean.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_findpeaks_simple_clean.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frclean
	$(MV) frclean $(L2_BIN_DIR)

frodo_red_trace: frodo_red_trace.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_trace.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frtrace
	$(MV) frtrace $(L2_BIN_DIR)

frodo_red_extract_simple: frodo_red_extract_simple.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_extract_simple.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frextract
	$(MV) frextract $(L2_BIN_DIR)

frodo_red_arcfit: frodo_red_arcfit.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_arcfit.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frarcfit
	$(MV) frarcfit $(L2_BIN_DIR)

frodo_red_arcfit2: frodo_red_arcfit2.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_arcfit2.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frarcfit2
	$(MV) frarcfit2 $(L2_BIN_DIR)

frodo_red_correct_throughput: frodo_red_correct_throughput.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_correct_throughput.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frcorrect
	$(MV) frcorrect $(L2_BIN_DIR)		

frodo_red_rebin: frodo_red_rebin.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_rebin.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frrebin
	$(MV) frrebin $(L2_BIN_DIR)

frodo_red_subsky: frodo_red_subsky.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_subsky.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frsubsky
	$(MV) frsubsky $(L2_BIN_DIR)

frodo_red_reformat: frodo_red_reformat.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_red_reformat.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frreformat
	$(MV) frreformat $(L2_BIN_DIR)

modhead: modhead.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) modhead.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o modhead
	$(MV) modhead $(L2_BIN_DIR)	

imcopy: imcopy.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) imcopy.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o imcopy
	$(MV) imcopy $(L2_BIN_DIR)

liststruc: liststruc.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) liststruc.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o liststruc
	$(MV) liststruc $(L2_BIN_DIR)	

frodo_aux_spec2tsv: frodo_aux_spec2tsv.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_aux_spec2tsv.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frspec2tsv
	$(MV) frspec2tsv $(L2_BIN_DIR)	

frodo_aux_numexts: frodo_aux_numexts.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_aux_numexts.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frnumexts
	$(MV) frnumexts $(L2_BIN_DIR)	

FITS2jpeg: FITS2jpeg.c jpegsubs
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) FITS2jpeg.c jpegsubs.o -lcfitsio -ljpeg -lm -o FITS2jpeg 
	$(MV) FITS2jpeg $(L2_BIN_DIR)
	$(RM) jpegsubs.o

frodo_aux_peakfinder: frodo_aux_peakfinder.c libfrodo
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) frodo_aux_peakfinder.c -lfrodo -lgsl -lgslcblas -lcfitsio -lm -o frpeakfinder
	$(MV) frpeakfinder $(L2_BIN_DIR)

jpegsubs: jpegsubs.c
	$(CC) $(CFLAGS) -c jpegsubs.c

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

