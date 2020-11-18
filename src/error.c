#include <stdio.h>
#include <stdlib.h>
#include "imMIRI.h"

const char *errmsg[] = { "",
/* 1 */ "please give input parameter filename,",
/* 2 */ "missing input file,",
/* 3 */ "empty input file!",
/* 4 */ "please give yes or no to verbose parameter",
/* 5 */ "please give either fits, airy, or disk input value",
/* 6 */ "incorrect number of variables for fits model",
/* 7 */ "incorrect wavelength",
/* 8 */ "incorrect oversample value",
/* 9 */ "incorrect telescope aperture value",\
/* 10 */ "incorrect telescope focal length",\
/* 11 */ "incorrect number of variables for airy model",
/* 12 */ "incorrect disk size given",\
/* 13 */ "incorrect number of variables for disk model",
/* 14 */ "incorrect value for starting photon number value",
/* 15 */ "incorrect value for output multiplication factor",
/* 16 */ "incomplete parameter file",
/* 17 */ "image header does not have OVERSAMP keyword. Defaulting to 1.0",
/* 18 */ "image header does not have WAVE1 keyword. Defaulting to 5.6.",
/* 19 */ "not enough memory. Please decrease image size or increase\n\
          subpixel size in imMIRI.h",
/* 20 */ "CUDA error!",
/* 21 */ "not enough CUDA device memory. Please decrease image size\n\
	  or increase subpixel size in imMIRI.h",
/* 22 */ "incorrect last output photon value",
/* 23 */ "incorrect photon count variables",
/* 24 */ "please give cuda device id(s) in parameter file",
/* 25 */ "incorrect value(s) for x and/or y offset",
	NULL }; 

void print_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\t./imMIRI [-h] {help} [-i] <init.param> {initializing parameter file}\n");

 return;
}

void exit_with_usage(int t)
{
 if ( t>0  )         fprintf(stderr,"\nError:\t%s\n\n",errmsg[t]);
 if ( t<=2 )         print_usage(stderr);
 if ( t )            exit(1);
 else                exit(0);
}

void print_warning(int t){
    fprintf(stderr,"\nWarning:\t%s\n",errmsg[t]);
}
                        