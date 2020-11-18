#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imMIRI.h"

void print_header(const char * Version,const char *VDate){
 int i,length;
 configdata *cfg = &cconfig;

 if ( cfg->verb ==0 ) return;

 length=strlen(Version)+strlen(VDate);

 printf("\n");
 printf("#########################################################################\n");
 printf("#                           Running GimMIRI                             #\n");
 printf("#                        Ver %s, %s",Version,VDate);
 for ( i=0 ; i<41-length ; i++ ) printf(" ");printf("#\n");
 printf("#                                                                       #\n");
 printf("# Please contact Andras Gaspar if there seems to be a problem at:       #\n");
 printf("# agaspar@as.arizona.edu                                                #\n");
 printf("#                                                                       #\n");
 printf("#########################################################################\n");
 printf("\n");

}

void printstat(long l){
    // Print status of calculation - good to know
    configdata *cfg = &cconfig;

    if ( cfg->verb ==0 ) return;

    printf("Model at photon: %ld (%.2e)\r",l,(float)l);
    fflush(stdout);
}

void printdone(){
    configdata *cfg = &cconfig;
    if ( cfg->verb ==0 ) return;

    printf("Program done.\n");
}

void printmem(size_t mem,int force){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 && force==0 ) return;

    printf("Available random-excess memory in the system: %.2f Gb\n",mem/1e9);
}

void printcudamem(size_t mem,int device,int force){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 && force==0 ) return;

    printf("Available global CUDA memory in device number %d: %.2f Gb\n",device,mem/1e9);
}

void printmemreq(size_t mem,int force){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 && force==0 ) return;

    printf("Requested system memory: %.2f Gb\n",mem/1e9);
    if ( force==0 ) printf("There is sufficient memory for the modeling!\n");
}

void printmemreqcuda(size_t mem,int device,int force){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 && force==0 ) return;

    printf("Requested CUDA memory for device number %d: %.2f Gb\n",device,mem/1e9);
}

void printinitcuda(int device){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 ) return;
    
    printf("Initialized random number generator states for device %d\n",device);

}

void printinitheader(int n){
    configdata *cfg = &cconfig;
    if ( cfg->verb == 0 ) return;
    
    if (n==1) printf("\nInitializing random number generator states on CUDA device\n");
    else      printf("\nInitializing random number generator states on CUDA devices!\n");
              printf("This will take a few seconds.\n");

}
        