/*****************************************************************************/
/**                             GimMIRI                                     **/
/**                                                                         **/
/** Code to model the internal scattering of light in detectors, using      **/
/** probabilities derived from quantum electrodynamics.                     **/
/** GPU version								    **/
/*****************************************************************************/
/** Release version & date:                                                 **/
/**/ const char *Version = "1.40";                                         /**/
/**/ const char *VDate   = "11.18.2019";                                   /**/
/** Author: Andras Gaspar (agaspar@as.arizona.edu)                          **/
/** Revision: 0.01 - Initial						    **/
/** Revision: 1.00 - Released version					    **/
/**           Includes Fresnel Eqs for total internal reflection at top     **/
/** Revision: 1.10 - Major physics fixes				    **/
/**           - No reflection at buried contact				    **/
/**           - Picking photons with sqrt of I instead of linear ...OMG	    **/
/**           - Start modeling at grid, since no more reflection at B.C.    **/
/** Revision: 1.20 - Some changes					    **/
/**           - New absorption and transmission				    **/
/** Revision: 1.30 - Revert to reflection                                   **/
/** Revision: 1.40 - Revert to no-reflection. New absorption and trans.     **/
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <curand.h>
#include <curand_kernel.h>
#include "fitsio.h"
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

// Global config parameters
configdata cconfig;
// constant values stored in GPUs - check their meanings in the header file
__constant__ int   method;
__constant__ int   nx1;
__constant__ int   nx2;
__constant__ int   nx1multnx2;
__constant__ int   subp;
__constant__ int   sizex;
__constant__ int   sizey;
__constant__ int   szxmultszy;
__constant__ float oamp;
__constant__ float rad;
__constant__ float alpha;
__constant__ float Kap;
__constant__ float Trans;
__constant__ float lam;
__constant__ float pxsc;
__constant__ float Dap;
__constant__ float nx1over2;
__constant__ float nx2over2;
__constant__ float pxoverosamp;
__constant__ float pxoversubp;
__constant__ float szxover2;
__constant__ float szyover2;
__constant__ float xoff;
__constant__ float yoff;

int main(int argc,char *argv[]){
    configdata   *cfg   = &cconfig;
    char	 *param = NULL;
    flux	  img_h,*img_d;
    input_h	  psf_h;
    input_d	 *psf_d;
    globalstruct *global;
    
    // Define parameter file with input arguments
    if ( argc==1 )                            exit_with_usage(0);
    for ( int i=1 ; i<argc ; i++ ){
        if ( strcmp(argv[i],"-h")==0 )
            exit_with_usage(0);
        else if ( strcmp(argv[i],"-i")==0 ){
            i++;if ( i==argc )                exit_with_usage(1);
            param=argv[i];
        }
        else                                  exit_with_usage(0);
    }
    if ( param == NULL )		      exit_with_usage(2);

    read_in_param(param);
    print_header(Version,VDate);
    
    // Defining device pointers
    img_d  = (flux         *)malloc(sizeof(flux)        *cfg->ndev);
    psf_d  = (input_d      *)malloc(sizeof(input_d)     *cfg->ndev);		// Even if not used
    global = (globalstruct *)malloc(sizeof(globalstruct)*cfg->ndev);

    // Reset GPUs
    reset_cuda();

    // If using a fits input, get image data
    if ( cfg->method == 0 ) getimdata();

    // Check if there is enough system and GPU memory for model
    checkmem();

    // Initialize image on host and device(s) and the fits psf (if used)
    init_img(&img_h,&img_d);
    init_psf(&psf_h);

    // If using a fits input, read in file and copy to device(s)
    if ( cfg->method == 0 ){
	fillinput(&psf_h);
	sortinput(&psf_h,&psf_d);
    }

    // Copy constant values needed to constant memory in device(s)
    fill_const();

    // Run model evolution
    evolve(psf_d,psf_h,img_d,img_h,&global,Version);
    printdone();

    // Free system and device memory
    freeall(img_h,img_d,psf_h,psf_d,global);

 return 0;
}
                                                                            