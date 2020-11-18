#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"


void freeall(flux img_h,flux *img_d,input_h psf_h,input_d *psf_d,globalstruct *global){
    configdata *cfg = &cconfig;

    free(img_h.phi_x);
    free(img_h.phi_y);
    free(img_h.N);
    if (!psf_h.F)  free(psf_h.F);
    if (!psf_h.F0) free(psf_h.F0);
    if (!psf_h.x)  free(psf_h.x);
    if (!psf_h.y)  free(psf_h.y);

    for ( int i=0 ; i<cfg->ndev ; i++){
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
 	gpuErrchk(cudaFree(img_d[i].phi_x));
	gpuErrchk(cudaFree(img_d[i].phi_y));
	gpuErrchk(cudaFree(img_d[i].N));
	gpuErrchk(cudaFree(global[i].states));
	if (!psf_d[i].F) gpuErrchk(cudaFree(psf_d[i].F));
        if (!psf_d[i].x) gpuErrchk(cudaFree(psf_d[i].x));
        if (!psf_d[i].y) gpuErrchk(cudaFree(psf_d[i].y));
    }
}