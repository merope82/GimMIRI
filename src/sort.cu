#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

typedef struct{
 double        F,F0;           /* Input fits file fluxes                               */
 int           x,y;            /* Input fits file coordinates                          */
}input_aos;

int cmp(const void *a, const void *b){
    double F1 = *((double *)a);
    double F2 = *((double *)b);

    if ( F1 >  F2 )      return -1;
    else if ( F1 == F2 ) return  0;
    else                 return  1;
}

void sortinput(input_h *rpsf_h,input_d **rpsf_d){
    int    	    i,tot;
    double 	  T = 0.0;
    configdata *cfg = &cconfig;

    input_h  psf_h = *rpsf_h;
    input_d *psf_d = *rpsf_d;

    tot = cfg->nx1*cfg->nx2;
    for ( i=0 ; i<tot ; i++ ) T+=psf_h.F[i];
    cfg->Tot=T;

    T=0;
    for ( i=0 ; i<tot ; i++ ) {
	psf_h.F[i] = sqrt(psf_h.F[i]);
	T+=psf_h.F[i];
    }
    for ( i=0 ; i<tot ; i++ ) psf_h.F[i]/=T;


    // Not elegant, but easier to sort AOS with qsort than SOA,
    // and CUDA is much faster with SOA
    input_aos *psf = (input_aos *)malloc(sizeof(input_aos)*tot);
    for ( i=0 ; i<tot ; i++ ){
	psf[i].F  = psf_h.F[i];
	psf[i].F0 = psf_h.F0[i];
	psf[i].x  = psf_h.x[i];
	psf[i].y  = psf_h.y[i];
    }
    qsort(psf,tot,sizeof(input_aos),cmp);
    for ( i=0 ; i<tot ; i++ ){
	psf_h.F[i]  = psf[i].F;
	psf_h.F0[i] = psf[i].F0;
	psf_h.x[i]  = psf[i].x;
	psf_h.y[i]  = psf[i].y;
    }
    free(psf);

    T=0.0;
    for ( i=0 ; i<tot ; i++ ){		  // Cummulative percent up to pixel
	T+=psf_h.F[i];
	psf_h.F[i]=T;
    }

    for ( i=0 ; i<cfg->ndev ; i++ ){
        gpuErrchk(cudaSetDevice(cfg->devs[i]));
	psf_d[i].F = NULL;
	psf_d[i].x = NULL;
	psf_d[i].y = NULL;
    	gpuErrchk(cudaMalloc(&psf_d[i].F,sizeof(double)*tot));
    	gpuErrchk(cudaMemcpy(psf_d[i].F,psf_h.F,sizeof(double)*tot,cudaMemcpyHostToDevice));
    	gpuErrchk(cudaMalloc(&psf_d[i].x,sizeof(int)*tot));
    	gpuErrchk(cudaMemcpy(psf_d[i].x,psf_h.x,sizeof(int)*tot,cudaMemcpyHostToDevice));
    	gpuErrchk(cudaMalloc(&psf_d[i].y,sizeof(int)*tot));
    	gpuErrchk(cudaMemcpy(psf_d[i].y,psf_h.y,sizeof(int)*tot,cudaMemcpyHostToDevice));
    }

    *rpsf_h = psf_h;
    *rpsf_d = psf_d;
}
                             