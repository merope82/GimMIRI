#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

void readbackfromdevices(flux img_h,flux *img_d){
    configdata *cfg = &cconfig;
    flux        tmp;

    int N = cfg->sizex*cfg->sizey*paths*IRlay;
    tmp.phi_x = (float *)malloc(sizeof(float)*N);
    tmp.phi_y = (float *)malloc(sizeof(float)*N);
    tmp.N     = (unsigned long long *)malloc(sizeof(unsigned long long)*N);

    // Reset host image
    for ( int i=0 ; i<N ; i++ ){
	img_h.phi_x[i] = 0.0;
    	img_h.phi_y[i] = 0.0;
    	img_h.N[i]     = 0;
    }

    for ( int i=0 ; i<cfg->ndev ; i++ ){
	for ( int n=0 ; n<N ; n++ ){ 
	    tmp.phi_x[n] = 0.0;
	    tmp.phi_y[n] = 0.0;
	    tmp.N[n]     = 0;
	}
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
	gpuErrchk(cudaMemcpy(tmp.phi_x,img_d[i].phi_x,sizeof(float)*N,cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(tmp.phi_y,img_d[i].phi_y,sizeof(float)*N,cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(tmp.N,    img_d[i].N,    sizeof(unsigned long long)*N,cudaMemcpyDeviceToHost));
	gpuErrchk(cudaDeviceSynchronize());
	for ( int n=0 ; n<N ; n++ ){ 
	    img_h.phi_x[n] += tmp.phi_x[n];
	    img_h.phi_y[n] += tmp.phi_y[n];
	    img_h.N[n]     += tmp.N[n];
	    tmp.phi_x[n] = 0.0;
	    tmp.phi_y[n] = 0.0;
	    tmp.N[n]     = 0;
	}	
    }

    free(tmp.phi_x);
    free(tmp.phi_y);
    free(tmp.N);
}
                                                                           