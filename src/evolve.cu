#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

void evolve(input_d *psf_d,input_h psf_h,flux *img_d,flux img_h,globalstruct **rglobal,const char *Version){
    configdata	*cfg	= &cconfig;
    long	j	= 0;
    long	Nwrite	= cfg->Nstart;

    int		 tpb	 = 128;
    int		 *chunk	 = (int *)malloc(sizeof(int)*cfg->ndev);
    int		 *blocks = (int *)malloc(sizeof(int)*cfg->ndev);
    globalstruct *global = *rglobal;

// Set thread-per-block and block count numbers per device
// Scaled to core count and by cycles variable
    for ( int i=0 ; i<cfg->ndev ; i++ ){
	chunk[i]	= cycles*cfg->cores[i]*cfg->MPcount[i];
	blocks[i]	= (chunk[i]-1.0)/tpb+1.0;
    }

// Init curand_states
// By initiating once per device and storing values, the code is much faster
    printinitheader(cfg->ndev);
    for ( int i=0 ; i<cfg->ndev ; i++ ){
        global[i].states = NULL;
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
        gpuErrchk(cudaMalloc(&global[i].states,sizeof(curandState_t)*chunk[i]));
	int offset = 0;
	int d=0;
	while ( d<i ){ offset+=chunk[d]; d++; }
        init_curand<<<blocks[i],tpb>>>(chunk[i],global[i].states,2.0,offset);
    }
    for ( int i=0 ; i<cfg->ndev ; i++ ){
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
	cudaError err = cudaGetLastError();
	if ( cudaSuccess != err ){
		print_warning(20);
		fprintf(stderr,"CUDA kernel execution error on device %d: %s\n",i,cudaGetErrorString(err));
		exit(-1);
	}
	gpuErrchk(cudaDeviceSynchronize());
	printinitcuda(cfg->devs[i]);
    }

    do{
	do{
	    for ( int i=0 ; i<cfg->ndev ; i++ ){
		gpuErrchk(cudaSetDevice(cfg->devs[i]));
        	evolve_photons<<<blocks[i],tpb>>>(chunk[i],global[i].states,img_d[i].phi_x,img_d[i].phi_y,img_d[i].N,
		    psf_d[i].x,psf_d[i].y,psf_d[i].F);
		j += (long)chunk[i];
    	    }
	    for ( int i=0 ; i<cfg->ndev ; i++ ){
		gpuErrchk(cudaSetDevice(cfg->devs[i]));
		cudaError err = cudaGetLastError();
		if ( cudaSuccess != err ){
		    print_warning(20);
		    fprintf(stderr,"CUDA kernel execution error on device %d: %s\n",i,cudaGetErrorString(err));
		    exit(-1);
		}
		gpuErrchk(cudaDeviceSynchronize());
	    }
	    printstat(j);
	}while(j<Nwrite);

    Nwrite = (long)(cfg->Nmult*Nwrite);
    readbackfromdevices(img_h,img_d);
    write_fits(img_h,psf_h,j,Version);
    }while ( j <= cfg->Nlast );

    *rglobal = global;
}
                                                           