#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

void reset_cuda(){
    configdata   *cfg = &cconfig;

    for ( int j=0 ; j<cfg->ndev ; j++){
	int dev=cfg->devs[j];
	gpuErrchk(cudaSetDevice(dev));
	gpuErrchk(cudaDeviceReset());	
	gpuErrchk(cudaDeviceSynchronize());
    }

}

inline int _ConvertSMVer2Cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct
    {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
	{ 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
        { 0x60, 64 }, // Pascal Generation (SM 6.0)
	{ 0x61, 128}, // Pascal Generation (SM 6.1)
        {   -1, -1 }
    };

    int index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }

        index++;
    }

    // If we don't find the values, we default use the previous one to run properly
    printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[index-1].Cores);
    return nGpuArchCoresPerSM[index-1].Cores;
}


void checkmem(){
    size_t mem;
    size_t memreq;
    cudaDeviceProp	prop;
    configdata   *cfg = &cconfig;

    cfg->subp  = ceil(pxsize/subpix);					// Integer number of subpixels per pixel
    cfg->sizex = imsizex * ceil(pxsize/subpix);
    cfg->sizey = imsizey * ceil(pxsize/subpix);
    
    memreq       = sizeof(flux)*cfg->sizex*cfg->sizey*paths*IRlay+sizeof(input_h)*cfg->nx1*cfg->nx2;
    cfg->MPcount = (int *)malloc(sizeof(int)*cfg->ndev);
    cfg->cores   = (int *)malloc(sizeof(int)*cfg->ndev);

    for ( int i=0 ; i<cfg->ndev ; i++ ){
	cudaSetDevice(cfg->devs[i]);
	cudaDeviceReset();
	cudaDeviceSynchronize();
	gpuErrchk(cudaGetDeviceProperties(&prop,cfg->devs[i]));
	cfg->MPcount[i] = prop.multiProcessorCount;
	cfg->cores[i]   = _ConvertSMVer2Cores(prop.major,prop.minor);	// cores per MP
	size_t local    = cycles*cfg->cores[i]*cfg->MPcount[i]*sizeof(curandState_t);	// Additional mem for random number generator states per thread
	if ( memreq+local > prop.totalGlobalMem ){
	    printcudamem(prop.totalGlobalMem,cfg->devs[i],1);
	    printmemreqcuda(memreq+local,cfg->devs[i],1);
	    exit_with_usage(21);
	}
	printcudamem(prop.totalGlobalMem,cfg->devs[i],0);
	printmemreqcuda(memreq+local,cfg->devs[i],0);
    }
    mem    = getMemorySize();						// System flux
    if ( memreq >= mem ){
	printmem(mem,1);
	exit_with_usage(19);
    }
    printmem(mem,0);
    printmemreq(memreq,0);
}

void init_img(flux *rimg_h,flux **rimg_d){
    int i;
    configdata   *cfg = &cconfig;

    flux   img_h =  *rimg_h;
    flux  *img_d =  *rimg_d;

    img_h.phi_x = NULL;
    img_h.phi_y = NULL;
    img_h.N = NULL;

    int N = cfg->sizex*cfg->sizey*paths*IRlay;
    if (!img_h.phi_x)	img_h.phi_x = (float              *)malloc(sizeof(              float)*N);
    if (!img_h.phi_y)	img_h.phi_y = (float              *)malloc(sizeof(              float)*N);
    if (!img_h.N)	img_h.N     = (unsigned long long *)malloc(sizeof(unsigned long long )*N);


    for ( i=0 ; i < N ; i++ ){
        img_h.phi_x[i] = 0.0;
        img_h.phi_y[i] = 0.0;
        img_h.N[i]     = 0;
    }

    for ( i=0 ; i<cfg->ndev ; i++ ){
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
	img_d[i].phi_x = NULL;
	img_d[i].phi_y = NULL;
	img_d[i].N     = NULL;
	gpuErrchk(cudaMalloc(&img_d[i].phi_x,sizeof(float)*N));
    	gpuErrchk(cudaMemcpy(img_d[i].phi_x,img_h.phi_x,sizeof(float)*N,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMalloc(&img_d[i].phi_y,sizeof(float)*N));
	gpuErrchk(cudaMemcpy(img_d[i].phi_y,img_h.phi_y,sizeof(float)*N,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMalloc(&img_d[i].N,sizeof(unsigned long long)*N));
	gpuErrchk(cudaMemcpy(img_d[i].N,img_h.N,sizeof(unsigned long long)*N,cudaMemcpyHostToDevice));
    }

    *rimg_h = img_h;
    *rimg_d = img_d;
}

void init_psf(input_h *rpsf){

    input_h psf = *rpsf;
    
    psf.F  = NULL;
    psf.F0 = NULL;
    psf.x  = NULL;
    psf.y  = NULL;

    *rpsf = psf;
}

void fill_const(){
    configdata   *cfg = &cconfig;
    int          mult;
    float	 multd;

    for ( int i=0 ; i < cfg->ndev ; i++ ){
	gpuErrchk(cudaSetDevice(cfg->devs[i]));
        gpuErrchk(cudaMemcpyToSymbol(method,       &cfg->method,   sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(nx1,          &cfg->nx1,      sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(nx2,          &cfg->nx2,      sizeof(int),  0,cudaMemcpyHostToDevice));
	mult = cfg->nx1*cfg->nx2;
        gpuErrchk(cudaMemcpyToSymbol(nx1multnx2,   &mult,          sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(subp,         &cfg->subp,     sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(sizex,        &cfg->sizex,    sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(sizey,        &cfg->sizey,    sizeof(int),  0,cudaMemcpyHostToDevice));
	mult = cfg->sizex*cfg->sizey;
        gpuErrchk(cudaMemcpyToSymbol(szxmultszy,   &mult,          sizeof(int),  0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(oamp,         &cfg->oversamp, sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(rad,          &cfg->r,        sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(alpha,        &cfg->alpha,    sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(Kap,          &cfg->Kap,      sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(Trans,        &cfg->Trans,    sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(lam,          &cfg->lam,      sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(pxsc,         &cfg->pxsc,     sizeof(float),0,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpyToSymbol(Dap,          &cfg->D,        sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)(cfg->nx1)*0.5;
        gpuErrchk(cudaMemcpyToSymbol(nx1over2,     &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)(cfg->nx2)*0.5;
        gpuErrchk(cudaMemcpyToSymbol(nx2over2,     &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)pxsize/cfg->oversamp;
        gpuErrchk(cudaMemcpyToSymbol(pxoverosamp,  &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)pxsize/((float)cfg->subp);
        gpuErrchk(cudaMemcpyToSymbol(pxoversubp,   &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)(cfg->sizex)*0.5;
        gpuErrchk(cudaMemcpyToSymbol(szxover2,     &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = (float)(cfg->sizey)*0.5;
        gpuErrchk(cudaMemcpyToSymbol(szyover2,     &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = cfg->xoff-(Top-0.5*Dact)*sin(offtheta)*cos(offgamma);
        gpuErrchk(cudaMemcpyToSymbol(xoff,         &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
	multd = cfg->yoff-(Top-0.5*Dact)*sin(offtheta)*sin(offgamma);
        gpuErrchk(cudaMemcpyToSymbol(yoff,         &multd,         sizeof(float),0,cudaMemcpyHostToDevice));
    }
}
                                