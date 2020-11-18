#ifndef __GIMMIRI_CUH__
#define __GIMMIRI_CUH__

/* Random states */
typedef struct{
    curandState_t *states;
}globalstruct;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert error: %s\nin code: %s\nline: %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ float getairyR(float);
__global__ void  evolve_photons(int,curandState_t *,float *,float *,unsigned long long *,int *,\
                               int *,double *);
__global__ void  init_curand(int,curandState_t *,unsigned int,int);

void   evolve(input_d *,input_h,flux *,flux,globalstruct **,const char *);
void   freeall(flux,flux *,input_h,input_d *,globalstruct *);

extern __constant__ int   method;		/* Describes flux input method	*/
extern __constant__ int   nx1;			/* PSF image size used		*/
extern __constant__ int   nx2;			/* with oversampling		*/
extern __constant__ int   nx1multnx2;
extern __constant__ int   subp;			/* Sub-pixel actual size	*/
extern __constant__ int   sizex;		/* Output image size		*/
extern __constant__ int   sizey;
extern __constant__ int   szxmultszy;
extern __constant__ float oamp;			/* Oversampling			*/
extern __constant__ float rad;			/* Radius of pinhole		*/
extern __constant__ float alpha;		/* Absorption coeff		*/
extern __constant__ float Kap;			/* 2pi/lambda			*/
extern __constant__ float Trans;		/* Transmission fraction	*/
extern __constant__ float lam;			/* Wavelength in detector	*/
extern __constant__ float pxsc;			/* Pixelscale			*/
extern __constant__ float Dap;			/* Aperture of Telescope	*/
extern __constant__ float nx1over2;		/* Half image size		*/
extern __constant__ float nx2over2;
extern __constant__ float pxoverosamp;		
extern __constant__ float pxoversubp;
extern __constant__ float szxover2;
extern __constant__ float szyover2;
extern __constant__ float xoff;
extern __constant__ float yoff;

#endif
                                                           