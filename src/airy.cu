#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>

#include "imMIRI.h"
#include "GimMIRI.cuh"

__device__ float getairyR(float F){
  float	      r,Tot,rmax,rmin,Phalf,c;

 c = M_PI * Dap * 1.0e6 * (float)pxsc * M_PI / (180.0*3600.0) / ( lam * nsil );
 if ( imsizex < imsizey ) rmax=((float)imsizex-2.0)*0.5;
 else			  rmax=((float)imsizey-2.0)*0.5;
 float                    up = c * rmax;

 Tot = 1.0-j0(up)*j0(up)-j1(up)*j1(up);

 r = rmax*0.5;
 up = c * r;
 Phalf = (1.0 - j0(up)*j0(up) - j1(up)*j1(up)) / Tot;

 rmin=0;
 do{
    if ( Phalf > F ) rmax=(rmax+rmin)/2.0;
    else	     rmin=(rmax+rmin)/2.0;
    r = (rmax+rmin)/2.0;
    up = c * r;
    Phalf = (1.0 - j0(up)*j0(up) - j1(up)*j1(up)) / Tot;
 }while ( powf(rmax - rmin,2.0) > 1.0e-8 );

 return(r*(float)pxsize);
}
                                                                   