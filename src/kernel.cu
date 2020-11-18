#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <curand.h>
#include <curand_kernel.h>
#include "fitsio.h"
#include "imMIRI.h"
#include "GimMIRI.cuh"
#include "externC.h"

__global__ void init_curand(int subj,curandState_t *global,unsigned int seed,int offset){

    int t = threadIdx.x + blockIdx.x * blockDim.x;
    if ( t < subj ){
	curand_init(seed,t,(offset+t)*1e20,&global[t]);
    }
}

__global__ void evolve_photons(int subj,curandState_t *global,float *phi_x,float *phi_y,\
                               unsigned long long *N,int *psf_x,int *psf_y,double *psf_F){

        int t = threadIdx.x + blockIdx.x * blockDim.x;

        if ( t < subj){
	    float xs,ys    = 0.0;
	    int   next     = 1;								// Start at reflection already!
	    int   pass     = 0;

	    float costheta = cosf(offtheta);
	    float sintheta = sinf(offtheta);
	    float cosgamma = cosf(offgamma);
	    float singamma = sinf(offgamma);
	    float r,z,temp = 1.0;
	    float r0       = 0.0;

            curandState_t state = global[t];

	    if ( method == 0 ){								// Fits file
		// binary search algorithm for PSF pixel
		int    is=-1,il=nx1multnx2-1,im=(int)((nx1multnx2-2)*0.5);
		double Fl=psf_F[il];
		double Fm=psf_F[im];
		double F;
		do { F=curand_uniform(&state); } while(F>=Fl);				// Just in case Fl 
											// is not rounded well at 1
		do{
		    if ( F <= Fm ){
			il=im;
			im=(int)((il+is)*0.5);
		    }
		    else{
			is=im;
			im=(int)((il+is)*0.5);
		    }
		Fm=psf_F[im];
		}while(il-is>1);
		
		xs = (curand_uniform(&state)+(float)(psf_x[il]-nx1over2))*pxoverosamp;	// Random location within found pixel
		temp *= curand_uniform(&state);						// Extra call. States seem to be correlated
	        ys = (curand_uniform(&state)+(float)(psf_y[il]-nx2over2))*pxoverosamp;
	    }
	    if ( method == 1 ){								// Airy profile  ... this is the slowest
		float F = sqrt(curand_uniform(&state));
	        float R = getairyR(F); 							// Returns Airy radial location from cummulative flux distribution
		temp *= curand_uniform(&state);						// Extra call. States seem to be correlated
	        float P = twopi*curand_uniform(&state);
	        xs = R * cosf(P);
	        ys = R * sinf(P);
	    }
	    if ( method == 2 ){								// Pinhole method
		float R = rad*sqrtf(curand_uniform(&state));
		temp *= curand_uniform(&state);						// Extra call. States seem to be correlated
	        float P = twopi*curand_uniform(&state);
	        xs = R * cosf(P);
	        ys = R * sinf(P);
	    }
	    xs+=xoff;
	    ys+=yoff;
    
	    do{
	    	if ( next==0 ){ // Evolve down in absorption layer

		pass++;

	        // Calculate distance for random absorption probability
		// curand_uniform does give 1.0, so inverting here.
	        r = -log(1.0-curand_uniform(&state))/alpha;
	        // That distance corresponds to the following z distance
	        z = Dact-r*costheta;
	        // if z < 0, light gets mirrored. Recalc r to contact surface.
		int l;
	        if ( z <= 0.0 ){
	    	    z = 0.0;
		    r = Dact/costheta;
	        }
	        else{ // Find layer center
		    l = floor(z*IRlay/Dact);
		    if ( l==IRlay ) l--;
	    	    z = (l+0.5)*(Dact/IRlay);
		    r = (l+0.5)*(Dact/IRlay)/costheta;
	        }
	        // Calculate new x and y positions at layer center.
		float xo = xs;
		float yo = ys;
	        xs += r*sintheta*cosgamma;
	        ys += r*sintheta*singamma;

		// Recentering will help with calculating the diffraction correctly.
		if ( z > 0 ){
		    xs = (floorf(xs/pxoversubp)+0.5)*pxoversubp; 
		    ys = (floorf(ys/pxoversubp)+0.5)*pxoversubp; 
		    r  = sqrtf(powf(xo-xs,2.0)+powf(yo-ys,2.0)+powf(Dact-z,2.0));
		}
	        // Find region in image. If off image, exit. If not, and absorbed, recenter
	        // to pixel center.

	        if ( xs < -halfx || xs >= halfx || ys < -halfy || ys >= halfy ){
		    global[t]=state;
		    return;
	        }

		r0 += r;

		if ( z==0.0 ) next  = 1;
		else{
		    // Find img id.
		    xs = floorf(xs/pxoversubp) + szxover2; 
		    ys = floorf(ys/pxoversubp) + szyover2; 
		    if ( pass <= paths ){
			int i = szxmultszy*paths*l+szxmultszy*(pass-1)+sizex*(int)ys+(int)xs;
			atomicAdd(&phi_x[i],(float)(cosf(r0*Kap)));
    			atomicAdd(&phi_y[i],(float)(sinf(r0*Kap)));
    			atomicAdd(&N[i],(unsigned long long)1);
		    }
		    global[t]=state;
		    return;
		    }    
		}

		if ( next==1 ){ // Reflecting from contact bottom
		// check if in gap
		if ( xs < floorf(xs/pxsize) * pxsize + 0.5*gap || 
                     xs > ceilf(xs/pxsize)  * pxsize - 0.5*gap || 
                     ys < floorf(ys/pxsize) * pxsize + 0.5*gap || 
                     ys > ceilf(ys/pxsize)  * pxsize - 0.5*gap ){ 
		    global[t]=state;
		    return;
		}
		else{
    		    costheta = curand_uniform(&state);				// New angle
		    sintheta = sqrt(1.0-costheta*costheta);
		    temp *= curand_uniform(&state);				// Extra call. States seem to be correlated
		    float gamma = twopi*curand_uniform(&state);
		    singamma = sinf(gamma);
		    cosgamma = cosf(gamma);
		    next     = 2;
		}
		}

		if ( next==2 ){ // Evolve up in absorption layer
		pass++;
		int l;

	        // Calculate distance for random absorption probability
		// curand_uniform does give 1.0, so inverting here.
	        r = -log(1.0-curand_uniform(&state))/alpha;
	        // That distance corresponds to the following z distance
	        z = r*costheta;
	        // if z > Dact, light either gets mirrored or goes intro transparent layer. Recalc r to transparent surface.
	        if ( z >= Dact ){
		    z = Dact;
		    r = Dact/costheta;
		}
		else{ // We are doing a single layer right now
		    l = floor(z*IRlay/Dact);
	    	    z = (l+0.5)*(Dact/IRlay);
		    r = (l+0.5)*(Dact/IRlay)/costheta;
		}
		// Calculate new x and y positions.
		float xo = xs;
		float yo = ys;
	        xs += r*sintheta*cosgamma;
	        ys += r*sintheta*singamma;

		// Recentering will help with calculating the diffraction correctly.
		if ( z < Dact ){
		    xs = (floorf(xs/pxoversubp)+0.5)*pxoversubp; 
		    ys = (floorf(ys/pxoversubp)+0.5)*pxoversubp; 
		    r  = sqrtf(powf(xo-xs,2.0)+powf(yo-ys,2.0)+powf(z,2.0));
		}
		// Find region in image. If off image, exit. If not, and absorbed, recenter
		// to pixel center.
		if ( xs < -halfx || xs >= halfx || ys < -halfy || ys >= halfy ){
		    global[t]=state;
		    return;
		}

	        r0 += r;
	        if ( z==Dact ){
		    if ( curand_uniform(&state) > Trans ){	//	Absorbed in contact
			global[t]=state;
			return;
		    }
		    else	next = 3;
	        }
		else{
		    xs = floorf(xs/pxoversubp) + szxover2; 
		    ys = floorf(ys/pxoversubp) + szyover2; 
		    if ( pass <= paths ){
			int i = szxmultszy*paths*l+szxmultszy*(pass-1)+sizex*(int)ys+(int)xs;
    			atomicAdd(&phi_x[i],(float)(cosf(r0*Kap)));
    			atomicAdd(&phi_y[i],(float)(sinf(r0*Kap)));
    			atomicAdd(&N[i],(unsigned long long)1);
		    }
		    global[t]=state;
		    return;
		}
		}

		if ( next==3 ){ // Evolve up in transparent layer

		// Check whether the photons exit or reflect. Will depend on 
		// internal incidence angle and index of refraction. 
		if ( sintheta >= 1.0/nsil ) next =  4;	// Definit total internal reflection
		else{ 
		    double costhetat = sqrt(1.0-nsil*nsil*sintheta*sintheta);
		    // Get random polarization angle
		    double polar = curand_uniform(&state) * 2.0 * M_PI;
		    double Eperp = cos(polar);
		    double Eparr = sin(polar);
		    double rperp = Eperp * (nsil*costheta - costhetat)/(nsil*costheta + costhetat);
		    double rparr = Eparr * (costheta - nsil*costhetat)/(nsil*costhetat + costheta);

		    double P = rperp * rperp + rparr * rparr;
		    if ( curand_uniform(&state) <= P ) next = 4;
		    else{ 
			global[t]=state;
			return; 	// Ray exits
		    }
		}

		// Calculate new x and y positions.
		r = (Top-Dact)/costheta;
	        xs += r*sintheta*cosgamma;
	        ys += r*sintheta*singamma;

		// Find region in image. If off image, exit. If not, and absorbed, recenter
		// to pixel center.
		if ( xs < -halfx || xs >= halfx || ys < -halfy || ys >= halfy ){
		    global[t]=state;
		    return;
		}
		r0   += r;
		}

		if ( next==4 ){ // Evolve down in transparent layer

		r   = (Top-Dact)/costheta;
	        xs += r*sintheta*cosgamma;
	        ys += r*sintheta*singamma;

		if ( xs < -halfx || xs >= halfx || ys < -halfy || ys >= halfy ){ 
		    global[t]=state; 
		    return;
		}

		r0 += r;

		if ( curand_uniform(&state) <= Trans )	next = 0; 	// Goes forward down
		else{ 		    global[t]=state;			// Absorbed in contact
				    return;
		}
		}

	    }while( next>=0 );

	} // End of main if clause

}	// End of kernel
                                                        