#ifndef __IMMIRI_H_INCLUDED
#define __IMMIRI_H_INCLUDED

/****************************************************************************************/
/*                                                                                      */
/* The constants used in the code. 							*/
/* (if constant is a float, use decimal places)						*/
/*											*/
/****************************************************************************************/
#define Dact     35.0   /* Width of active layer in microns				*/
#define Top      500.0  /* Location of the top						*/
#define gap      1.3    /* Assuming a 1.3 micron gap between pixel contacts		*/
#define pxsize   25.0   /* micron pitch of MIRI detector pixels				*/
#define subpix	 1.0	/* Sub-pixeling in microns in pixels 				*/
			/* Sub-pixeling will be rounded to nearest integer 		*/
			/* to yield subpix closest to one defined here			*/
			/* Should be <= than wavelength in material			*/
#define paths	 10	/* Number of photon path-types to trace				*/
			/* The path-types are ordered by rough likelihood               */
			/* First straight down path not included in this number,	*/
			/* that is automatic						*/
#define IRlay	 1	/* Number of IR sub-layers to be considered			*/
#define nsil	 3.4127 /* Refractive index of silicon					*/
#define imsizex  200	/* Size of final image - # of columns - give even number	*/
#define imsizey  200	/* Size of final image - # of rows - give even number		*/
#define offtheta 7.8e-4	/* Detector off axis in theta (radial axis angle) [rad]		*/
#define offgamma -0.24  /* Detector off axis in gamma (CCW from x=0) [rad]		*/
#define buffsize 1000   /* read in fits image in small chunks				*/
#define cycles	 320	/* Thread count multiplier for cycles				*/

#define halfx    imsizex*pxsize*0.5
#define halfy    imsizey*pxsize*0.5
#define twopi	 6.28318530718

/****************************************************************************************/
/*                                                                                      */
/* The structures used in the code. 							*/
/*											*/
/****************************************************************************************/

/* Pixels and layer									*/
typedef struct{
 float 	  *phi_x,*phi_y;	/* Phaser amplitudes; index=sx*sy*l*p+sx*sy*p+x*sy+y	*/
 unsigned long long  *N;	/* Number of photons captured per pixel			*/
}flux;

/* On host										*/
typedef struct{
 double 	*F,*F0;		/* Input fits file fluxes; gets replaced by CDF		*/
 int		*x,*y;		/* Input fits file coordinates				*/
}input_h;

/* On device(s)										*/
typedef struct{
 double 	*F;		/* Input fits file fluxes in CDF;			*/
 int		*x,*y;		/* Input fits file coordinates				*/
}input_d;

typedef struct{
 int		verb;		/* verbose? 0=no, 1=yes					*/
 int		method;		/* 0=fits, 1=airy, 2=disk				*/
 int		def;		/* 0=default/header, 1=force				*/
 char		*input;		/* if fits: filename w/ opt. ext. []			*/	
 char		*output;	/* output fits filename					*/
 float		lam;		/* Observation wavelength				*/
 float		oversamp;	/* Oversamp rate in input fits file			*/
 float		D;		/* if airy: Tel. aper (m)				*/
 float		pxsc;		/* if airy: pixelscale (''/px)				*/
 float		r;		/* if disk: radius (in mu)				*/
 long 		Nstart;		/* First output at					*/
 float		Nmult;		/* output factor					*/
 long		Nlast;		/* last	output file at					*/
 double		Tot;		/* Total original flux in system			*/
 int            sizex,sizey;	/* subpix image size					*/
 int		subp;		/* Number of subpixels per pixel			*/
 float 		alpha;		/* Absorption in silicon active layer			*/
 float		Kap;		/* Kappa 						*/
 float		Trans;		/* Transparency factor of buried contact		*/
 int		*MPcount;	/* number of MP in device				*/
 int		*cores;		/* number of cores per MP in device			*/
 int		nx1,nx2;	/* input PSF file dimensions				*/
 int		ndev;		/* number of CUDA devices				*/
 int		*devs;		/* array of CUDA device numbers				*/
 float		xoff,yoff;	/* Given offset in microns				*/
} configdata;

extern	   configdata 	cconfig;

void   checkmem();
void   reset_cuda();
void   init_psf(input_h *);
void   init_img(flux *,flux **);
void   sortinput(input_h *,input_d **);
void   readbackfromdevices(flux,flux *);
void   fill_const();

#endif
                                                            