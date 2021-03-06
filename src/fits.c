#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fitsio.h"
#include "imMIRI.h"

void print_warning(int);

double Fairy(double x1,double x2,double y1,double y2){	// coordinates in pixels
 configdata	*cfg	= &cconfig;
 double c,R,Iz;
 double Fx1y1,Fx1y2,Fx2y1,Fx2y2;

 c = M_PI * cfg->D * 1.0e6 * cfg->pxsc * M_PI / (180.0*3600.0) / ( cfg->lam * nsil );

 Iz = M_PI*pow(cfg->D*1.0e6*0.5*M_PI*cfg->pxsc/cfg->lam/nsil/pxsize/180.0/3600.0,2.0);

 R = sqrt(x1*x1+y1*y1);
 if (R!=0) Fx1y1 = Iz*pow(2.0*j1(c*R)/(c*R),2.0);
 else      Fx1y1 = Iz;
 R = sqrt(x2*x2+y1*y1);
 if (R!=0) Fx2y1 = Iz*pow(2.0*j1(c*R)/(c*R),2.0);
 else      Fx2y1 = Iz;
 R = sqrt(x1*x1+y2*y2);
 if (R!=0) Fx1y2 = Iz*pow(2.0*j1(c*R)/(c*R),2.0);
 else      Fx1y2 = Iz;
 R = sqrt(x2*x2+y2*y2);
 if (R!=0) Fx2y2 = Iz*pow(2.0*j1(c*R)/(c*R),2.0);
 else      Fx2y2 = Iz;

 // return(pow(2.0*j1(c*R)/(c*R),2.0)*(c*pxsize*c*pxsize));
 return(0.25*(Fx1y1+Fx2y1+Fx1y2+Fx2y2)*pxsize*pxsize*0.01*0.01);
}

void getimdata(){
    fitsfile   *fptr;       		/* pointer to the FITS file, defined in fitsio.h */
    configdata *cfg;
    float     dtemp = 0;
    int status;
    char comment[FLEN_COMMENT];   	/* standard string lengths defined in fitsioc.h */

    cfg = &cconfig;

    status = 0;

    fits_open_file(&fptr,cfg->input, READONLY, &status);
    fits_report_error(stderr, status);
    fits_read_key(fptr,TINT,"NAXIS1",&cfg->nx1,comment,&status);
    fits_report_error(stderr, status);
    fits_read_key(fptr,TINT,"NAXIS2",&cfg->nx2,comment,&status);
    fits_report_error(stderr, status);

    if ( cfg->def==0 ){ 
	fits_read_key(fptr,TFLOAT,"DET_SAMP",&dtemp,comment,&status);
        fits_report_error(stderr, status);
	if ( status==202 || status==204 ) print_warning(17);
	else cfg->oversamp=dtemp;
    }
    else ;
    if ( cfg->def==0 ){ fits_read_key(fptr,TFLOAT,"WAVELEN",&dtemp,comment,&status);
        if ( status==202 || status==204 ) print_warning(18);
	else{
	 cfg->lam=(dtemp)*1.0e6;
//	 cfg->alpha = 0.0102 * pow(cfg->lam/7.0,2.0);
	 cfg->alpha = 0.0112 * pow(cfg->lam/7.0,2.0) * ( -0.0003 * pow(cfg->lam,3.0) + 0.0164 * pow(cfg->lam,2.0) - 0.3057 * cfg->lam + 2.9069 );
	 cfg->lam  /= nsil;
	 cfg->Kap   = 2.0*M_PI/cfg->lam;
	}
    }
    else ;
    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);
}

void fillinput(input_h *rpsf){
    // Crop image? Or only get flux from pixels not cropped?
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,anynull;
    long fpixel,nbuffer,i;
    int npix,k,x,y,xs,xe,ys,ye;
    int nx1_orig,nx2_orig;
    double nullval,buffer[buffsize];

    configdata *cfg = &cconfig;
    input_h psf     = *rpsf;
    status          = 0;
    nx1_orig        = cfg->nx1;
    nx2_orig        = cfg->nx2;
    
    if ((int)(imsizex*cfg->oversamp)<cfg->nx1 ) cfg->nx1 = (int)(imsizex*cfg->oversamp);
    if ((int)(imsizey*cfg->oversamp)<cfg->nx2 ) cfg->nx2 = (int)(imsizey*cfg->oversamp);
    npix = nx1_orig*nx2_orig;
    if (!psf.x)  psf.x  = (int *)malloc(sizeof(int)*cfg->nx1*cfg->nx2);
    if (!psf.y)  psf.y  = (int *)malloc(sizeof(int)*cfg->nx1*cfg->nx2);
    if (!psf.F)  psf.F  = (double *)malloc(sizeof(double)*cfg->nx1*cfg->nx2);
    if (!psf.F0) psf.F0 = (double *)malloc(sizeof(double)*cfg->nx1*cfg->nx2);

    xs = (int)((nx1_orig-cfg->nx1)*0.5);
    ys = (int)((nx2_orig-cfg->nx2)*0.5);
    xe = xs + cfg->nx1;
    ye = ys + cfg->nx2;

    fits_open_file(&fptr, cfg->input, READONLY, &status);

    fpixel   = 1;
    nullval  = 0;

    while (npix > 0)
    {
      nbuffer = npix;
      if (npix > buffsize) nbuffer = buffsize;

      fits_read_img(fptr,TDOUBLE,fpixel,nbuffer,&nullval,buffer,&anynull,&status);

      // naxis1 - number of columns (x elements); naxis2 - number of rows (y elements)
      // I use my own standard of starting pixels at zero (sorry, C rules)
      for (i = 0; i < nbuffer; i++){
	x = i+fpixel-1 - (long int)((i+fpixel-1)/nx1_orig)*nx1_orig;
	y = (long int)((i+fpixel-1)/nx1_orig);
	if ( x>=xs && x<xe && y>=ys && y<ye ){	
	    k = cfg->nx1*(y-ys)+(x-xs);
	    psf.x[k]  = x-xs;
	    psf.y[k]  = y-ys;
	    psf.F[k]  = buffer[i];
	    psf.F0[k] = buffer[i];
	}
      }
      npix    -= nbuffer;
      fpixel  += nbuffer;
    }
    fits_close_file(fptr, &status);    

    *rpsf = psf;
}

void write_fits(flux img,input_h psf_h,long j,const char * Version){
    fitsfile *fptr;					/* pointer to the FITS file; defined in fitsio.h */
    configdata *cfg = &cconfig;
    char *fname;
    char *ver;
    int status,i,x,y,l,p,xs,ys,jlen;
    long naxis = 2;
    long naxes[2];
    long fpixel = 1;
    unsigned long long Ntot;
    float lam;
    double Ptot;
    double *F[imsizey];
    double R,Flux;
    double pathlaytot;
    double x1,y1;
    float  ftemp;
    int	   itemp;
    status = 0; 					/* initialize status before calling fitsio routines */

    jlen = ( j == 0 ? 1 : (int)(log10(j)+1));
    fname=(char *)malloc(strlen(cfg->output)+jlen+7);
    ver=(char *)malloc(strlen(Version)+1);
    strcpy(ver,Version);
    sprintf(fname,"%s_%ld.fits",cfg->output,j);
    naxes[0]=imsizex;
    naxes[1]=imsizey;

    remove(fname);

    fits_create_file(&fptr,fname, &status); 		/* create new file */
    fits_report_error(stderr, status);                  /* print out any error messages */
    fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    fits_report_error(stderr, status);                  /* print out any error messages */
    
    F[0] = (double *)malloc(sizeof(double)*imsizex*imsizey);
    for ( y=1 ; y < imsizey ; y++ ) F[y] = F[y-1] + imsizex;
    Ptot=0.0;

    for ( y=0 ; y<imsizey ; y++ ){
        for ( x=0 ; x<imsizex ; x++){
		F[y][x] = 0;
    }}

    double abs = cfg->Trans * exp(-cfg->alpha*Dact);

    for ( l=0 ; l < IRlay ; l++ ){
	// get path/layer total first
	for ( p=0 ; p<paths ; p++ ){
	    pathlaytot = 0.0;
	    Ntot   = 0;
	    for ( ys=0 ; ys<imsizey*cfg->subp ; ys++ ){
		for ( xs=0 ; xs<imsizex*cfg->subp ; xs++ ){
		    i = cfg->sizex * cfg->sizey * paths * l + cfg->sizex * cfg->sizey * p + cfg->sizex * ys + xs;
		    pathlaytot += pow((double)img.phi_x[i],2.0)+pow((double)img.phi_y[i],2.0);
		    Ntot   += img.N[i];
		}
	    }
	    if (Ntot!=0){
	    for ( y=0 ; y<imsizey ; y++ ){
		for ( x=0 ; x<imsizex ; x++){
		    for ( ys=y*cfg->subp ; ys<(y+1)*cfg->subp ; ys++ ){
			for ( xs=x*cfg->subp ; xs<(x+1)*cfg->subp ; xs++ ){
			    i = cfg->sizex * cfg->sizey * paths * l +cfg->sizex * cfg->sizey * p + cfg->sizex * ys + xs;
		    	    F[y][x] += ((double)Ntot)/((double)j/abs)*(pow((double)img.phi_x[i],2.0)+pow((double)img.phi_y[i],2.0))/pathlaytot;
			}
		    }
		}
	    }
	    }
	}
    }

    // Add additional flux absorbed in first passing

    abs = cfg->Trans * (1.0-exp(-cfg->alpha*Dact));

    double xoff = 0.5*(imsizex-cfg->nx1/cfg->oversamp);
    double yoff = 0.5*(imsizex-cfg->nx1/cfg->oversamp);

    if ( cfg->method == 0 ){
	for ( i=0 ; i < cfg->nx1*cfg->nx2 ; i++){
	    x = (int)(psf_h.x[i]/cfg->oversamp+xoff);
	    y = (int)(psf_h.y[i]/cfg->oversamp+yoff);
	    F[y][x] += psf_h.F0[i]*abs;
	}
	for ( y=0 ; y<imsizey ; y++ ){
	    for ( x=0 ; x<imsizex ; x++){
		F[y][x]*=cfg->Tot;
	    }
	}
    }

    if ( cfg->method == 1 ){
	for ( y=0 ; y<imsizey ; y++ ){
	    for ( x=0 ; x<imsizex ; x++ ){
		// Need to subpixel for smooth image - Going by 0.01 px;
		for ( xs=0 ; xs<100 ; xs++ ){
		    for ( ys=0 ; ys<100 ; ys++ ){
			// Get Airy Flux
		    	x1=(double)(x+xs*0.01);
		    	y1=(double)(y+ys*0.01);
			F[y][x] += abs*Fairy(x1-imsizex/2,x1+0.01-imsizex/2,y1-imsizey/2,y1-imsizey/2+0.01);
		    }
		}
	    }
	}
    }

    if ( cfg->method == 2 ){
	abs*=pxsize*pxsize*0.01*0.01/pow(cfg->r,2.0)/M_PI;
	for ( y=0 ; y<imsizey ; y++ ){
	    for ( x=0 ; x<imsizex ; x++){
		for ( xs=0 ; xs<100 ; xs++ ){
		    for ( ys=0 ; ys<100 ; ys++ ){
		    	x1=(double)(x+xs*0.01);
		    	y1=(double)(y+ys*0.01);
			// Find R
			R = sqrt(pow(((double)x1+0.005)*pxsize-halfx,2.0)+pow(((double)y1+0.005)*pxsize-halfy,2.0));
			// If in pin, add flux
			if ( R <= cfg->r ) F[y][x] += abs;
		    }
		}
	    }
	}
    }

    // Add up total flux

    for ( y=0 ; y<imsizey ; y++ ){
	for ( x=0 ; x<imsizex ; x++){
	    Ptot+=F[y][x];
	}
    }

    lam = cfg->lam * nsil;
    ftemp = (float)pxsize;
    fits_update_key(fptr, TFLOAT, "PXSIZE",&ftemp,"[micron] Pixel pitch", &status);
    fits_update_key(fptr, TFLOAT, "WAVELENG", &lam,"[micron] Wavelength of dataset", &status);
    fits_update_key(fptr, TDOUBLE, "FLUX", &Ptot,"Total flux detected", &status);
    fits_update_key(fptr, TLONG, "NUMPHOT", &j,"Number of photons in model", &status);
    if ( cfg->method ==0 ){
	fits_update_key(fptr, TSTRING, "METHOD", "Fits input","Method used for model", &status);
	fits_update_key(fptr, TSTRING, "ORIGFILE", cfg->input,"Original input file", &status);
	fits_update_key(fptr, TFLOAT, "OVERSAMP",&cfg->oversamp,"Oversampling in original input file", &status);
	fits_update_key(fptr, TDOUBLE, "ORIGFLUX",&cfg->Tot,"Total flux in original (cropped) image", &status);
	}
    if ( cfg->method ==1 ){
	fits_update_key(fptr, TSTRING, "METHOD", "Airy disk","Method used for model", &status);
	fits_update_key(fptr, TFLOAT, "DIAMETER",&cfg->D,"[m] Diameter of input telescope", &status);
	fits_update_key(fptr, TFLOAT, "PIXELSCA",&cfg->pxsc,"[''/px] pixel scale of instrument setup", &status);
    }
    if ( cfg->method ==2 ){
	fits_update_key(fptr, TSTRING, "METHOD", "Pinhole disk","Method used for model", &status);
	fits_update_key(fptr, TFLOAT, "RADIUS",&cfg->r,"[micron] Radius of input pinhole", &status);
    }
    ftemp = (float)offtheta;
    fits_update_key(fptr, TFLOAT, "OFFTHETA",&ftemp,"[rad] Off-axis theta value", &status);
    ftemp = (float)offgamma;
    fits_update_key(fptr, TFLOAT, "OFFGAMMA",&ftemp,"[rad] Off-axis gamma value", &status);
    fits_update_key(fptr, TFLOAT, "XOFF",&cfg->xoff,"[px] x direction dither", &status);
    fits_update_key(fptr, TFLOAT, "YOFF",&cfg->yoff,"[px] y direction dither", &status);
    ftemp = (float)Dact;
    fits_update_key(fptr, TFLOAT, "DACT",&ftemp,"[micron] IR-active layer width in model", &status);
    ftemp = (float)Top;
    fits_update_key(fptr, TFLOAT, "TOP",&ftemp,"[micron] Total width of the detector wafer", &status);
    ftemp = (float)gap;
    fits_update_key(fptr, TFLOAT, "GAP",&ftemp,"[micron] Gap size between metallic contacts", &status);
    fits_update_key(fptr, TINT, "SUBPIX",&cfg->subp,"Number of subpixels per pixel in model", &status);
    itemp = (int)paths;
    fits_update_key(fptr, TINT, "PATHS",&itemp,"Number of photon paths used by model", &status);
    ftemp = (float)cfg->Trans;
    fits_update_key(fptr, TFLOAT, "TRANS",&ftemp,"Transmission efficiency of transparent contact", &status);
    ftemp = (float)nsil;
    fits_update_key(fptr, TFLOAT, "NREFRACT",&ftemp,"Refractive index of silicon used in model", &status);
    fits_update_key(fptr, TSTRING, "VERSION", ver,"G-imMIRI version", &status);
    fits_write_history(fptr,"Created with G-imMIRI",&status);
    fits_write_date(fptr,&status);

    fits_write_img(fptr, TDOUBLE, fpixel, imsizex*imsizey, F[0], &status);
    fits_report_error(stderr, status);                  /* print out any error messages */
    fits_close_file(fptr, &status); 			/* close the file */
    fits_report_error(stderr, status);                  /* print out any error messages */
    free(F[0]);
    free(ver);

    if ( cfg->verb!=0 ) printf("\nFile %s successfully written\n",fname);
}
                                                                       