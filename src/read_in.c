
/******************************************************************************/
/* Configuration file processing program                                      */
/* Much thanks to Andras Pal for source codes on character processing         */
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "imMIRI.h"

void exit_with_usage(int);

int char_is_space(int c)
{
 if ( c==64 || c==32 || c==13 || c==10 || c==9 )        return(1);
 else                                   return(0);
}

void remove_quotes(char *buff)
{
 int k;
 while ( *buff )
  {     for ( k=0 ; buff[k]=='"' ; )    k++;
        if ( k )        memmove(buff,buff+k,strlen(buff)+1-k);
        else            buff++;
  }
}

int tokenize_spaces(char *buff,char **tokens,int max)
{
 int    intoken,inquota,n;
 char   **tsave;

 tsave=tokens;

 intoken=0,inquota=0;n=0;
 while ( *buff && n<max )
  {     if ( ( ! char_is_space(*buff) ) && ! intoken )
         {      *tokens=buff;
                intoken=!0,inquota=0;n++;
                if ( *buff=='"' )       inquota=!0;
                tokens++,buff++;
         }
        else if ( intoken && ( (char_is_space(*buff) && inquota) || (!char_is_space(*buff)) ) )
         {      if ( *buff=='"' )       inquota=!inquota;
                buff++;
         }
        else if ( intoken && ! inquota && char_is_space(*buff) )
         {      *buff=0,buff++;
                intoken=0;
         }
        else    buff++;
  };
 *tokens=NULL;

 while ( *tsave != NULL )
  {     remove_quotes(*tsave);
        tsave++;
  };

 return(n);
}

int read_in_param(char *param){
 FILE        *fr;
 int         i,j,t,itemp;
 long	     ltemp;
 float	     dtemp;
 char        buff[512],*dat[64];
 configdata  *cfg;

 cfg=&cconfig;

 fr=fopen(param,"rb");
 if ( fr==NULL )        exit_with_usage(3);

 i=0;
 while ( ! feof(fr) ){
    if ( fgets(buff,512,fr)==NULL )     break;
    if ( buff[0]=='#' )                 continue;
    t=tokenize_spaces(buff,dat,9);
    if ( t==0 )				continue;
    if ( i==0 ){											// verbose
	if ( strcmp(dat[0],"yes")==0 || strcmp(dat[0],"Y")==0 || strcmp(dat[0],"YES")==0
	    || strcmp(dat[0],"y")==0 || strcmp(dat[0],"true")==0 || strcmp(dat[0],"TRUE")==0 )	
		cfg->verb=1;
	else if ( strcmp(dat[0],"no")==0 || strcmp(dat[0],"N")==0 || strcmp(dat[0],"NO")==0
	    || strcmp(dat[0],"n")==0 || strcmp(dat[0],"false")==0 || strcmp(dat[0],"FALSE")==0 )
		cfg->verb=0;
	else exit_with_usage(4);
    }
    if ( i==1 ){											// Light source
	if ( strcmp(dat[0],"fits")==0 || strcmp(dat[0],"FITS")==0 ){
	    cfg->method=0;
	    cfg->input=(char *)malloc(strlen(dat[1])+1);
	    strcpy(cfg->input,dat[1]);
	    cfg->D=0.0;
	    cfg->pxsc=0.0;
	    cfg->r=0.0;
	    if ( t==2){
	        cfg->def=0;
	        cfg->lam=5.6;			// Will get overwritten when file is read in
	        cfg->oversamp=1.0;		// if in header
	    }
	    else if ( t==4 ){
	        cfg->def=1;
	        sscanf(dat[2],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(7);
	        else cfg->lam=dtemp;
	        sscanf(dat[3],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(8);
	        else cfg->oversamp=dtemp;
	    }
	    else exit_with_usage(6);
    	}
	else if ( strcmp(dat[0],"airy")==0 || strcmp(dat[0],"AIRY")==0 ){
	    cfg->method=1;
	    cfg->input=NULL;
	    cfg->r=0;
	    cfg->oversamp=1.0;
	    if ( t==1){
	        cfg->def=0;
	        cfg->lam=5.6;
	        cfg->D=6.5;
	        cfg->pxsc=0.11;
	    }
	    else if ( t==4 ){
	        cfg->def=1;
	        sscanf(dat[1],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(9);
	        else cfg->D=dtemp;
	        sscanf(dat[2],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(10);
	        else cfg->pxsc=dtemp;
	        sscanf(dat[3],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(7);
	        else cfg->lam=dtemp;
	    }
	    else exit_with_usage(11);
    	}
	else if ( strcmp(dat[0],"disk")==0 || strcmp(dat[0],"DISK")==0 ){
	    cfg->method=2;
	    cfg->input=NULL;
	    cfg->oversamp=1.0;
	    cfg->pxsc=0.0;
	    cfg->D=0.0;
	    cfg->def=1;
	    if ( t==1){
	        cfg->lam=5.6;
	        cfg->r=25.0;
	    }
	    else if ( t==3){
	        sscanf(dat[1],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(12);
	        else cfg->r=dtemp;
	        sscanf(dat[2],"%f",&dtemp);
	        if ( dtemp <=0 ) exit_with_usage(7);
	        else cfg->lam=dtemp;
	    }
	    else exit_with_usage(13);
    	}
	else exit_with_usage(5);
    }
    if ( i==2 ){
	if ( t==2 ){
	    cfg->xoff=0.0;
	    sscanf(dat[0],"%f",&dtemp);
	    cfg->xoff=dtemp;
	    cfg->yoff=0.0;
	    sscanf(dat[1],"%f",&dtemp);
	    cfg->yoff=dtemp;
	}
	else exit_with_usage(25);
    }
    if ( i==3 ){ 
	if ( t==3 ){
	    sscanf(dat[0],"%f",&dtemp); 
	    if ( dtemp <=0 ) exit_with_usage(14);
	    else cfg->Nstart=(long)dtemp;
	    sscanf(dat[1],"%f",&dtemp); 
	    if ( dtemp <=0 ) exit_with_usage(15);
	    else cfg->Nmult=dtemp;
	    sscanf(dat[2],"%f",&dtemp); 
	    if ( dtemp <=0 ) exit_with_usage(22);
	    else cfg->Nlast=(long)dtemp;
	}
	    else exit_with_usage(23);
    }
    if ( i==4 ){											// stubs
	cfg->output=(char *)malloc(strlen(dat[0])+1);
	strcpy(cfg->output,dat[0]);
    }    
    if ( i==5 ){
	cfg->ndev = t;
	cfg->devs = (int *)malloc(sizeof(int)*t);
	for ( j=0 ; j<t ; j++ ){
	    sscanf(dat[j],"%d",&itemp);
	    if ( itemp <0 ) exit_with_usage(24);
	    else cfg->devs[j]=itemp;	    
	}
    }    
 i++;
 }
 if ( i<5 )            exit_with_usage(16);
 
// cfg->alpha = 0.0102 * pow(cfg->lam/7.0,2.0); // v1
// cfg->alpha = 0.0112 * pow(cfg->lam/7.0,2.0) * ( -0.0003 * pow(cfg->lam,3.0) + 0.0164 * pow(cfg->lam,2.0) - 0.3057 * cfg->lam + 2.9069 ); // v2
 cfg->alpha = 0.0890*(-0.0000840284*pow(cfg->lam,3.0)+0.007094981254*pow(cfg->lam,2.0)-0.035111701146*cfg->lam+0.093887179734);
// cfg->Trans = exp(-1.0e4*cfg->alpha*0.00045)-0.01; // v1
// cfg->Trans = 0.67 + 0.33*exp(-1.0e4*cfg->alpha*0.001409); // v2
 cfg->Trans = pow((0.01*(0.0010429704*pow(cfg->lam,3.0)-0.0824166173*pow(cfg->lam,2.0)+0.7497968952*cfg->lam+97.6969896581)),2.8);
 cfg->lam  /= nsil;
 cfg->Kap   = 2.0*M_PI/cfg->lam;
 cfg->nx1   = 0;
 cfg->nx2   = 0;
 cfg->subp  = 1;
 cfg->sizex = 0;
 cfg->sizey = 0;
 cfg->Tot   = 0;

 fclose(fr);
 
 return 0;
}
                                                  