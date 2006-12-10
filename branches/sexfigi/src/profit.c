 /*
 				profit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Fit an arbitrary profile combination to a detection.
*
*	Last modify:	10/12/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"lmfit/lmmin.h"
#include	"fft.h"
#include	"check.h"
#include	"psf.h"
#include	"profit.h"


static double	prof_interpolate(profstruct *prof, double *posin);
static double	interpolate_pix(double *posin, double *pix, int *naxisn,
		interpenum interptype);

static void	make_kernel(double pos, double *kernel, interpenum interptype);

/*------------------------------- variables ---------------------------------*/

int		interp_kernwidth[5]={1,2,4,6,8};
int theniter;
/* "Local" global variables; it seems dirty but it simplifies a lot */
/* interfacing to the LM routines */
static objstruct	*the_obj;
static picstruct	*the_field, *the_wfield;
profitstruct		*theprofit;

/****** profit_init ***********************************************************
PROTO	profitstruct profit_init(psfstruct *psf,
			proftypenum *profcode, int nprof)
PURPOSE	Allocate and initialize a new profile-fitting structure.
INPUT	Pointer to PSF structure,
	pointer to the list of profile component type codes,
	number of profile types.
OUTPUT	A pointer to an allocated profit structure.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2006
 ***/
profitstruct	*profit_init(psfstruct *psf, proftypenum *profcode, int nprof)
  {
   profitstruct		*profit;
   int			p;

  QCALLOC(profit, profitstruct, 1);
  profit->nprof = nprof;
  profit->nparam = 11;
  QMALLOC(profit->param, double, profit->nparam);
  QMALLOC(profit->paraminit, double, profit->nparam);
  QMALLOC(profit->parammin, double, profit->nparam);
  QMALLOC(profit->parammax, double, profit->nparam);
  QMALLOC(profit->prof, profstruct *, nprof);

  for (p=0; p<nprof; p++)
    profit->prof[p] = prof_init(profit, profcode[p]);

  profit->psf = psf;
  profit->psfdft = NULL;

  return profit;
  }  


/****** profit_end ************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile-fitting structure.
INPUT	Prof structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/12/2006
 ***/
void	profit_end(profitstruct *profit)
  {
   int	p;

  for (p=0; p<profit->nprof; p++)
    prof_end(profit->prof[p]);
  free(profit->prof);
  free(profit->param);
  free(profit->paraminit);
  free(profit->parammin);
  free(profit->parammax);
  free(profit);

  return;
  }


/****** prof_fit ************************************************************
PROTO	void prof_fit(profitstruct *profit, picstruct *field, picstruct *wfield,
		objstruct *obj, obj2struct *obj2)
PURPOSE	Fit profile(s) convolved with the PSF to a detected object.
INPUT	Array of profile structures,
	Number of profiles,
	Pointer to the profile-fitting structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Pointer to an allocated fit structure (containing details about the
	fit).
NOTES	It is a modified version of the lm_minimize() of lmfit.
AUTHOR	E. Bertin (IAP)
VERSION	09/12/2006
 ***/
void	prof_fit(profitstruct *profit,
		picstruct *field, picstruct *wfield,
		objstruct *obj, obj2struct *obj2)
  {
    psfstruct		*psf;
    lm_control_type	control;
    checkstruct		*check;
    double		*diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4,
			psf_fwhm, scaling;
    int			*ipvt,
			ix,iy, m,n,p;


  if (profit->psfdft)
    {
    QFREE(profit->psfdft);
    }

  psf = profit->psf;

/* Compute the local PSF */
  psf_build(psf);


/* Create pixmaps at image resolution */
  psf_fwhm = psf->masksize[0]*psf->pixstep;
  profit->objnaxisn[0] = (int)((obj->xmax - obj->xmin)*2.0 + psf_fwhm + 0.499);
  profit->objnaxisn[1] = (int)((obj->ymax - obj->ymin)*2.0 + psf_fwhm + 0.499);
  ix = (int)(obj2->posx+0.49999);
  iy = (int)(obj2->posy+0.49999);
  QMALLOC(profit->objpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  QMALLOC(profit->lmodpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  profit->nresi = profit_copyobjpix(profit, field, ix,iy);
  if (profit->nresi >= profit->nparam)
    {
    QMALLOC(profit->resi, double, profit->nresi);
    }
  else
    {
    for (p=0; p<profit->nparam; p++)
      obj2->prof_vector[p] = 0.0;
    obj2->prof_niter = 0;
    return;
    }

/* Create pixmap at PSF resolution */
  profit->modnaxisn[0] = (int)(profit->objnaxisn[0] / psf->pixstep +0.4999); 
  profit->modnaxisn[1] = (int)(profit->objnaxisn[1] / psf->pixstep +0.4999); 
if (profit->modnaxisn[1] < profit->modnaxisn[0])
profit->modnaxisn[1] = profit->modnaxisn[0] = 256;
else
profit->modnaxisn[0] = profit->modnaxisn[1] = 256;

  QCALLOC(profit->modpix, double, profit->modnaxisn[0]*profit->modnaxisn[1]);
  scaling = psf->pixstep / profit->prof[0]->typscale;

/* Set initial guesses and boundaries */
  profit->paraminit[0] = 0.0;
  profit->parammin[0] = -6.0*field->backsig;
  profit->parammax[0] = 6.0*field->backsig;

  profit->paraminit[1] = 0.0;
  profit->parammin[1] = (double)(obj->xmin - ix);
  profit->parammax[1] = (double)(obj->xmax - ix);

  profit->paraminit[2] = 0.0;
  profit->parammin[2] = (double)(obj->ymin - iy);
  profit->parammax[2] = (double)(obj->ymax - iy);

  profit->paraminit[3] = obj->peak/2.0;
  profit->parammin[3] = 0.0;
  profit->parammax[3] = 1000.0*obj->peak;

  profit->paraminit[4] = obj->a*scaling;
  profit->parammin[4] = 0.0;
  profit->parammax[4] = 3.0*profit->paraminit[4];

  profit->paraminit[5] = obj->b*scaling;
  profit->parammin[5] = 0.0;
  profit->parammax[5] = 3.0*profit->paraminit[5];

  profit->paraminit[6] = obj->theta;
  profit->parammin[6] = 0.0;
  profit->parammax[6] = 180.0;

  profit->paraminit[7] = obj->peak/2.0;
  profit->parammin[7] = 0.0;
  profit->parammax[7] = 1000.0*obj->peak;

  profit->paraminit[8] = obj->b*scaling;
  profit->parammin[8] = 0.0;
  profit->parammax[8] = 3.0*profit->paraminit[8];

  profit->paraminit[9] = obj->b*scaling;
  profit->parammin[9] = 0.0;
  profit->parammax[9] = 3.0*profit->paraminit[9];

  profit->paraminit[10] = obj->theta;
  profit->parammin[10] = 0.0;
  profit->parammax[10] = 180.0;
/*
  profit->paraminit[11] = 4.0;
  profit->parammin[11] = 0.2;
  profit->parammax[11] = 8.2;
*/
  profit->sigma = field->backsig;

/* Use (dirty) global variables to interface with lmfit */
  the_field = field;
  the_wfield = wfield;
  the_obj = obj;

/* Allocate work space */
  n = profit->nparam;
  m = profit->nresi;

  control.ftol =      1.0e-7;
  control.xtol =      1.0e-7;
  control.gtol =      1.0e-7;
  control.maxcall =   PROFIT_MAXITER;
  control.epsilon =   1.0e-2;
  control.stepbound = 100.0;

  QMALLOC(diag, double,n);
  QMALLOC(qtf, double, n);
  QMALLOC(fjac, double,n*m);
  QMALLOC(wa1, double, n);
  QMALLOC(wa2, double, n);
  QMALLOC(wa3, double, n);
  QMALLOC(wa4, double,   m);
  QMALLOC(ipvt, int,   n);

/* Perform fit */
  control.info = 0;
  control.nfev = 0;

/* This goes through the modified legacy interface */
  lm_lmdif( m, n, profit->paraminit, profit->resi,
	control.ftol, control.xtol, control.gtol,
	control.maxcall*(n+1), control.epsilon, diag, 1,
	control.stepbound, &(control.info),
	&(control.nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
	profit_evaluate, profit_printout, NULL);

  control.fnorm = lm_enorm(m, profit->resi);
  if (control.info < 0 )
    control.info = 10;
printf("%d\n", control.info);
/* CHECK-Images */
  if ((check = prefs.check[CHECK_SUBPROFILES]))
    {
    profit_residuals(profit,field,wfield,obj,profit->paraminit,profit->resi);
    addcheck(check, profit->lmodpix, profit->objnaxisn[0],profit->objnaxisn[1],
		ix,iy, -1.0);
    }
  if ((check = prefs.check[CHECK_PROFILES]))
    {
    profit_residuals(profit,field,wfield,obj,profit->paraminit,profit->resi);
    addcheck(check, profit->lmodpix, profit->objnaxisn[0],profit->objnaxisn[1],
		ix,iy, 1.0);
    }
  if (FLAG(obj2.prof_vector))
    {
    for (p=0; p<profit->nparam; p++)
      obj2->prof_vector[p]= profit->paraminit[p];
    obj2->prof_vector[4] *= profit->prof[0]->typscale / psf->pixstep;
    obj2->prof_vector[5] *= profit->prof[0]->typscale / psf->pixstep;
    obj2->prof_vector[8] *= profit->prof[1]->typscale / psf->pixstep;
    obj2->prof_vector[9] *= profit->prof[1]->typscale / psf->pixstep;
    obj2->prof_niter = control.nfev;
    }

/* clean up. */
  free(diag);
  free(qtf); 
  free(fjac);
  free(wa1); 
  free(wa2); 
  free(wa3 );
  free(wa4); 
  free(ipvt);
  free(profit->modpix);
  free(profit->lmodpix);
  free(profit->objpix);
  free(profit->resi);

  return;
  }


/****** profit_printout *******************************************************
PROTO	void profit_printout(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev )
PURPOSE	Provide a function to print out results to lmfit.
INPUT	Number of fitted parameters,
	pointer to the vector of parameters,
	number of data points,
	pointer to the vector of residuals (output),
	pointer to the data structure (unused),
	0 (init) 1 (outer loop) 2(inner loop) -1(terminated),
	outer loop counter,
	number of calls to evaluate().
OUTPUT	-.
NOTES	Input arguments are there only for compatibility purposes (unused)
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2006
 ***/
void	profit_printout(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev )
  {
/*
printf("%d:	%g	%g	%g	%g	%g	%g	%g\n",
	iter, par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
*/
  if (0)
    lm_print_default(n_par, par, m_dat, fvec, data, iflag, iter, nfev);

  return;
  }


/****** profit_evaluate *******************************************************
PROTO	void profit_evaluate(double *par, int m_dat, double *fvec,
		void *data, int *info)
PURPOSE	Provide a function returning residuals to lmfit.
INPUT	Pointer to the vector of parameters,
	number of data points,
	pointer to the vector of residuals (output),
	pointer to the data structure (unused),
	pointer to the info structure (unused).
OUTPUT	-.
NOTES	Input arguments are there only for compatibility purposes (unused)
AUTHOR	E. Bertin (IAP)
VERSION	05/12/2006
 ***/
void	profit_evaluate(double *par, int m_dat, double *fvec,
			void *data, int *info)
  {
  profit_residuals(theprofit, the_field, the_wfield, the_obj,
		par, fvec);

  return;
  }


/****** profit_residuals ******************************************************
PROTO	double *prof_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj, double *param, double *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	pointer to the field,
	pointer to the field weight,
	pointer to the obj,
	pointer to the model parameters (output),
	pointer to the computed residuals (output).
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	08/12/2006
 ***/
double	*profit_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj, double *param, double *resi)
  {
   int		p;

  memset(profit->modpix, 0,
	profit->modnaxisn[0]*profit->modnaxisn[1]*sizeof(double));
  for (p=0; p<profit->nparam; p++)
    profit->param[p] = param[p];
  for (p=0; p<profit->nprof; p++)
    prof_add(profit->prof[p], profit);
  profit_convolve(profit);
  profit_resample(profit);
  profit_compresi(profit, field, wfield, obj, resi);

  return resi;
  }


/****** profit_compresi ******************************************************
PROTO	double *prof_compresi(profitstruct *profit,
			picstruct *field, picstruct *wfield, objstruct *obj,
			double *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	pointer to the field,
	pointer to the field weight,
	pointer to the obj,
	vector of residuals (output).
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	09/12/2006
 ***/
double	*profit_compresi(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj, double *resi)
  {
   double	*resit,
		invsig;
   PIXTYPE	*objpix, *objweight, *lmodpix,
		val;
  int		i;
  
/* Compute vector of residuals */
  resit = resi;
  objpix = profit->objpix;
  objweight = profit->objweight;
  lmodpix = profit->lmodpix;
  invsig = 1.0/profit->sigma;
  for (i=profit->objnaxisn[0]*profit->objnaxisn[1]; i--; lmodpix++)
    if ((val=*(objpix++))>-BIG)
      *(resit++) = (double)(val - *lmodpix)*invsig;

  return resi;
  }


/****** profit_resample ******************************************************
PROTO	PIXTYPE *prof_resample(profitstruct *profit)
PURPOSE	Resample the current full resolution model to image resolution.
INPUT	Profile-fitting structure.
OUTPUT	Resampled pixmap.
AUTHOR	E. Bertin (IAP)
VERSION	08/12/2006
 ***/
PIXTYPE	*profit_resample(profitstruct *profit)
  {
   double	posin[2], posout[2], dnaxisn[2],
		invpixstep, back;
   PIXTYPE	*pixout;
   int		ixcout,iycout, ixcin,iycin,
		d,i;

  ixcout = profit->objnaxisn[0]/2;
  iycout = profit->objnaxisn[1]/2;
  ixcin = profit->modnaxisn[0]/2;
  iycin = profit->modnaxisn[1]/2;
  invpixstep = 1.0/profit->psf->pixstep;
  
/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 1.0;					/* FITS standard */
    dnaxisn[d] = profit->objnaxisn[d] + 0.99999;
    }

/* Remap each pixel */
  pixout = profit->lmodpix;
  back = profit->param[0];
  for (i=profit->objnaxisn[0]*profit->objnaxisn[1]; i--;)
    {
    posin[0] = (posout[0] - ixcout)*invpixstep + ixcin;
    posin[1] = (posout[1] - iycout)*invpixstep + iycin;
    (*(pixout++) = (PIXTYPE)(interpolate_pix(posin, profit->modpix,
		profit->modnaxisn, INTERP_BILINEAR)) + back);
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 1.0;
    }

  return profit->lmodpix;
  }


/****** profit_convolve *******************************************************
PROTO	void profit_convolve(profitstruct *profit)
PURPOSE	Convolve the composite profile model with PSF.
INPUT	Pointer to the profit structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/12/2006
 ***/
void	profit_convolve(profitstruct *profit)
  {
  if (!profit->psfdft)
    profit_makedft(profit);

  fft_conv(profit->modpix, profit->psfdft, profit->modnaxisn);

  return;
  }


/****** profit_makedft *******************************************************
PROTO	void profit_makedft(profitstruct *profit)
PURPOSE	Create the Fourier transform of the descrambled PSF component.
INPUT	Pointer to the profit structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/12/2006
 ***/
void	profit_makedft(profitstruct *profit)
  {
   psfstruct	*psf;
   double      *mask,*maskt, *ppix;
   double       dx,dy, r,r2,rmin,rmin2,rmax,rmax2,rsig,invrsig2;
   int          width,height,npix,offset, psfwidth,psfheight,psfnpix,
                cpwidth, cpheight,hcpwidth,hcpheight, i,j,x,y;

  if (!(psf=profit->psf))
    return;

  psfwidth = psf->masksize[0];
  psfheight = psf->masksize[1];
  psfnpix = psfwidth*psfheight;
  width = profit->modnaxisn[0];
  height = profit->modnaxisn[1];
  npix = width*height;
  QCALLOC(mask, double, npix);

  cpwidth = (width>psfwidth)?psfwidth:width;
  hcpwidth = cpwidth>>1;
  cpwidth = hcpwidth<<1;
  offset = width - cpwidth;
  cpheight = (height>psfheight)?psfheight:height;
  hcpheight = cpheight>>1;
  cpheight = hcpheight<<1;

/* Frame and descramble the PSF data */
  ppix = psf->maskloc + (psfheight/2)*psfwidth + psfwidth/2;
  maskt = mask;
  for (j=hcpheight; j--; ppix+=psfwidth)
    {
    for (i=hcpwidth; i--;)
      *(maskt++) = *(ppix++);      
    ppix -= cpwidth;
    maskt += offset;
    for (i=hcpwidth; i--;)
      *(maskt++) = *(ppix++);      
    }

  ppix = psf->maskloc + ((psfheight/2)-hcpheight)*psfwidth + psfwidth/2;
  maskt += width*(height-cpheight);
  for (j=hcpheight; j--; ppix+=psfwidth)
    {
    for (i=hcpwidth; i--;)
      *(maskt++) = *(ppix++);      
    ppix -= cpwidth;
    maskt += offset;
    for (i=hcpwidth; i--;)
      *(maskt++) = *(ppix++);      
    }

/* Truncate to a disk that has diameter = (box width) */
  rmax = cpwidth - 1.0 - hcpwidth;
  if (rmax > (r=hcpwidth))
    rmax = r;
  if (rmax > (r=cpheight-1.0-hcpheight))
    rmax = r;
  if (rmax > (r=hcpheight))
    rmax = r;
  if (rmax<1.0)
    rmax = 1.0;
  rmax2 = rmax*rmax;
  rsig = psf->fwhm/psf->pixstep;
  invrsig2 = 1/(2*rsig*rsig);
  rmin = rmax - (3*rsig);     /* 3 sigma annulus (almost no aliasing) */
  rmin2 = rmin*rmin;

  maskt = mask;
  dy = 0.0;
  for (y=hcpheight; y--; dy+=1.0)
    {
    dx = 0.0;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
    dx = -hcpwidth;
    maskt += offset;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
    }
  dy = -hcpheight;
  maskt += width*(height-cpheight);
  for (y=hcpheight; y--; dy+=1.0)
    {
    dx = 0.0;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
    dx = -hcpwidth;
    maskt += offset;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
    }

/* Finally move to Fourier space */
  profit->psfdft = fft_rtf(mask, profit->modnaxisn);

  free(mask);

  return;
  }


/****** profit_copyobjpix *****************************************************
PROTO	int profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy)
PURPOSE	Copy a piece of the input field image to a profit structure.
INPUT	Pointer to the profit structure,
	Pointer to the field structure,
	integer position in X,
	integer position in Y.
OUTPUT	The number of valid pixels copied.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
int	profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy)
  {
   PIXTYPE	*pixin, *pixout;
   int		i,x,y, xmin,xmax,ymin,ymax, w,h,w2,dw, npix;

/* First put the image background to -BIG */
  pixout = profit->objpix;
  w = profit->objnaxisn[0];
  h = profit->objnaxisn[1];
  for (i=w*h; i--;)
    *(pixout++) = -BIG;

/* Don't go further if out of frame!! */
  if (ix<0 || ix>=field->width || iy<field->ymin || iy>=field->ymax)
    return 0;

/* Set the image boundaries */
  pixout = profit->objpix;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<field->ymin)
    {
    pixout += (field->ymin-ymin)*w;
    ymin = field->ymin;
    }
  if (ymax>field->ymax)
    ymax = field->ymax;

  xmin = ix-w/2;
  xmax = xmin + w;
  w2 = w;
  if (xmax>field->width)
    {
    w2 -= xmax-field->width;
    xmax = field->width;
    }
  if (xmin<0)
    {
    pixout -= xmin;
    w2 += xmin;
    xmin = 0;
    }

/* Copy the right pixels to the destination */
  dw = w - w2;
  npix = 0;
  for (y=ymin; y<ymax; y++, pixout+=dw)
    {
    pixin = &PIX(field, xmin, y);
    for (x=w2; x--;)
      if ((*(pixout++) = *(pixin++))>-BIG)
        npix++;
    }

  return npix;
  }


/****** prof_init *************************************************************
PROTO	profstruct prof_init(profitstruct *profit, proftypenum profcode)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	Pointer to the profile-fittinh structure,
	profile type.
OUTPUT	A pointer to an allocated prof structure.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
profstruct	*prof_init(profitstruct *profit, proftypenum profcode)
  {
   profstruct	*prof;
   double	*pix,
		rmax2, rh, re2, dy2,r2, scale, zero, k,n, hinvn;
   int		width,height, ixc,iyc, ix,iy, nsub,
		d,s;

  QCALLOC(prof, profstruct, 1);
  prof->code = profcode;
  switch(profcode)
    {
    case PROF_SERSIC:
      prof->naxis = 3;
      width =  prof->naxisn[0] = PROFIT_PROFRES;
      height = prof->naxisn[1] = PROFIT_PROFRES;
      nsub = prof->naxisn[2] = PROFIT_PROFSRES;
      QCALLOC(prof->pix, double, width*height*nsub);
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      re2 = prof->typscale = width/64.0;
      re2 *= re2;
      zero = prof->extrazero[0] = 0.2;
      scale = prof->extrascale[0]= 8.0/PROFIT_PROFSRES;
      pix = prof->pix;
      for (s=0; s<nsub; s++)
        {
        n = s*scale + zero;
        hinvn = 0.5/n;
        k = -1.0/3.0 + 2.0*n + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
		+ 131.0/(1148175*n*n*n);
        for (iy=0; iy<height; iy++)
          {
          dy2 = (iy-iyc)*(iy-iyc);
          for (ix=0; ix<width; ix++)
            {
            r2 = dy2 + (ix-ixc)*(ix-ixc);
            *(pix++) = (r2<rmax2)? exp(-k*pow(r2/re2,hinvn)) : 0.0;
            }
          }
        }
      prof->x[0] = &profit->param[1];
      prof->x[1] = &profit->param[2];
      prof->amp = &profit->param[7];
      prof->scale[0] = &profit->param[8];
      prof->scale[1] = &profit->param[9];
      prof->posangle = &profit->param[10];
      prof->extra[0] = &profit->param[11];
      break;
    case PROF_DEVAUCOULEURS:
      prof->naxis = 2;
      width =  prof->naxisn[0] = PROFIT_PROFRES;
      height = prof->naxisn[1] = PROFIT_PROFRES;
      QCALLOC(prof->pix, double, width*height);
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      rh = prof->typscale = width/128.0;
      rh *= rh;
      pix = prof->pix;
      for (iy=0; iy<height; iy++)
        {
        dy2 = (iy-iyc)*(iy-iyc);
        for (ix=0; ix<width; ix++)
          {
          r2 = dy2 + (ix-ixc)*(ix-ixc);
          *(pix++) = (r2<rmax2)? exp(-pow(r2/rh,0.125)) : 0.0;
          }
        }
      prof->x[0] = &profit->param[1];
      prof->x[1] = &profit->param[2];
      prof->amp = &profit->param[7];
      prof->scale[0] = &profit->param[8];
      prof->scale[1] = &profit->param[9];
      prof->posangle = &profit->param[10];
      break;
    case PROF_EXPONENTIAL:
      prof->naxis = 2;
      width =  prof->naxisn[0] = PROFIT_PROFRES;
      height = prof->naxisn[1] = PROFIT_PROFRES;
      QCALLOC(prof->pix, double, width*height);
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      rh = prof->typscale = width/128.0;
      pix = prof->pix;
      for (iy=0; iy<height; iy++)
        {
        dy2 = (iy-iyc)*(iy-iyc);
        for (ix=0; ix<width; ix++)
          {
          r2 = dy2 + (ix-ixc)*(ix-ixc);
          *(pix++) = (r2<rmax2)? exp(-sqrt(r2)/rh) : 0.0;
          }
        }
      prof->x[0] = &profit->param[1];
      prof->x[1] = &profit->param[2];
      prof->amp = &profit->param[3];
      prof->scale[0] = &profit->param[4];
      prof->scale[1] = &profit->param[5];
      prof->posangle = &profit->param[6];
      break;
    default:
      break;
    }

  prof->kernelnlines = 1;
  for (d=0; d<prof->naxis; d++)
    {
    prof->interptype[d] = INTERP_BILINEAR;
    prof->kernelnlines *= 
	(prof->kernelwidth[d] = interp_kernwidth[prof->interptype[d]]);
    }
  prof->kernelnlines /= prof->kernelwidth[0];
  QMALLOC(prof->kernelbuf, double, prof->kernelnlines);

  return prof;
  }  


/****** prof_end **************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile structure.
INPUT	Prof structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
void	prof_end(profstruct *prof)
  {
  free(prof->pix);
  free(prof->kernelbuf);

  free(prof);

  return;
  }


/****** prof_add **************************************************************
PROTO	void prof_add(profstruct *prof, profitstruct *profit)
PURPOSE	Add a model profile to an image.
INPUT	Profile structure,
	profile-fitting structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
   double	posin[PROFIT_MAXEXTRA], posout[2], dnaxisn[2],
		*pixout,
		amp,ctheta,stheta,cd11,cd12,cd21,cd22, dx1,dx2, x1,x2,
		x1cin,x2cin, x1cout,x2cout, xscale,yscale;
   int		npix,
		d,e,i;

/* Initialize multi-dimensional counters */
  for (d=0; d<prof->naxis; d++)
    {
    posout[d] = 1.0;	/* FITS standard */
    dnaxisn[d] = profit->modnaxisn[d] + 0.99999;
    }

  for (e=0; e<prof->naxis - 2; e++)
    {
    d = 2+e;
/*-- Compute position along axis */
    posin[d] = (*prof->extra[e]-prof->extrazero[e])/prof->extrascale[e]+1.0;
/*-- Keep position within boundaries and let interpolation do the rest */
    if (posin[d] < 0.99999)
      {
      if (prof->extracycleflag[e])
        posin[d] += (double)prof->naxisn[d];
      else
        posin[d] = 1.0;
      }
    else if (posin[d] > (double)prof->naxisn[d])
      {
      if (prof->extracycleflag[e])
        posin[d] = (prof->extracycleflag[e])?
		  fmod(posin[d], (double)prof->naxisn[d])
		: (double)prof->naxisn[d];
      }
    }

  amp = *prof->amp;
/* Compute Profile CD matrix */
  ctheta = cos(*prof->posangle*DEG);
  stheta = sin(*prof->posangle*DEG);
  xscale = (*prof->scale[0] == 0.0)? 0.0 : 1.0 / *prof->scale[0];
  yscale = (*prof->scale[1] == 0.0)? 0.0 : 1.0 / *prof->scale[1];
  cd11 = xscale*ctheta;
  cd12 = xscale*stheta;
  cd21 =-yscale*stheta;
  cd22 = yscale*ctheta;

  x1cout = (double)(profit->modnaxisn[0]/2);
  x2cout = (double)(profit->modnaxisn[1]/2);
  x1cin = (double)(prof->naxisn[0]/2);
  x2cin = (double)(prof->naxisn[1]/2);

  dx1 = *prof->x[0];
  dx2 = *prof->x[1];

/* Remap each pixel */
  pixout = profit->modpix;
  npix = profit->modnaxisn[0]*profit->modnaxisn[1];
  for (i=npix; i--;)
    {
    x1 = posout[0] - x1cout - dx1;
    x2 = posout[1] - x2cout - dx2;
    posin[0] = cd11*x1 + cd12*x2 + x1cin;
    posin[1] = cd21*x1 + cd22*x2 + x2cin;
    *(pixout++) += amp*prof_interpolate(prof, posin);
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 1.0;
    }

  return;
  }


/****** prof_interpolate ******************************************************
PROTO	double	prof_interpolate(profstruct *prof, double *posin)
PURPOSE	Interpolate a multidimensional model profile at a given position.
INPUT	Profile structure,
	input position vector.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
static double	prof_interpolate(profstruct *prof, double *posin)
  {
   double		dpos[2+PROFIT_MAXEXTRA],
			kernel_vector[INTERP_MAXKERNELWIDTH],
			*kvector, *pixin,*pixout,
			val;
   long			step[2+PROFIT_MAXEXTRA],
			start, fac;
   int			linecount[2+PROFIT_MAXEXTRA],
			*naxisn,
			i,j,n, ival, nlines, kwidth,width, badpixflag, naxis;

  naxis = prof->naxis;
  naxisn = prof->naxisn;
  start = 0;
  fac = 1;
  for (n=0; n<naxis; n++)
    {
    val = *(posin++);
    width = *(naxisn++);
/*-- Get the integer part of the current coordinate or nearest neighbour */
    ival = (prof->interptype[n]==INTERP_NEARESTNEIGHBOUR)?
                                        (int)(val-0.50001):(int)val;
/*-- Store the fractional part of the current coordinate */
    dpos[n] = val - ival;
/*-- Check if interpolation start/end exceed image boundary... */
    kwidth = prof->kernelwidth[n];
    ival-=kwidth/2;
    if (ival<0 || ival+kwidth<=0 || ival+kwidth>width)
      return 0.0;

/*-- Update starting pointer */
    start += ival*fac;
/*-- Update step between interpolated regions */
    step[n] = fac*(width-kwidth);
    linecount[n] = 0.0;
    fac *= width;
    }

/* Update Interpolation kernel vectors */
  make_kernel(*dpos, kernel_vector, prof->interptype[0]);
  kwidth = prof->kernelwidth[0];
  nlines = prof->kernelnlines;
/* First step: interpolate along NAXIS1 from the data themselves */
  badpixflag = 0;
  pixin = prof->pix+start;
  pixout = prof->kernelbuf;
  for (j=nlines; j--;)
    {
    val = 0.0;
    kvector = kernel_vector;
    for (i=kwidth; i--;)
       val += *(kvector++)**(pixin++);
    *(pixout++) = val;
    for (n=1; n<naxis; n++)
      {
      pixin+=step[n-1];
      if (++linecount[n]<prof->kernelwidth[n])
        break;
      else
        linecount[n] = 0;       /* No need to initialize it to 0! */
      }
    }

/* Second step: interpolate along other axes from the interpolation buffer */
  for (n=1; n<naxis; n++)
    {
    make_kernel(dpos[n], kernel_vector, prof->interptype[n]);
    kwidth = prof->kernelwidth[n];
    pixout = pixin = prof->kernelbuf;
    for (j = (nlines/=kwidth); j--;)
      {
      val = 0.0;
      kvector = kernel_vector;
      for (i=kwidth; i--;)
        val += *(kvector++)**(pixin++);
      *(pixout++) = val;
     }
    }

  return prof->kernelbuf[0];
  }


/****** interpolate_pix ******************************************************
PROTO	void interpolate_pix(double *posin, double *pix, int naxisn,
		interpenum interptype)
PURPOSE	Interpolate a model profile at a given position.
INPUT	Profile structure,
	input position vector,
	input pixmap dimension vector,
	interpolation type.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2006
 ***/
static double	interpolate_pix(double *posin, double *pix, int *naxisn,
		interpenum interptype)
  {
   double	buffer[INTERP_MAXKERNELWIDTH],
		kernel[INTERP_MAXKERNELWIDTH], dpos[2],
		*kvector, *pixin, *pixout,
		val;
   int		fac, ival, kwidth, start, width, step,
		i,j, n;

  kwidth = interp_kernwidth[interptype];
  start = 0;
  fac = 1;
  for (n=0; n<2; n++)
    {
    val = *(posin++);
    width = naxisn[n];
/*-- Get the integer part of the current coordinate or nearest neighbour */
    ival = (interptype==INTERP_NEARESTNEIGHBOUR)?
                                        (int)(val-0.50001):(int)val;
/*-- Store the fractional part of the current coordinate */
    dpos[n] = val - ival;
/*-- Check if interpolation start/end exceed image boundary... */
    ival-=kwidth/2;
    if (ival<0 || ival+kwidth<=0 || ival+kwidth>width)
      return 0.0;
/*-- Update starting pointer */
    start += ival*fac;
/*-- Update step between interpolated regions */
    fac *= width;
    }

/* First step: interpolate along NAXIS1 from the data themselves */
  make_kernel(dpos[0], kernel, interptype);
  step = naxisn[0]-kwidth;
  pixin = pix+start;
  pixout = buffer;
  for (j=kwidth; j--;)
    {
    val = 0.0;
    kvector = kernel;
    for (i=kwidth; i--;)
      val += *(kvector++)**(pixin++);
    *(pixout++) = val;
    pixin += step;
    }

/* Second step: interpolate along NAXIS2 from the interpolation buffer */
  make_kernel(dpos[1], kernel, interptype);
  pixin = buffer;
  val = 0.0;
  kvector = kernel;
  for (i=kwidth; i--;)
    val += *(kvector++)**(pixin++);

  return val;
  }


/****** make_kernel **********************************************************
PROTO   void make_kernel(double pos, double *kernel, interpenum interptype)
PURPOSE Conpute interpolation-kernel data
INPUT   Position,
        Pointer to the output kernel data,
        Interpolation method.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 18/04/2003
 ***/
static void	make_kernel(double pos, double *kernel, interpenum interptype)
  {
   double       x, val, sinx1,sinx2,sinx3,sinx4;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR)
    {
    *(kernel++) = 1.0-pos;
    *kernel = pos;
    }
  else if (interptype == INTERP_LANCZOS2)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/2.0*(pos+1.0);
      val = (*(kernel++) = (sinx1=sin(x))/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -(sinx2=sin(x))/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0;
      val += (*kernel = sinx2/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS3)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/3.0*(pos+2.0);
      val = (*(kernel++) = (sinx1=sin(x))/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx2=-sin(x))/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx3=sin(x))/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/3.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS4)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/4.0*(pos+3.0);
      val = (*(kernel++) = (sinx1=sin(x))/(x*x));
      x += PI/4.0;
      val +=(*(kernel++) = -(sinx2=sin(x))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = (sinx3=sin(x))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -(sinx4=sin(x))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx3/(x*x));
      x += PI/4.0;
      val += (*kernel = sinx4/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else
    error(EXIT_FAILURE, "*Internal Error*: Unknown interpolation type in ",
                "make_kernel()");

  return;
  }

