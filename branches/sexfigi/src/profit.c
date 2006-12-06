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
*	Last modify:	06/12/2006
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
#include	"check.h"
#include	"image.h"
#include	"profit.h"
#include	"psf.h"

/*------------------------------- variables ---------------------------------*/

int		interp_kernwidth[5]={1,2,4,6,8};

/* "Local" global variables; it seems dirty but it simplifies a lot */
/* interfacing to the LM routines */
static objstruct	*the_obj;
static picstruct	*the_field, *the_wfield;
static profitstruct	*profit;

/****** profit_init ***********************************************************
PROTO	profitstruct profit_init(psfstruct *psf,
				proftypestruct *profcode, int nprof)
PURPOSE	Allocate and initialize a new profile-fitting structure.
INPUT	-.
OUTPUT	A pointer to an allocated profit structure.
AUTHOR	E. Bertin (IAP)
VERSION	06/12/2006
 ***/
profitstruct	*profit_init(psfstruct *psf,
				proftypestruct *profcode, int nprof)
  {
   profitstruct		*profit;
   int			p;

  QCALLOC(profit, profitstruct, 1);
  profit->nparam = 6;
  QMALLOC(profit->param, double, profit->nparam);
  QMALLOC(profit->parammin, double, profit->nparam);
  QMALLOC(profit->parammax, double, profit->nparam);
  QMALLOC(profit->prof, profstruct *, nprof);
  profit->nprof = nprof;
  for (p=0; p<nprof; p++)
    profit->prof[p] = prof_init(profcode[p]);
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
  free(profit->psfdft);
  free(profit->param);
  free(profit->parammin);
  free(profit->parammax);
  free(profit->objpix);
  free(profit->resi);
  free(profit);

  return;
  }


/****** prof_fit ************************************************************
PROTO	fitstruct *prof_fit(profstruct *prof, int nprof, psfstruct *psf,
		picstruct *field, picstruct *wfield, objstruct *obj)
PURPOSE	Fit profile(s) convolved with the PSF to a detected object.
INPUT	Array of profile structures,
	Number of profiles,
	Pointer to a PSF structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Pointer to an allocated fit structure (containing details about the
	fit).
NOTES	It is a modified version of the lm_minimize() of lmfit.
AUTHOR	E. Bertin (IAP)
VERSION	10/11/2006
 ***/
fitstruct	*prof_fit(profstruct *prof, int nprof, psfstruct *psf,
		picstruct *field, picstruct *wfield, objstruct *obj)
  {
    profstruct		*prof;
    profitstruct	*profit;
    double		*diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
    int			*ipvt;

/* Use (dirty) global variables to interface with lmfit */
  the_profit = profit = profit_init(psf, EXPONENTIAL, 1);
  profit->objnaxisn[0] = (obj->xmax - obj->xmin + 1) * 2;
  profit->objnaxisn[1] = (obj->ymax - obj->ymin + 1) * 2;
  QMALLOC(profit->objpix, double, profit->objnaxisn[0]*profit->objnaxisn[1]);
  prof->nresi = profit_copyimage(profit, field,
			(int)(obj->mx+0.49999),(int)(obj->my+0.49999));
  if (prof->nresi >= prof->nparam)
    {
    QMALLOC(prof->resi, double, prof->nresi);
    }
  else
    {
    profit_end(profit);
    return NULL;
    }

  the_field = field;
  the_wfield = wfield;
  the_obj = obj;

/* Allocate work space */
  n = profit->nparam;
  m = profit->nresi;

  lm_control_type	control;
  lm_initialize_control(&control);

  QMALLOC(diag, double,n);
  QMALLOC(qtf, double, n);
  QMALLOC(fjac, double,n*m);
  QMALLOC(wa1, double, n);
  QMALLOC(wa2, double, n);
  QMALLOC(wa3, double, n);
  QMALLOC(wa4, double,   m);
  QMALLOC(ipvt, int,   n);

/* Perform fit */

  control->info = 0;
  control->nfev = 0;

/* This goes through the modified legacy interface */
  lm_lmdif( m, n, profit->param, profit->resi,
	control->ftol, control->xtol, control->gtol,
	control->maxcall*(n+1), control->epsilon, diag, 1,
	control->stepbound, &(control->info),
	&(control->nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
	profit_evaluate, profit_printout, NULL);

  control->fnorm = lm_enorm(m, profit->resi);
  if (control->info < 0 )
    control->info = 10;

/* clean up. */
  profit_end(profit);
  free(diag);
  free(qtf); 
  free(fjac);
  free(wa1); 
  free(wa2); 
  free(wa3 );
  free(wa4); 
  free(ipvt);

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
VERSION	05/12/2006
 ***/
void	profit_evaluate(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev )
  {

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
  profit_residuals(the_profit, the_field, the_wfield, the_obj);

  return;
  }


/****** profit_residuals ******************************************************
PROTO	double *prof_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	15/11/2006
 ***/
double	*profit_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj)
  {
   int		p;

  memset(profit->pix, 0, profit->naxisn[0]*profit->naxisn[1]*sizeof(double));
  for (p=0; p<profit->nprof; p++)
    prof_add(profit->prof[p], profit);
  prof_convolve(profit);
  prof_compresi(profit, field, wfield, obj);

  return profit->resi;
  }


/****** profit_compresi ******************************************************
PROTO	double *prof_compresi(profitstruct *profit,
			picstruct *field, picstruct *wfield, objstruct *obj)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	30/11/2006
 ***/
double	*profit_compresi(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj)
  {
   double	posout[2], dnaxisn[2],
		*pixout, *weighout,
		invpixstep, back;
   int		ixcout,iycout, ixcin,iycin,
		d,i;
  
  ixcout = profit->naxisn[0]/2;
  iycout = profit->naxisn[1]/2;
  ixcin = profit->modnaxisn[0]/2;
  iycin = profit->modnaxisn[1]/2;
  invpixstep = 1.0/profit->psf->pixstep;
  
/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 0.0;
    dnaxisn[d] = profit->objnaxisn[d] - 0.00001;
    }

/* Remap each pixel */
  pixout = profit->objpix;
  weightout = profit->objweight;
  back = profit->param[0];
  for (i=npix; i--;)
    {
    posin[0] = (posout[0] - ixcout)*invpixstep + ixcin;
    posin[1] = (posout[1] - iycout)*invpixstep + iycin;
    *(resi++) = (*(pixout++) - interpolate_pix(posin, prof->modpix,
		profit->modnaxisn, INTERP_LANCZOS3)) - back;
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 0.0;
    }

  return profit->resi;
  }


/****** profit_convolve *******************************************************
PROTO	void profit_convolve(profitstruct *profit)
PURPOSE	Convolve the composite profile model with PSF.
INPUT	Pointer to the profit structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/11/2006
 ***/
void	profit_convolve(profitstruct *profit)
  {
   double	*psfbuf;
   int		p;

  psf = profit->psf;
  if (!prof->psfdft)
    {
    QCALLOC(psfbuf, double, profit->naxisn[0]*profit->naxisn[1]);
    psfw = psf->masksize[0];
    idx = (idw = profit->naxisn[0] - psfw)/2;
/*-- Recenter slightly if PSF size odd and profile size even */
    if (psfw&2 && !(profit->naxisn[0]&2))
      idx++;
    psfh = psf->masksize[1];
    idy = (profit->naxisn[1] - psfh)/2;
/*-- Recenter slightly if PSF size odd and profile size even */
    if (psfh&2 && !(profit->naxisn[1]&2))
      idy++;
    psfbuft = psfbuf + idy*profit->naxisn[1] + idx;
    for (iy=psfh; iy--; psfbuft += idw)
      for (ix=psfw; ix--;)
        *(psfbuft++) = *(psfin++);
    profit->psfdft = fft_rtf(psfbuf, profit->naxisn);
    free(psfbuf);
    }

  fft_conv(profit->pix, profit->psfdft, profit->naxisn)

  return;
  }


/****** profit_copyobjpix *****************************************************
PROTO	int profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy, int w, int h)
PURPOSE	Copy a piece of the input field image to a profit structure.
INPUT	Pointer to the profit structure,
	Pointer to the field structure.
OUTPUT	The number of valid pixels copied.
AUTHOR	E. Bertin (IAP)
VERSION	12/06/2006
 ***/
int	profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy, int w, int h)
  {
   double	*pixout;
   PIXTYPE	*pixin;
   int		i,x,y, xmin,xmax,ymin,ymax, w2, npix;

/* First put the image background to -BIG */
  pixout = profit->objpix;
  for (i=w*h; i--;)
    *(pixout++) = -BIG;

/* Don't go further if out of frame!! */
  if (ix<0 || ix>=field->width || iy<field->ymin || iy>=field->ymax)
    return 0;

/* Set the image boundaries */
  w2 = w;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<field->ymin)
    {
    dest += (field->ymin-ymin)*w;
    ymin = field->ymin;
    }
  if (ymax>field->ymax)
    ymax = field->ymax;

  xmin = ix-w/2;
  xmax = xmin + w;
  if (xmax>field->width)
    {
    w2 -= xmax-field->width;
    xmax = field->width;
    }
  if (xmin<0)
    {
    dest += -xmin;
    w2 -= -xmin;
    xmin = 0;
    }

/* Copy the right pixels to the destination */
  pixout = profit->objpix;
  npix = 0;
  for (y=ymin; y<ymax; y++, dest += w)
    {
    pixin = &PIX(field, xmin, y);
    for (x=w2; x--;)
       if ((*(pixout++) = (double)*(pixin++))>-BIG)
         npix++;
    }

  return npix;
  }


/****** prof_init *************************************************************
PROTO	profstruct prof_init(void)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	-.
OUTPUT	A pointer to an allocated prof structure.
AUTHOR	E. Bertin (IAP)
VERSION	06/12/2006
 ***/
profstruct	*prof_init(proftypestruct *profcode, profitstruct *profit)
  {
   profstruct	*prof;
   double	rmax2, rh, dy2,r2;
   int		width,height, ixc,iyc, ix,iy;

  QCALLOC(prof, profstruct, 1);
  prof->code = profcode;
  switch(profcode)
    {
    case PROF_SERSIC:
      break;
    case PROF_DEVAUCOULEURS:
      break;
    case PROF_EXPONENTIAL:
      prof->naxis = 2;
      width =  prof->naxisn[0] = PROFIT_RESOL;
      height = prof->naxisn[1] = PROFIT_RESOL;
      QCALLOC(prof->pix, double, width*height);
      prof->comppix = prof->pix;
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      rh = width/10.0;
      for (iy=0; iy<height; iy++)
        {
        dy2 = (iy-iyc)*(iy-iyc);
        for (ix=0; ix<width; ix++)
          {
          r2 = dy2 + (ix-ixc)*(ix-ixc);
          if (r2<rmax2)
            *(pixt++) = exp(-sqrt(r2)/rh);
          }
        }
      prof->x[0] = &profit->param[1];
      prof->x[1] = &profit->param[2];
      prof->scale[0] = &profit->param[3];
      prof->scale[1] = &profit->param[4];
      prof->posangle = &profit->param[5];
      break;
    default:
      break;
    }
  

  return prof;
  }  


/****** prof_end **************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile structure.
INPUT	Prof structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/11/2006
 ***/
void	prof_end(profstruct *prof)
  {
  free(prof->pix);
  if (prof->comppix != prof->pix)
    free(prof->comppix);
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
VERSION	06/12/2006
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
   double	posin[2], posout[2], dnaxisn[2],
		*pixout,
		amp,ctheta,stheta,cd11,cd12,cd21,cd22, dx1,dx2, x1,x2;
   int		npix,nextra,
		d,i;

  nextra = prof->naxis-2;
  if (nextra>0)
/*-- Create an interpolated image of model */
    prof_interpextra(prof);

  amp = prof->amp;
/* Compute Profile CD matrix */
  ctheta = cos(*prof->posangle*DEG);
  stheta = sin(*prof->posangle*DEG);
  cd11 = *prof->scale[0]*ctheta;
  cd12 =-*prof->scale[1]*stheta;
  cd21 = *prof->scale[0]*stheta;
  cd22 = *prof->scale[1]*ctheta;
  dx1 = prof->x[0];
  dx2 = prof->x[1];

/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 0.0;
    dnaxisn[d] = prof->naxisn[d] - 0.00001;
    }

  pixout = profit->modpix;
  npix = profit->modnaxisn[0]*prof->modnaxisn[1];
/* Remap each pixel */
  for (i=npix; i--;)
    {
    x1 = posout[0] - dx1;
    x2 = posout[1] - dx2;
    posin[0] = cd11*x1 + cd12*x2;
    posin[1] = cd21*x1 + cd22*x2;
    *(pixout++) += amp*interpolate_pix(posin, prof->comppix, prof->naxisn,
		prof->interptype);
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 0.0;
    }

  return;
  }


/****** interpolate_pix ******************************************************
PROTO	void interpolate_pix(double *posin, double *pix, int naxisn,
		interpernum interptype)
PURPOSE	Interpolate a model profile at a given position.
INPUT	Profile structure,
	Position vector.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/11/2006
 ***/
static double	interpolate_pix(double *posin, double *pix, int naxisn,
		interpernum interptype)
  {
   double	buffer[INTERP_MAXKERNELWIDTH],
		kernel[INTERP_MAXKERNELWIDTH], dpos[2],
		*kvector, *pixin, *pixout,
		val;
   int		fac, ival, kwidth, start, width,
		i,j, n;

  kwidth = interp_kernwidth[interptype]);
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
  make_kernel(dpos[0], kernel, prof->interptype);
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


/****** prof_interpextra ******************************************************
PROTO	void prof_interpextra(profstruct *prof)
PURPOSE	Interpolate the "extra" components of a model profile.
INPUT	Profile structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/12/2006
 ***/
void	prof_interpextra(profstruct *prof)
  {
   double	profrac[PROFIT_MAXEXTRA],
		*pixpoint, *pixint,
		pos,fac;
   int		pixstep[PROFIT_MAXEXTRA], count[PROFIT_MAXEXTRA],
		npix,nextra,inc,extralen,proint,
		e,i,n;

  npix = prof->naxisn[0]*prof->naxisn[1];
  nextra = prof->naxis-2;
/* Interpolate pixmap over all extra axes */
  inc = npix;
  pixpoint = prof->pix;
  for (e=0; e<nextra; e++)
    {
    extralen = prof->naxisn[e+2];
/*-- Compute position along axis */
    pos = prof->extrapos[e]*prof->extrascale[e]+prof->extraoffset[e];
    proint = (int)pos;
/*-- Keep position within boundaries and let interpolation do the rest */
    if (prof->extracycleflag)
      {
      proint = proint % extralen;
      if (proint < 0)
        proint += extralen;
      }
    else
      {
      if (proint < 0)
        proint = 0;
      else if (proint >= extralen)
        proint = extralen - 1;
      }
    profrac[e] = pos - proint;
    pixpoint += proint*inc;
    pixstep[e] = inc;
    inc *= extralen;
    count[e] = 0;
    }
  nprof = 2 << nextra;	/* Nb of pixmaps involved in interpolation */
  for (n=0; n<nprof; n++)
    {
    fac = 0.0;
    pixint = pixpoint;
    for (e=0; e<nextra; e++)
      {
      fac *= count[e]? 1.0 - profac[e]: profac[e];
      pixint += count[e]*pixstep[e];
      }
/*-- The actual multi-linear interpolation */
    pixout = prof->comppix;
    for (i=npix; i--;)
      *(pixout++) += fac**(pixint++);
    for (e=0; e<nextra; e++)
      if ((++count[e])<2)
        break;
      else
        count[e] = 0;
    }

  return;
  }


