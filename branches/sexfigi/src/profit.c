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
*	Last modify:	18/07/2007
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
#include	"levmar/lm.h"
#include	"lmfit/lmmin.h"
#include	"fft.h"
#include	"check.h"
#include	"psf.h"
#include	"profit.h"


static double	gammln(double);
static double	prof_interpolate(profstruct *prof, double *posin);
static double	interpolate_pix(double *posin, double *pix, int *naxisn,
		interpenum interptype);

static void	make_kernel(double pos, double *kernel, interpenum interptype);

/*------------------------------- variables ---------------------------------*/

int		interp_kernwidth[5]={1,2,4,6,8};
int theniter, the_x,the_y, the_gal;
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
VERSION	30/03/2007
 ***/
profitstruct	*profit_init(psfstruct *psf, proftypenum *profcode, int nprof)
  {
   profitstruct		*profit;
   int			p;

  QCALLOC(profit, profitstruct, 1);
  profit->psf = psf;
  profit->psfdft = NULL;

  profit->nprof = nprof;
  profit->nparam = 0;
  QMALLOC(profit->prof, profstruct *, nprof);

  for (p=0; p<nprof; p++)
    profit->prof[p] = prof_init(profit, profcode[p]);

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
  free(profit);

  return;
  }


/****** profit_fit ************************************************************
PROTO	void profit_fit(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj, obj2struct *obj2)
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
VERSION	18/07/2007
 ***/
void	profit_fit(profitstruct *profit,
		picstruct *field, picstruct *wfield,
		objstruct *obj, obj2struct *obj2)
  {
    psfstruct		*psf;
    checkstruct		*check;
    double		*param, *oldparaminit,
			psf_fwhm, oldchi2;
    int			ix,iy, p, oldniter, flags;


  if (profit->psfdft)
    {
    QFREE(profit->psfdft);
    }

  psf = profit->psf;

/* Compute the local PSF */
  psf_build(psf);

/* Create pixmaps at image resolution */
  psf_fwhm = psf->masksize[0]*psf->pixstep;
  profit->objnaxisn[0] = (((int)((obj->xmax-obj->xmin+1) + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
  profit->objnaxisn[1] = (((int)((obj->ymax-obj->ymin+1) + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
  ix = (int)(obj->mx + 0.49999);	/* internal convention: 1st pix = 0 */
  iy = (int)(obj->my + 0.49999);	/* internal convention: 1st pix = 0 */

/* Use (dirty) global variables to interface with lmfit */
  the_field = field;
  the_wfield = wfield;
  theprofit = profit;
  the_obj = obj;
  the_x = ix;
  the_y = iy;

  QMALLOC(profit->objpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  QMALLOC(profit->lmodpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  profit->nresi = profit_copyobjpix(profit, field, ix,iy);
  if (profit->nresi < profit->nparam)
    {
    for (p=0; p<profit->nparam; p++)
      obj2->prof_vector[p] = 0.0;
    obj2->prof_niter = 0;
    return;
    }

  QCALLOC(profit->resi, double, profit->nresi);

/* Create pixmap at PSF resolution */
  profit->modnaxisn[0] = ((int)(profit->objnaxisn[0]/psf->pixstep +0.4999)/2)*2; 
  profit->modnaxisn[1] = ((int)(profit->objnaxisn[1]/psf->pixstep +0.4999)/2)*2; 
  if (profit->modnaxisn[1] < profit->modnaxisn[0])
    profit->modnaxisn[1] = profit->modnaxisn[0];
  else
    profit->modnaxisn[0] = profit->modnaxisn[1];
/* Allocate memory for the complete model */
  QCALLOC(profit->modpix, double, profit->modnaxisn[0]*profit->modnaxisn[1]);
/* Allocate memory for the partial model */
  QMALLOC(profit->pmodpix, float, profit->modnaxisn[0]*profit->modnaxisn[1]);

/* Set initial guesses and boundaries */
  obj2->prof_flag = 0;
  profit->sigma = obj->sigbkg;
  profit_resetparams(profit, obj, obj2);

for (p=0; p<profit->nparam; p++)
printf("%g ", profit->paraminit[p]);
printf("\n");

the_gal++;
  profit->niter = profit_minimize(profit, PROFIT_MAXITER);
  QMEMCPY(profit->paraminit, oldparaminit, double, profit->nparam);
  if (profit_resetparam(profit, PARAM_ARMS_PITCH, 150.0, 120.0, 170.0)==RETURN_OK)
    {
//    profit_resetparam(profit, PARAM_ARMS_FLUX, obj->peak, 0.0, 1000.0*obj->peak);
//    profit_resetparam(profit, PARAM_ARMS_SCALE, 1.0, 0.5, 10.0);
//    profit_resetparam(profit, PARAM_ARMS_START, 0.5, 0.0,
//			5.0*obj2->hl_radius/1.67835);
//    profit_resetparam(profit, PARAM_ARMS_POSANG, 0.0, -3600.0, 3600.0);
//    profit_resetparam(profit, PARAM_ARMS_WIDTH, 0.7, 0.0, 1.0);
    oldchi2 = profit->chi2;
    oldniter = profit->niter;
    profit_resetparams(profit, obj, obj2);
    profit_resetparam(profit, PARAM_ARMS_PITCH, 150.0, 120.0, 170.0);
    profit->niter = profit_minimize(profit, PROFIT_MAXITER);
    if (profit->chi2 > oldchi2)
      {
      memcpy(profit->paraminit, oldparaminit, profit->nparam*sizeof(double));
      profit->chi2 = oldchi2;
      profit->niter = oldniter;
      }
    else
      obj2->prof_flag |= PROFIT_FLIPPED;
    }  

printf("--> ");
for (p=0; p<profit->nparam; p++)
printf("%g     %g %g\n", profit->paraminit[p], profit->parammin[p], profit->parammax[p]);
printf("(%d)\n", profit->niter);


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

/* Fill measurement parameters */
  if (FLAG(obj2.prof_vector))
    {
    for (p=0; p<profit->nparam; p++)
      obj2->prof_vector[p]= profit->param[p];
    }
  obj2->prof_niter = profit->niter;
  if ((param=profit->paramlist[PARAM_X]))
    obj2->x_prof = ix + 1.0 + *param;		/* FITS convention */
  if ((param=profit->paramlist[PARAM_Y]))
    obj2->y_prof = iy + 1.0 + *param;		/* FITS convention */
  obj2->flux_prof = profit->flux;
  obj2->prof_chi2 = profit->chi2;

/* clean up. */
  free(profit->modpix);
  free(profit->pmodpix);
  free(profit->lmodpix);
  free(profit->objpix);
  free(profit->resi);
  free(oldparaminit);

  return;
  }


/****** profit_minimize *******************************************************
PROTO	void profit_minimize(profitstruct *profit)
PURPOSE	Provide a function returning residuals to lmfit.
INPUT	Pointer to the profit structure involved in the fit,
	maximum number of iterations.
OUTPUT	Number of iterations used.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/06/2007
 ***/
int	profit_minimize(profitstruct *profit, int niter)
  {
   lm_control_type	control;
   double		lm_opts[5],lm_info[LM_INFO_SZ],
			*diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
   int			*ipvt,
			m,n;

/* Allocate work space */
  n = profit->nparam;
  m = profit->nresi;

  profit_boundtounbound(profit, profit->paraminit);

/* Perform fit */
  if (1)
    {
    lm_opts[0] = 1.0e-6;
    lm_opts[1] = 1.0e-12;
    lm_opts[2] = 1.0e-12;
    lm_opts[3] = 1.0e-12;
    lm_opts[4] = 1.0e-6;

    niter = dlevmar_dif(profit_evaluate2, profit->paraminit, profit->resi,
	n, m, niter,  lm_opts, NULL, NULL, NULL, profit);
    }
  else
    {
    control.ftol =      1.0e-12;
    control.xtol =      1.0e-12;
    control.gtol =      1.0e-12;
    control.maxcall =   niter;
    control.epsilon =   1.0e-6;
    control.stepbound = 100.0;
    control.info = 0;
    control.nfev = 0;
//  control.fnorm = lm_enorm(m, profit->resi);
//  if (control.info < 0 )
//    control.info = 10;
    QMALLOC(diag, double,n);
    QMALLOC(qtf, double, n);
    QMALLOC(fjac, double,n*m);
    QMALLOC(wa1, double, n);
    QMALLOC(wa2, double, n);
    QMALLOC(wa3, double, n);
    QMALLOC(wa4, double,   m);
    QMALLOC(ipvt, int,   n);
    lm_lmdif(profit->nresi, profit->nparam, profit->paraminit, profit->resi,
	control.ftol, control.xtol, control.gtol,
	control.maxcall*(n+1), control.epsilon, diag, 1,
	control.stepbound, &(control.info),
	&(control.nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
	profit_evaluate, profit_printout, profit);
    niter = control.nfev;
/*-- clean up. */
    free(diag);
    free(qtf); 
    free(fjac);
    free(wa1); 
    free(wa2); 
    free(wa3 );
    free(wa4); 
    free(ipvt);
    }

  profit_unboundtobound(profit, profit->paraminit);


  return niter;
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
   checkstruct	*check;
   profitstruct	*profit;
   char		filename[256];
   static int	itero;

  profit = (profitstruct *)data;

  if (0 && (iter!=itero || iter<0))
    {
    if (iter<0)
      itero++;
    else
      itero = iter;
    sprintf(filename, "check_%d_%04d.fits", the_gal, itero);
    check=initcheck(filename, CHECK_PROFILES, 0);
    reinitcheck(the_field, check);
    addcheck(check, profit->lmodpix, profit->objnaxisn[0],profit->objnaxisn[1],
		the_x,the_y, 1.0);

    reendcheck(the_field, check);
    endcheck(check);
//    lm_print_default(n_par, par, m_dat, fvec, data, iflag, iter, nfev);
    }

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
VERSION	18/05/2007
 ***/
void	profit_evaluate(double *par, int m_dat, double *fvec,
			void *data, int *info)
  {
   profitstruct *profit;

  profit = (profitstruct *)data;
  profit_unboundtobound(profit, par);
  profit_residuals(theprofit, the_field, the_wfield, the_obj,
		par, fvec);
  profit_boundtounbound(profit, par);

  return;
  }


/****** profit_evaluate2 ******************************************************
PROTO	void profit_evaluate2(double *par, double *fvec, int m, int n,
		void *adata)
PURPOSE	Provide a function returning residuals to levmar.
INPUT	Pointer to the vector of parameters,
	pointer to the vector of residuals (output),
	number of model parameters,
	number of data points,
	pointer to a data structure (unused).
OUTPUT	-.
NOTES	Input arguments are there only for compatibility purposes (unused)
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2007
 ***/
void	profit_evaluate2(double *par, double *fvec, int m, int n,
			void *adata)
  {
   profitstruct	*profit;

  profit = (profitstruct *)adata;
  profit_unboundtobound(profit, par);
  profit_residuals(profit, the_field, the_wfield, the_obj, par, fvec);
  profit_boundtounbound(profit, par);
  profit_printout(m, par, n, fvec, adata, 0, -1, 0 );
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
VERSION	12/07/2007
 ***/
double	*profit_compresi(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj, double *resi)
  {
   double	*resit,
		invsig, error;
   PIXTYPE	*objpix, *objweight, *lmodpix,
		val,val2;
   int		i, npix;
  
/* Compute vector of residuals */
  npix = profit->objnaxisn[0]*profit->objnaxisn[1];
  resit = resi;
  objpix = profit->objpix;
  objweight = profit->objweight;
  lmodpix = profit->lmodpix;
  invsig = 1.0/(PROFIT_DYNPARAM*profit->sigma);
  error = 0.0;
  for (i=npix; i--; lmodpix++)
    if ((val=*(objpix++))>-BIG)
      {
      val2 = (double)(val - *lmodpix)*invsig;
      val2 = val2>0.0? log(1.0+val2) : -log(1.0-val2);
      *(resit++) = val2;
      error += val2*val2;
      }

  profit->chi2 = PROFIT_DYNPARAM*PROFIT_DYNPARAM*error/npix;

  return resi;
  }


/****** profit_resample ******************************************************
PROTO	PIXTYPE *prof_resample(profitstruct *profit)
PURPOSE	Resample the current full resolution model to image resolution.
INPUT	Profile-fitting structure.
OUTPUT	Resampled pixmap.
AUTHOR	E. Bertin (IAP)
VERSION	22/06/2007
 ***/
PIXTYPE	*profit_resample(profitstruct *profit)
  {
   double	posin[2], posout[2], dnaxisn[2],
		*dx,*dy,
		xcout,ycout, xcin,ycin, invpixstep, flux;
   PIXTYPE	*pixout;
   int		d,i;

  xcout = (double)(profit->objnaxisn[0]/2) + 1.0;	/* FITS convention */
  if ((dx=(profit->paramlist[PARAM_X])))
    xcout += *dx;
  ycout = (double)(profit->objnaxisn[1]/2) + 1.0;	/* FITS convention */
  if ((dy=(profit->paramlist[PARAM_Y])))
    ycout += *dy;
  xcin = (profit->modnaxisn[0]/2) + 1.0;		/* FITS convention */
  ycin = (profit->modnaxisn[1]/2) + 1.0;		/* FITS convention */
  invpixstep = 1.0/profit->psf->pixstep;

/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 1.0;					/* FITS convention */
    dnaxisn[d] = profit->objnaxisn[d]+0.5;
    }

/* Remap each pixel */
  pixout = profit->lmodpix;
  flux = 0.0;
  for (i=profit->objnaxisn[0]*profit->objnaxisn[1]; i--;)
    {
    posin[0] = (posout[0] - xcout)*invpixstep + xcin;
    posin[1] = (posout[1] - ycout)*invpixstep + ycin;
    flux += ((*(pixout++) = (PIXTYPE)(interpolate_pix(posin, profit->modpix,
		profit->modnaxisn, INTERP_LANCZOS3))));
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 1.0;
    }

  profit->flux = flux;

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
	integer position in X (SExtractor convention),
	integer position in Y (SExtractor convention).
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


/****** profit_spiralindex ****************************************************
PROTO	double profit_spiralindex(profitstruct *profit, objstruct *obj,
			obj2struct *obj2)
PURPOSE	Compute the spiral index of a galaxy image (positive for arms
	extending counter-clockwise and negative for arms extending CW, 0 for
	no spiral pattern).
INPUT	Profile-fitting structure,
	pointer to the obj,
	pointer to the obj2.
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	24/06/2007
 ***/
double profit_spiralindex(profitstruct *profit, objstruct *obj,
			obj2struct *obj2)
  {
   double	*dx,*dy, *fdx,*fdy, *gdx,*gdy, *gdxt,*gdyt, *pix,
		fwhm, invtwosigma2, hw,hh, ohw,ohh, x,y,xstart, tx,ty,txstart,
		gx,gy, r2, spirindex, invsig, val, sep;
   PIXTYPE	*fpix;
   int		i,j, npix;

  npix = profit->objnaxisn[0]*profit->objnaxisn[1];

/* Compute simple derivative vectors at a fraction of the object scale */
  fwhm = obj2->hl_radius * 2.0 / 4.0;
  if (fwhm < 2.0)
    fwhm = 2.0;
  sep = 2.0;

  invtwosigma2 = -(2.35*2.35/(2.0*fwhm*fwhm));
  hw = (double)(profit->objnaxisn[0]/2);
  ohw = profit->objnaxisn[0] - hw;
  hh = (double)(profit->objnaxisn[1]/2);
  ohh = profit->objnaxisn[1] - hh;
  txstart = -hw;
  ty = -hh;
  QMALLOC(dx, double, npix);
  pix = dx;
  for (j=profit->objnaxisn[1]; j--; ty+=1.0)
    {
    tx = txstart;
    y = ty < -0.5? ty + hh : ty - ohh;
    for (i=profit->objnaxisn[0]; i--; tx+=1.0)
      {
      x = tx < -0.5? tx + hw : tx - ohw;
      *(pix++) = exp(invtwosigma2*((x+sep)*(x+sep)+y*y))
		- exp(invtwosigma2*((x-sep)*(x-sep)+y*y));
      }
    }
  QMALLOC(dy, double, npix);
  pix = dy;
  ty = -hh;
  for (j=profit->objnaxisn[1]; j--; ty+=1.0)
    {
    tx = txstart;
    y = ty < -0.5? ty + hh : ty - ohh;
    for (i=profit->objnaxisn[0]; i--; tx+=1.0)
      {
      x = tx < -0.5? tx + hw : tx - ohw;
      *(pix++) = exp(invtwosigma2*(x*x+(y+sep)*(y+sep)))
		- exp(invtwosigma2*(x*x+(y-sep)*(y-sep)));
      }
    }

  QMALLOC(gdx, double, npix);
  gdxt = gdx;
  fpix = profit->objpix;
  invsig = npix/profit->sigma;
  for (i=npix; i--; fpix++)
    {
    val = *fpix > -1e29? *fpix*invsig : 0.0;
    *(gdxt++) = (val>0.0? log(1.0+val) : -log(1.0-val));
    }
  QMEMCPY(gdx, gdy, double, npix);
  fdx = fft_rtf(dx, profit->objnaxisn);
  fft_conv(gdx, fdx, profit->objnaxisn);
  fdy = fft_rtf(dy, profit->objnaxisn);
  fft_conv(gdy, fdy, profit->objnaxisn);

/* Compute estimator */
  invtwosigma2 = -1.18*1.18/(2.0*obj2->hl_radius*obj2->hl_radius);
  xstart = -hw - obj->mx + (int)(obj->mx+0.49999);;
  y = -hh -  obj->my + (int)(obj->my+0.49999);;
  spirindex = 0.0;
  gdxt = gdx;
  gdyt = gdy;
  for (j=profit->objnaxisn[1]; j--; y+=1.0)
    {
    x = xstart;
    for (i=profit->objnaxisn[0]; i--; x+=1.0)
      {
      gx = *(gdxt++);
      gy = *(gdyt++);
      if ((r2=x*x+y*y)>0.0)
        spirindex += (x*y*(gx*gx-gy*gy)+gx*gy*(y*y-x*x))/r2
			* exp(invtwosigma2*r2);
      }
    }

  free(dx);
  free(dy);
  free(fdx);
  free(fdy);
  free(gdx);
  free(gdy);

  return spirindex;
  }


/****** profit_addparam *******************************************************
PROTO	void profit_addparam(profitstruct *profit, paramenum paramindex,
		double **param)
PURPOSE	Add a profile parameter to the list of fitted items.
INPUT	Pointer to the profit structure,
	Parameter index,
	Pointer to the parameter pointer.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/03/2007
 ***/
void	profit_addparam(profitstruct *profit, paramenum paramindex,
		double **param)
  {
/* Check whether the parameter has already be registered */
  if (profit->paramlist[paramindex])
/*-- Yes */
    *param = profit->paramlist[paramindex];
  else
/*-- No */
    *param = profit->paramlist[paramindex] = &profit->param[profit->nparam++];

  return;
  }


/****** profit_resetparams ****************************************************
PROTO	void profit_resetparams(profitstruct *profit, objstruct *obj,
		*obj2struct *obj)
PURPOSE	Set the initial, lower and upper boundary values of profile parameters.
INPUT	Pointer to the profit structure,
	Parameter index,
	Initial parameter guess,
	Lower boundary to the parameter,
	Upper boundary to the parameter.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	19/07/2007
 ***/
void	profit_resetparams(profitstruct *profit, objstruct *obj,
		obj2struct *obj2)
  {
   double	*paramptr,
		param, parammin,parammax;
   int		p, index;


  for (p=0; p<PARAM_NPARAM; p++)
/*-- Check whether the parameter has already be registered */
    {
    param = parammin = parammax = 0.0;	/* Avoid gcc -Wall warnings*/
    switch((paramenum)p)
      {
      case PARAM_BACK:
        param = 0.0;
        parammin = -6.0*obj->sigbkg;
        parammax =  6.0*obj->sigbkg;
        break;
      case PARAM_X:
        param = obj->mx - (int)(obj->mx+0.49999);
        parammin = -obj2->hl_radius*4;
        parammax =  obj2->hl_radius*4;
        break;
      case PARAM_Y:
        param = obj->my - (int)(obj->my+0.49999);
        parammin = -obj2->hl_radius*4;
        parammax =  obj2->hl_radius*4;
        break;
      case PARAM_DEVAUC_FLUX:
      case PARAM_SERSIC_FLUX:
        param = obj2->flux_auto/2.0;
        parammin = -obj2->flux_auto;
        parammax = 2*obj2->flux_auto;
        break;
      case PARAM_EXPO_FLUX:
        param = obj2->flux_auto/2.0;
        parammin = -obj2->flux_auto;
        parammax = 2*obj2->flux_auto;
        break;
      case PARAM_DEVAUC_MAJ:
      case PARAM_SERSIC_MAJ:
        param = obj2->hl_radius;
        parammin = 0.1;
        parammax = param * 4.0;
        break;
      case PARAM_EXPO_MAJ:
        param = obj2->hl_radius/1.67835;	/* From scalelength to Re */
        parammin = param / 4.0;
        parammax = param * 4.0;
        break;
      case PARAM_DEVAUC_ASPECT:
      case PARAM_SERSIC_ASPECT:
        param = obj->b/obj->a;
        parammin = 0.5;
        parammax = 2.0;
        break;
      case PARAM_EXPO_ASPECT:
        param = obj->b/obj->a;
        parammin = 0.01;
        parammax = 100.0;
        break;
      case PARAM_DEVAUC_POSANG:
      case PARAM_EXPO_POSANG:
      case PARAM_SERSIC_POSANG:
        param = obj->theta;
        parammin = -3600.0;
        parammax =  3600.0;
        break;
      case PARAM_SERSIC_N:
        param = 2.0;
        parammin = 1.0;
        parammax = 10.0;
        break;
      case PARAM_ARMS_FLUX:
        param = obj2->flux_auto/2.0;
        parammin = 0.0;
        parammax = obj2->flux_auto*2.0;
        break;
      case PARAM_ARMS_QUADFRAC:
        param = 0.5;
        parammin = 0.0;
        parammax = 1.0;
        break;
      case PARAM_ARMS_SCALE:
        param = 1.0;
        parammin = 0.5;
        parammax = 10.0;
        break;
      case PARAM_ARMS_START:
        param = 1.0;
        parammin = 0.0;
        parammax = 5.0;
        break;
      case PARAM_ARMS_PITCH:
        param = 30.0;
        parammin = 10.0;
        parammax = 60.0;
//        if ((profit->spirindex=profit_spiralindex(profit, obj, obj2)) > 0.0)
//          {
//          param = -param;
//          parammin = -parammax;
//          parammax = -parammin;
//          }
//        printf("spiral index: %g  \n", profit->spirindex);
        break;
      case PARAM_ARMS_POSANG:
        param = 45.0;
        parammin = -3600.0;
        parammax = 3600.0;
        break;
      case PARAM_ARMS_WIDTH:
        param = 2.0;
        parammin = 1.0;
        parammax = 4.0;
        break;
      case PARAM_BAR_FLUX:
        param = obj2->flux_auto/10.0;
        parammin = 0.0;
        parammax = 2.0*obj2->flux_auto;
        break;
      case PARAM_BAR_ASPECT:
        param = 0.3;
        parammin = 0.2;
        parammax = 0.5;
        break;
      case PARAM_BAR_POSANG:
        param = 0.0;
        parammin = -3600.0;
        parammax = 3600.0;
        break;
      case PARAM_INRING_FLUX:
        param = obj2->flux_auto/10.0;
        parammin = 0.0;
        parammax = 2.0*obj2->flux_auto;
        break;
      case PARAM_INRING_WIDTH:
        param = 0.2;
        parammin = 0.0;
        parammax = 1.0;
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown profile parameter in ",
		"profit_resetparams()");
        break;
      }
    profit_resetparam(profit, (proftypenum)p, param, parammin,parammax);
    }

  return;
  }


/****** profit_resetparam ****************************************************
PROTO	void profit_resetparam(profitstruct *profit, proftypenum proftype,
		double param, double parammin, double parammax)
PURPOSE	Set the initial, lower and upper boundary values of a profile parameter.
INPUT	Pointer to the profit structure,
	Parameter index,
	Initial parameter guess,
	Lower boundary to the parameter,
	Upper boundary to the parameter.
OUTPUT	RETURN_OK if the parameter is registered, RETURN_ERROR otherwise.
AUTHOR	E. Bertin (IAP)
VERSION	25/06/2007
 ***/
int	profit_resetparam(profitstruct *profit, proftypenum proftype,
		double param, double parammin, double parammax)
  {
   double	*paramptr;
   int		index;

/* Check whether the parameter has already be registered */
  if ((paramptr=profit->paramlist[(int)proftype]))
    {
    index = paramptr - profit->param;
    profit->paraminit[index] = param;
    profit->parammin[index] = parammin;
    profit->parammax[index] = parammax;
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }

  
/****** profit_boundtounbound *************************************************
PROTO	void profit_boundtounbound(profitstruct *profit)
PURPOSE	Convert parameters from bounded to unbounded space.
INPUT	Pointer to the profit structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/05/2007
 ***/
void	profit_boundtounbound(profitstruct *profit, double *param)
  {
   double	num,den;
   int		p;

  for (p=0; p<profit->nparam; p++)
    if (profit->parammin[p]!=profit->parammax[p])
      {
      num = param[p] - profit->parammin[p];
      den = profit->parammax[p] - param[p];
      param[p] = num>1e-100? (den>1e-100? log(num/den): 200.0) : -200.0;
      }

  return;

  }


/****** profit_unboundtobound *************************************************
PROTO	void profit_unboundtobound(profitstruct *profit)
PURPOSE	Convert parameters from unbounded to bounded space.
INPUT	Pointer to the profit structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2007
 ***/
void	profit_unboundtobound(profitstruct *profit, double *param)
  {
   int		p;

  for (p=0; p<profit->nparam; p++)
    if (profit->parammin[p]!=profit->parammax[p])
      param[p] = (profit->parammax[p] - profit->parammin[p])
		/ (1.0 + exp(-(param[p]>200.0? 200.0 : param[p])))
		+ profit->parammin[p];

  return;
  }


/****** prof_init *************************************************************
PROTO	profstruct prof_init(profitstruct *profit, proftypenum profcode)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	Pointer to the profile-fitting structure,
	profile type.
OUTPUT	A pointer to an allocated prof structure.
AUTHOR	E. Bertin (IAP)
VERSION	20/07/2007
 ***/
profstruct	*prof_init(profitstruct *profit, proftypenum profcode)
  {
   profstruct	*prof;
   double	*pix,
		rmax2, re2, dy2,r2, scale, zero, k,n, hinvn;
   int		width,height, ixc,iyc, ix,iy, nsub,
		d,s;

  QCALLOC(prof, profstruct, 1);
  prof->code = profcode;
  switch(profcode)
    {
    case PROF_BACK:
      prof->naxis = 2;
      prof->pix = NULL;
      profit_addparam(profit, PARAM_BACK, &prof->flux);
      prof->typscale = 1.0;
      break;
    case PROF_SERSIC:
      prof->naxis = 3;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_SERSIC_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_SERSIC_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_SERSIC_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SERSIC_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_SERSIC_N, &prof->extra[0]);
      break;
    case PROF_DEVAUCOULEURS:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DEVAUC_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_DEVAUC_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_DEVAUC_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DEVAUC_POSANG, &prof->posangle);
      break;
    case PROF_EXPONENTIAL:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_EXPO_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_EXPO_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_EXPO_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_EXPO_POSANG, &prof->posangle);
      break;
    case PROF_ARMS:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_EXPO_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_EXPO_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_EXPO_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_ARMS_QUADFRAC, &prof->featfrac);
      profit_addparam(profit, PARAM_ARMS_SCALE, &prof->featscale);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_ARMS_PITCH, &prof->featpitch);
      profit_addparam(profit, PARAM_ARMS_POSANG, &prof->featposang);
      profit_addparam(profit, PARAM_ARMS_WIDTH, &prof->featwidth);
      break;
    case PROF_BAR:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_EXPO_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_EXPO_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_EXPO_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_BAR_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_BAR_ASPECT, &prof->feataspect);
      profit_addparam(profit, PARAM_ARMS_POSANG, &prof->featposang);
      break;
    case PROF_INRING:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_EXPO_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_EXPO_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_EXPO_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_INRING_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_INRING_WIDTH, &prof->featwidth);
      break;
    case PROF_SERSIC_TABEX:	/* An example of tabulated profile */
      prof->naxis = 3;
      width =  prof->naxisn[0] = PROFIT_PROFRES;
      height = prof->naxisn[1] = PROFIT_PROFRES;
      nsub = prof->naxisn[2] = PROFIT_PROFSRES;
      QCALLOC(prof->pix, double, width*height*nsub);
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      re2 = width/64.0;
      prof->typscale = re2;
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
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_SERSIC_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_SERSIC_MAJ, &prof->scale);
      profit_addparam(profit, PARAM_SERSIC_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SERSIC_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_SERSIC_N, &prof->extra[0]);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown profile in ",
		"prof_init()");
      break;
    }

  prof->scaling = profit->psf->pixstep / prof->typscale;
  if (prof->pix)
    {
    prof->kernelnlines = 1;
    for (d=0; d<prof->naxis; d++)
      {
      prof->interptype[d] = INTERP_BILINEAR;
      prof->kernelnlines *= 
	(prof->kernelwidth[d] = interp_kernwidth[prof->interptype[d]]);
      }
    prof->kernelnlines /= prof->kernelwidth[0];
    QMALLOC(prof->kernelbuf, double, prof->kernelnlines);
    }

  return prof;
  }  


/****** prof_end **************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile structure.
INPUT	Prof structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/04/2007
 ***/
void	prof_end(profstruct *prof)
  {
  if (prof->pix)
    {
    free(prof->pix);
    free(prof->kernelbuf);
    }
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
VERSION	23/07/2007
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
   double	posin[PROFIT_MAXEXTRA], posout[2], dnaxisn[2],
		*pixout,
		flux,fluxfac;
   float	*pixin,
		amp,ctheta,stheta,cd11,cd12,cd21,cd22, dcd11,dcd21, dx1,dx2,
		x1,x10,x2, x1cin,x2cin, x1cout,x2cout, xscale,yscale, saspect,
		x1in,x2in, odx2, ostep,
		n,k, hinvn, x1t,x2t, ca,sa, u,
		armamp,arm2amp, armrdphidr, posang, width, invwidth2,
		r, r2, rmin, r2min, r2minxin, r2minxout, rmax, r2max, invr2xdif,
		val, theta, thresh, ra,rb;
   int		npix, noversamp,
		d,e,i, ix1,ix2, idx1,idx2;

  npix = profit->modnaxisn[0]*profit->modnaxisn[1];

  if (prof->code==PROF_BACK)
    {
    amp = fabs(*prof->flux);
    pixout = profit->modpix;
    for (i=npix; i--;)
      *(pixout++) += amp;
    return;
    }

/* Compute Profile CD matrix */
  ctheta = cos(*prof->posangle*DEG);
  stheta = sin(*prof->posangle*DEG);
  saspect = sqrt(fabs(*prof->aspect));
  xscale = (*prof->scale==0.0)?
			0.0 : fabs(saspect*prof->scaling
				/ (*prof->scale*prof->typscale));
  yscale = (*prof->scale*saspect == 0.0)?
		0.0 : fabs(prof->scaling
				/ (*prof->scale*prof->typscale*saspect));
  cd11 = xscale*ctheta;
  cd12 = xscale*stheta;
  cd21 =-yscale*stheta;
  cd22 = yscale*ctheta;
  dcd11 = cd11*ostep;
  dcd21 = cd21*ostep;

  dx1 = 0.0;	/* Shifting operations have been moved to profit_resample() */
  dx2 = 0.0;	/* Shifting operations have been moved to profit_resample() */

  x1cout = (double)(profit->modnaxisn[0]/2);
  x2cout = (double)(profit->modnaxisn[1]/2);

  switch(prof->code)
    {
    case PROF_SERSIC:
      n = fabs(*prof->extra[0]);
      k = 1.0/3.0 - 2.0*n - 4.0/(405.0*n) - 46.0/(25515.0*n*n)
		- 131.0/(1148175*n*n*n);
      hinvn = 0.5/n;
/*---- The consequence of sampling on flux is compensated by PSF normalisation*/
      x10 = -x1cout - dx1 + 0.5*(ostep-1.0);
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1 = x10;
        for (ix1=profit->modnaxisn[0]; ix1--; x1+=1.0)
          {
          x1in = cd12*x2 + cd11*x1;
          x2in = cd22*x2 + cd21*x1;
          ra = x1in*x1in+x2in*x2in;
          val = expf(k*PROFIT_POWF(ra,hinvn));
          noversamp  = (int)(val*PROFIT_OVERSAMP+0.1);
          if (noversamp < 2)
            *(pixin++) = val;
          else
            {
            ostep = 1.0/noversamp;
            odx2 = 0.5*(ostep-1.0);
            val = 0.0;
            for (idx2=noversamp; idx2--; odx2+=ostep)
              {
              x1in = cd12*(x2+odx2) + cd11*x1;
              x2in = cd22*(x2+odx2) + cd21*x1;
              for (idx1=noversamp; idx1--;)
                {
                ra = x1in*x1in+x2in*x2in;
                val += expf(k*PROFIT_POWF(ra,hinvn));
                x1in += dcd11;
                x2in += dcd21;
                }
              }
            *(pixin++) = val*ostep*ostep;
            }
          }
        }
      break;
    case PROF_DEVAUCOULEURS:
/*---- The consequence of sampling on flux is compensated by PSF normalisation*/
      x1 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1in = cd12*x2 + cd11*x1;
        x2in = cd22*x2 + cd21*x1;
        for (ix1=profit->modnaxisn[0]; ix1--;)
          {
          ra = x1in*x1in+x2in*x2in;
          *(pixin++) = expf(-7.6692f*PROFIT_POWF(ra,0.125));
          x1in += cd11;
          x2in += cd21;
          }
        }
      break;
    case PROF_EXPONENTIAL:
      x1 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1in = cd12*x2 + cd11*x1;
        x2in = cd22*x2 + cd21*x1;
        for (ix1=profit->modnaxisn[0]; ix1--;)
          {
          *(pixin++) = exp(-sqrt(x1in*x1in+x2in*x2in));
          x1in += cd11;
          x2in += cd21;
          }
        }
      break;
    case PROF_ARMS:
      r2min = *prof->featstart**prof->featstart;
      r2minxin = r2min * (1.0 - PROFIT_BARXFADE) * (1.0 - PROFIT_BARXFADE);
      r2minxout = r2min * (1.0 + PROFIT_BARXFADE) * (1.0 + PROFIT_BARXFADE);
      if ((invr2xdif = (r2minxout - r2minxin)) > 0.0)
        invr2xdif = 1.0 / invr2xdif;
      else
        invr2xdif = 1.0;
      arm2amp = *prof->featfrac;
      armamp = 1.0 - arm2amp;
      armrdphidr = 1.0/tan(*prof->featpitch*DEG);
      posang = *prof->featposang*DEG;
      width = fabs(*prof->featwidth);
width = 2.0;
      x1 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=profit->modnaxisn[0]; ix1--;)
          {
          r2 = x1t*x1t+x2t*x2t;
          if (r2>r2minxin)
            {
            u = logf(r2 + 0.00001);
            theta = armrdphidr*u+posang;
            ca = cosf(theta);
            sa = sinf(theta);
            x1in = (x1t*ca - x2t*sa);
            x2in = (x1t*sa + x2t*ca);
            amp = expf(-sqrtf(x1t*x1t+x2t*x2t));
            if (r2<r2minxout)
              amp *= (r2 - r2minxin)*invr2xdif;
            ra = x1in*x1in/r2;
            rb = x2in*x2in/r2;
            *(pixin++) = amp * (armamp*PROFIT_POWF(ra,width)
				+ arm2amp*PROFIT_POWF(rb,width));
            }
          else
            *(pixin++) = 0.0;
          x1t += cd11;
          x2t += cd21; 
         }
        }
      break;
    case PROF_BAR:
      r2min = *prof->featstart**prof->featstart;
      r2minxin = r2min * (1.0 - PROFIT_BARXFADE) * (1.0 - PROFIT_BARXFADE);
      r2minxout = r2min * (1.0 + PROFIT_BARXFADE) * (1.0 + PROFIT_BARXFADE);
      if ((invr2xdif = (r2minxout - r2minxin)) > 0.0)
        invr2xdif = 1.0 / invr2xdif;
      else
        invr2xdif = 1.0;
      invwidth2 = fabs(1.0 / (*prof->featstart**prof->feataspect));
      posang = *prof->featposang*DEG;
      ca = cosf(posang);
      sa = sinf(posang);
      x1 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=profit->modnaxisn[0]; ix1--;)
          {
          r2 = x1t*x1t+x2t*x2t;
          if (r2<r2minxout)
            {
            x1in = x1t*ca - x2t*sa;
            x2in = invwidth2*(x1t*sa + x2t*ca);
            *(pixin++) = (r2>r2minxin) ?
				(r2minxout - r2)*invr2xdif*expf(-x2in*x2in)
				: expf(-x2in*x2in);
            }
          else
            *(pixin++) = 0.0;
          x1t += cd11;
          x2t += cd21;
          }
        }
      break;
    case PROF_INRING:
      rmin = *prof->featstart;
      r2minxin = *prof->featstart-4.0**prof->featwidth;
      if (r2minxin < 0.0)
        r2minxin = 0.0;
      r2minxin *= r2minxin;
      r2minxout = *prof->featstart+4.0**prof->featwidth;
      r2minxout *= r2minxout;
      invwidth2 = 0.5 / (*prof->featwidth**prof->featwidth);
      x1 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = profit->pmodpix;
      for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=profit->modnaxisn[0]; ix1--;)
          {
          r2 = x1t*x1t+x2t*x2t;
          if (r2>r2minxin && r2<r2minxout)
            {
            r = sqrt(r2) - rmin;
            *(pixin++) = expf(-invwidth2*r*r);
            }
          else
            *(pixin++) = 0.0;
          x1t += cd11;
          x2t += cd21;
          }
        }
      break;
    default:
/*---- Tabulated profile: remap each pixel */
/*---- Initialize multi-dimensional counters */
      for (d=0; d<2; d++)
        {
        posout[d] = 1.0;	/* FITS convention */
        dnaxisn[d] = profit->modnaxisn[d] + 0.99999;
        }

      for (e=0; e<prof->naxis - 2; e++)
        {
        d = 2+e;
/*------ Compute position along axis */
        posin[d] = (*prof->extra[e]-prof->extrazero[e])/prof->extrascale[e]+1.0;
/*------ Keep position within boundaries and let interpolation do the rest */
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
      x1cin = (double)(prof->naxisn[0]/2);
      x2cin = (double)(prof->naxisn[1]/2);
      pixin = profit->pmodpix;
      for (i=npix; i--;)
        {
        x1 = posout[0] - x1cout - 1.0 - dx1;
        x2 = posout[1] - x2cout - 1.0 - dx2;
        posin[0] = cd11*x1 + cd12*x2 + x1cin + 1.0;
        posin[1] = cd21*x1 + cd22*x2 + x2cin + 1.0;
        *(pixin++) = prof_interpolate(prof, posin);
        for (d=0; d<2; d++)
          if ((posout[d]+=1.0) < dnaxisn[d])
            break;
          else
            posout[d] = 1.0;
        }
    break;
    }

/* Now find truncation threshold */
/* Find the shortest distance to a vignet border */
  rmax = profit->modnaxisn[0] - x1cout;
  if (rmax > (r = profit->modnaxisn[1] - x2cout))
    rmax = r;
  rmax -= 0.99;
  if (rmax<1.0)
    rmax = 1.0;
  r2max = rmax*rmax;
  rmin = rmax - 1.0;
  r2min = rmin*rmin;

/* Find best threshold (the max around the circle with radius rmax */
  dx2 = -x2cout;
  pixin = profit->pmodpix;
  if (*prof->flux>=0.0)
    {
    thresh = -BIG;
    for (ix2=profit->modnaxisn[1]; ix2--; dx2 += 1.0)
      {
      dx1 = -x1cout;
      for (ix1=profit->modnaxisn[0]; ix1--; dx1 += 1.0)
        if ((val=*(pixin++))>thresh && (r=dx1*dx1+dx2*dx2)>r2min && r2<r2max)
          thresh = val;
      }
    }
  else
    {
    thresh = BIG;
    for (ix2=profit->modnaxisn[1]; ix2--; dx2 += 1.0)
      {
      dx1 = -x1cout;
      for (ix1=profit->modnaxisn[0]; ix1--; dx1 += 1.0)
        if ((val=*(pixin++))<thresh && (r=dx1*dx1+dx2*dx2)>r2min && r2<r2max)
          thresh = val;
      }
    }

/* Threshold and measure the flux */
  flux = 0.0;
  pixin = profit->pmodpix;
  for (n=npix; n--; pixin++)
    if (*pixin >= thresh)
      flux += *pixin;
    else
      *pixin = 0.0;

/* Correct final flux */
  fluxfac = fabs(flux)>0.0? *prof->flux / flux : 1.0;
  pixin = profit->pmodpix;
  pixout = profit->modpix;
  for (n=npix; n--;)
    *(pixout++) += fluxfac * *(pixin++);

  return;
  }


/****i* gammln ***************************************************************
PROTO	double gammln(double xx)
PURPOSE	Returns the log of the Gamma function (from Num. Recipes in C, p.168).
INPUT	A double.
OUTPUT	Log of the Gamma function.
NOTES	-.
AUTHOR	E. Bertin (IAP
VERSION	29/10/97
*/
static double	gammln(double xx)

  {
   double               x,tmp,ser;
   static double        cof[6]={76.18009173,-86.50532033,24.01409822,
                        -1.231739516,0.120858003e-2,-0.536382e-5};
   int                  j;

  tmp=(x=xx-1.0)+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<6;j++)
    ser += cof[j]/(x+=1.0);

  return log(2.50662827465*ser)-tmp;
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
    ival = (interptype==INTERP_NEARESTNEIGHBOUR)? (int)(val-0.50001):(int)val;
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

