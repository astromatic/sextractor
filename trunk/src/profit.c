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
*	Last modify:	20/03/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#define _GNU_SOURCE
#include <math.h>
#endif

#ifdef HAVE_LOGF
#define	LOGF	logf
#else
#define	LOGF	log
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"levmar/lm.h"
#include	"fft.h"
#include	"fitswcs.h"
#include	"check.h"
#include	"pattern.h"
#include	"psf.h"
#include	"profit.h"

static double	prof_interpolate(profstruct *prof, double *posin);
static double	interpolate_pix(double *posin, double *pix, int *naxisn,
		interpenum interptype);

static void	make_kernel(double pos, double *kernel, interpenum interptype);

/*------------------------------- variables ---------------------------------*/

char		profname[][32]={"background offset", "Sersic spheroid",
		"De Vaucouleurs spheroid", "exponential disk", "spiral arms",
		"bar", "inner ring", "outer ring", "tabulated model",
		""};

int		interp_kernwidth[5]={1,2,4,6,8};
int theniter, the_gal;
/* "Local" global variables; it seems dirty but it simplifies a lot */
/* interfacing to the LM routines */
static picstruct	*the_field, *the_wfield;
profitstruct		*theprofit;

/****** profit_init ***********************************************************
PROTO	profitstruct profit_init(psfstruct *psf)
PURPOSE	Allocate and initialize a new profile-fitting structure.
INPUT	Pointer to PSF structure.
OUTPUT	A pointer to an allocated profit structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/04/2008
 ***/
profitstruct	*profit_init(psfstruct *psf)
  {
   profitstruct		*profit;
   int			p, nprof,
			backflag, spheroidflag, diskflag, barflag, armsflag;

  QCALLOC(profit, profitstruct, 1);
  profit->psf = psf;
  profit->psfdft = NULL;

  profit->nparam = 0;
  QMALLOC(profit->prof, profstruct *, PROF_NPROF);
  backflag = spheroidflag = diskflag = barflag = armsflag = 0;
  nprof = 0;
  for (p=0; p<PROF_NPROF; p++)
    if (!backflag && FLAG(obj2.prof_offset_flux))
      {
      profit->prof[p] = prof_init(profit, PROF_BACK);
      backflag = 1;
      nprof++;
      }
    else if (!spheroidflag && FLAG(obj2.prof_spheroid_flux))
      {
      profit->prof[p] = prof_init(profit,
	FLAG(obj2.prof_spheroid_sersicn)? PROF_SERSIC : PROF_DEVAUCOULEURS);
      spheroidflag = 1;
      nprof++;
      }
    else if (!diskflag && FLAG(obj2.prof_disk_flux))
      {
      profit->prof[p] = prof_init(profit, PROF_EXPONENTIAL);
      diskflag = 1;
      nprof++;
      }
    else if (diskflag && !barflag && FLAG(obj2.prof_bar_flux))
      {
      profit->prof[p] = prof_init(profit, PROF_BAR);
      barflag = 1;
      nprof++;
      }
    else if (barflag && !armsflag && FLAG(obj2.prof_arms_flux))
      {
      profit->prof[p] = prof_init(profit, PROF_ARMS);
      armsflag = 1;
      nprof++;
      }

  QMALLOC(profit->covar, double, profit->nparam*profit->nparam);
  profit->nprof = nprof;

  return profit;
  }  


/****** profit_end ************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile-fitting structure.
INPUT	Prof structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/04/2008
 ***/
void	profit_end(profitstruct *profit)
  {
   int	p;

  for (p=0; p<profit->nprof; p++)
    prof_end(profit->prof[p]);
  free(profit->prof);
  free(profit->covar);
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
VERSION	20/03/2009
 ***/
void	profit_fit(profitstruct *profit,
		picstruct *field, picstruct *wfield,
		objstruct *obj, obj2struct *obj2)
  {
    profitstruct	pprofit;
    patternstruct *pattern;
    psfstruct		*psf;
    checkstruct		*check;
    double		*oldparaminit,
			psf_fwhm, oldchi2, a , cp,sp, emx2,emy2,emxy, dchi2;
    int			i,j,p, oldniter, nparam, ncomp;

  nparam = profit->nparam;
  if (profit->psfdft)
    {
    QFREE(profit->psfdft);
    }

  psf = profit->psf;
  profit->pixstep = psf->pixstep;

/* Create pixmaps at image resolution */
  psf_fwhm = psf->masksize[0]*psf->pixstep;
  profit->objnaxisn[0] = (((int)((obj->xmax-obj->xmin+1) + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
  profit->objnaxisn[1] = (((int)((obj->ymax-obj->ymin+1) + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
  profit->ix = (int)(obj->mx + 0.49999);/* internal convention: 1st pix = 0 */
  profit->iy = (int)(obj->my + 0.49999);/* internal convention: 1st pix = 0 */

  if (profit->objnaxisn[1]<profit->objnaxisn[0])
    profit->objnaxisn[1] = profit->objnaxisn[0];
  else
    profit->objnaxisn[0] = profit->objnaxisn[1];

/* Use (dirty) global variables to interface with lmfit */
  the_field = field;
  the_wfield = wfield;
  theprofit = profit;
  profit->obj = obj;
  profit->obj2 = obj2;

  QMALLOC(profit->objpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  QMALLOC(profit->objweight, PIXTYPE,profit->objnaxisn[0]*profit->objnaxisn[1]);
  QMALLOC(profit->lmodpix, PIXTYPE, profit->objnaxisn[0]*profit->objnaxisn[1]);
  profit->nresi = profit_copyobjpix(profit, field, wfield);
  if (profit->nresi < nparam)
    {
    if (FLAG(obj2.prof_vector))
      for (p=0; p<nparam; p++)
        obj2->prof_vector[p] = 0.0;
    obj2->prof_niter = 0;
    return;
    }

  QMALLOC(profit->resi, double, profit->nresi);

/* Create pixmap at PSF resolution */
  profit->modnaxisn[0] =
	((int)(profit->objnaxisn[0]/profit->pixstep +0.4999)/2+1)*2; 
  profit->modnaxisn[1] =
	((int)(profit->objnaxisn[1]/profit->pixstep +0.4999)/2+1)*2; 
  if (profit->modnaxisn[1] < profit->modnaxisn[0])
    profit->modnaxisn[1] = profit->modnaxisn[0];
  else
    profit->modnaxisn[0] = profit->modnaxisn[1];

/* Allocate memory for the complete model */
  QCALLOC(profit->modpix, double, profit->modnaxisn[0]*profit->modnaxisn[1]);
  QMALLOC(profit->psfpix, double, profit->modnaxisn[0]*profit->modnaxisn[1]);
/* Allocate memory for the partial model */
  QMALLOC(profit->pmodpix, float, profit->modnaxisn[0]*profit->modnaxisn[1]);

/* Compute the local PSF */
  profit_psf(profit);

/* Set initial guesses and boundaries */
  obj2->prof_flag = 0;
  profit->sigma = obj->sigbkg;

  profit_resetparams(profit);

the_gal++;

/* Actual minimisation */
  profit->niter = profit_minimize(profit, PROFIT_MAXITER);
  profit_residuals(profit,field,wfield, 0.0, profit->param,profit->resi);

/*
  QMEMCPY(profit->paraminit, oldparaminit, double, nparam);
  if (profit_setparam(profit, PARAM_ARMS_PITCH, 160.0, 130.0, 175.0)==RETURN_OK)
    {
    oldchi2 = profit->chi2;
    oldniter = profit->niter;
    profit_resetparams(profit);
    profit_setparam(profit, PARAM_ARMS_PITCH, 160.0, 130.0, 175.0);
    profit->niter = profit_minimize(profit, PROFIT_MAXITER);
    if (profit->chi2 > oldchi2)
      {
      memcpy(profit->paraminit, oldparaminit, nparam*sizeof(double));
      profit->chi2 = oldchi2;
      profit->niter = oldniter;
      }
    else
      obj2->prof_flag |= PROFIT_FLIPPED;
    }  
*/

/* Convert covariance matrix to bound space */
  profit_covarunboundtobound(profit);
  for (p=0; p<nparam; p++)
    profit->paramerr[p]= sqrt(profit->covar[p*(nparam+1)]);

/* Equate param and paraminit vectors to avoid confusion later on */
  for (p=0; p<profit->nparam; p++)
    profit->param[p] = profit->paraminit[p];

/* CHECK-Images */
  if ((check = prefs.check[CHECK_SUBPROFILES]))
    {
    profit_residuals(profit,field,wfield, 0.0, profit->param,profit->resi);
    addcheck(check, profit->lmodpix, profit->objnaxisn[0],profit->objnaxisn[1],
		profit->ix,profit->iy, -1.0);
    }
  if ((check = prefs.check[CHECK_PROFILES]))
    {
    profit_residuals(profit,field,wfield, 0.0, profit->param,profit->resi);
    addcheck(check, profit->lmodpix, profit->objnaxisn[0],profit->objnaxisn[1],
		profit->ix,profit->iy, 1.0);
    }
/* Fill measurement parameters */
  if (FLAG(obj2.prof_vector))
    {
    for (p=0; p<nparam; p++)
      obj2->prof_vector[p]= profit->param[p];
    }
  if (FLAG(obj2.prof_errvector))
    {
    for (p=0; p<nparam; p++)
      obj2->prof_errvector[p]= profit->paramerr[p];
    }

  obj2->prof_niter = profit->niter;
  obj2->flux_prof = profit->flux;
  obj2->prof_chi2 = (profit->nresi > profit->nparam)?
		profit->chi2 / (profit->nresi - profit->nparam) : 0.0;

  if (FLAG(obj2.x_prof))
    {
    i = profit->paramindex[PARAM_X];
    j = profit->paramindex[PARAM_Y];
/*-- Model coordinates follow the FITS convention (first pixel at 1,1) */
    if (profit->paramlist[PARAM_X])
      {
      obj2->x_prof = profit->ix + *profit->paramlist[PARAM_X] + 1.0;
      obj2->poserrmx2_prof = emx2 = profit->covar[i*(nparam+1)];
      }
    else
      emx2 = 0.0;
    if (profit->paramlist[PARAM_Y])
      {
      obj2->y_prof = profit->iy + *profit->paramlist[PARAM_Y] + 1.0;
      obj2->poserrmy2_prof = emy2 = profit->covar[j*(nparam+1)];
      }
    else
      emy2 = 0.0;
    if (profit->paramlist[PARAM_X] && profit->paramlist[PARAM_Y])
      obj2->poserrmxy_prof = emxy = profit->covar[i+j*nparam];
    else
      emxy = 0.0;

/*-- Error ellipse parameters */
    if (FLAG(obj2.poserra_prof))
      {
       double	pmx2,pmy2,temp,theta;

      if (fabs(temp=emx2-emy2) > 0.0)
        theta = atan2(2.0 * emxy,temp) / 2.0;
      else
        theta = PI/4.0;

      temp = sqrt(0.25*temp*temp+ emxy*emxy);
      pmy2 = pmx2 = 0.5*(emx2+emy2);
      pmx2+=temp;
      pmy2-=temp;

      obj2->poserra_prof = (float)sqrt(pmx2);
      obj2->poserrb_prof = (float)sqrt(pmy2);
      obj2->poserrtheta_prof = theta*180.0/PI;
      }

    if (FLAG(obj2.poserrcxx_prof))
      {
       double	temp;

      obj2->poserrcxx_prof = (float)(emy2/(temp=emx2*emy2-emxy*emxy));
      obj2->poserrcyy_prof = (float)(emx2/temp);
      obj2->poserrcxy_prof = (float)(-2*emxy/temp);
      }
    }


  if (FLAG(obj2.prof_mx2))
    {
    memset(profit->modpix, 0,
	profit->modnaxisn[0]*profit->modnaxisn[1]*sizeof(double));
    for (p=0; p<profit->nprof; p++)
      prof_add(profit->prof[p], profit);
    profit_moments(profit);
    }

/* Bulge */
  if (FLAG(obj2.prof_spheroid_flux))
    {
    obj2->prof_spheroid_flux = *profit->paramlist[PARAM_SPHEROID_FLUX];
    obj2->prof_spheroid_fluxerr =
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_FLUX]];
    obj2->prof_spheroid_reff = *profit->paramlist[PARAM_SPHEROID_REFF];
    obj2->prof_spheroid_refferr = 
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_REFF]];
    obj2->prof_spheroid_aspect = *profit->paramlist[PARAM_SPHEROID_ASPECT];
    obj2->prof_spheroid_aspecterr = 
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_ASPECT]];
    obj2->prof_spheroid_theta =
			fmod_m90_p90(*profit->paramlist[PARAM_SPHEROID_POSANG]);
    obj2->prof_spheroid_thetaerr = 
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_POSANG]];
    if (FLAG(obj2.prof_spheroid_sersicn))
      {
      obj2->prof_spheroid_sersicn = *profit->paramlist[PARAM_SPHEROID_SERSICN];
      obj2->prof_spheroid_sersicnerr = 
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_SERSICN]];
      }
    }

/* Disk */
  if (FLAG(obj2.prof_disk_flux))
    {
    obj2->prof_disk_flux = *profit->paramlist[PARAM_DISK_FLUX];
    obj2->prof_disk_fluxerr =
		profit->paramerr[profit->paramindex[PARAM_DISK_FLUX]];
    obj2->prof_disk_scale = *profit->paramlist[PARAM_DISK_SCALE];
    obj2->prof_disk_scaleerr =
		profit->paramerr[profit->paramindex[PARAM_DISK_SCALE]];
    obj2->prof_disk_aspect = *profit->paramlist[PARAM_DISK_ASPECT];
    obj2->prof_disk_aspecterr =
		profit->paramerr[profit->paramindex[PARAM_DISK_ASPECT]];
    obj2->prof_disk_theta = fmod_m90_p90(*profit->paramlist[PARAM_DISK_POSANG]);
    obj2->prof_disk_thetaerr =
		profit->paramerr[profit->paramindex[PARAM_DISK_POSANG]];
    if (FLAG(obj2.prof_disk_inclination))
      {
      obj2->prof_disk_inclination = acos(obj2->prof_disk_aspect) / DEG;
      if (FLAG(obj2.prof_disk_inclinationerr))
        {
        a = sqrt(1.0-obj2->prof_disk_aspect*obj2->prof_disk_aspect);
        obj2->prof_disk_inclinationerr = obj2->prof_disk_aspecterr
					/(a>0.1? a : 0.1)/DEG;
        }
      }

/* Disk pattern */
    if (prefs.pattern_flag)
      {
      profit_residuals(profit,field,wfield, PROFIT_DYNPARAM,
			profit->param,profit->resi);
      pattern = pattern_init(profit, prefs.pattern_type,
		prefs.prof_disk_patternncomp);
      pattern_fit(pattern, profit);
      if (FLAG(obj2.prof_disk_patternspiral))
        obj2->prof_disk_patternspiral = pattern_spiral(pattern);
      if (FLAG(obj2.prof_disk_patternvector))
        {
        ncomp = pattern->size[2];
        for (p=0; p<ncomp; p++)
          obj2->prof_disk_patternvector[p] = (float)pattern->coeff[p];
        }
      if (FLAG(obj2.prof_disk_patternmodvector))
        {
        ncomp = pattern->ncomp*pattern->nfreq;
        for (p=0; p<ncomp; p++)
          obj2->prof_disk_patternmodvector[p] = (float)pattern->mcoeff[p];
        }
      if (FLAG(obj2.prof_disk_patternargvector))
        {
        ncomp = pattern->ncomp*pattern->nfreq;
        for (p=0; p<ncomp; p++)
          obj2->prof_disk_patternargvector[p] = (float)pattern->acoeff[p];
        }
      pattern_end(pattern);
      }

/* Bar */
    if (FLAG(obj2.prof_bar_flux))
      {
      obj2->prof_bar_flux = *profit->paramlist[PARAM_BAR_FLUX];
      obj2->prof_bar_fluxerr =
		profit->paramerr[profit->paramindex[PARAM_BAR_FLUX]];
      obj2->prof_bar_length = *profit->paramlist[PARAM_ARMS_START]
				**profit->paramlist[PARAM_DISK_SCALE];
      obj2->prof_bar_lengtherr = *profit->paramlist[PARAM_ARMS_START]
		  * profit->paramerr[profit->paramindex[PARAM_DISK_SCALE]]
		+ *profit->paramlist[PARAM_DISK_SCALE]
		  * profit->paramerr[profit->paramindex[PARAM_ARMS_START]];
      obj2->prof_bar_aspect = *profit->paramlist[PARAM_BAR_ASPECT];
      obj2->prof_bar_aspecterr =
		profit->paramerr[profit->paramindex[PARAM_BAR_ASPECT]];
      obj2->prof_bar_posang = 
			fmod_m90_p90(*profit->paramlist[PARAM_ARMS_POSANG]);
      obj2->prof_bar_posangerr =
		profit->paramerr[profit->paramindex[PARAM_ARMS_POSANG]];
      if (FLAG(obj2.prof_bar_theta))
        {
        cp = cos(obj2->prof_bar_posang*DEG);
        sp = sin(obj2->prof_bar_posang*DEG);
        a = obj2->prof_disk_aspect;
        obj2->prof_bar_theta = fmod_m90_p90(atan2(a*sp,cp)/DEG
				+ obj2->prof_disk_theta);
        obj2->prof_bar_thetaerr = obj2->prof_bar_posangerr*a/(cp*cp+a*a*sp*sp);
        }

/* Arms */
      if (FLAG(obj2.prof_arms_flux))
        {
        obj2->prof_arms_flux = *profit->paramlist[PARAM_ARMS_FLUX];
        obj2->prof_arms_fluxerr =
		profit->paramerr[profit->paramindex[PARAM_ARMS_FLUX]];
        obj2->prof_arms_pitch =
		fmod_m90_p90(*profit->paramlist[PARAM_ARMS_PITCH]);
        obj2->prof_arms_pitcherr =
		profit->paramerr[profit->paramindex[PARAM_ARMS_PITCH]];
        obj2->prof_arms_start = *profit->paramlist[PARAM_ARMS_START]
				**profit->paramlist[PARAM_DISK_SCALE];
        obj2->prof_arms_starterr = *profit->paramlist[PARAM_ARMS_START]
		  * profit->paramerr[profit->paramindex[PARAM_DISK_SCALE]]
		+ *profit->paramlist[PARAM_DISK_SCALE]
		  * profit->paramerr[profit->paramindex[PARAM_ARMS_START]];
        obj2->prof_arms_quadfrac = *profit->paramlist[PARAM_ARMS_QUADFRAC];
        obj2->prof_arms_quadfracerr =
		profit->paramerr[profit->paramindex[PARAM_ARMS_QUADFRAC]];
        obj2->prof_arms_posang =
			fmod_m90_p90(*profit->paramlist[PARAM_ARMS_POSANG]);
        obj2->prof_arms_posangerr =
		profit->paramerr[profit->paramindex[PARAM_ARMS_POSANG]];
        }
      }
    }

  if (FLAG(obj2.prof_class_star))
    {
    pprofit = *profit;
    memset(pprofit.paramindex, 0, PARAM_NPARAM*sizeof(int));
    memset(pprofit.paramlist, 0, PARAM_NPARAM*sizeof(double *));
    pprofit.nparam = 0;
    QMALLOC(pprofit.prof, profstruct *, 1);
    pprofit.prof[0] = prof_init(&pprofit, PROF_DIRAC);
    QMALLOC(pprofit.covar, double, pprofit.nparam*pprofit.nparam);
    pprofit.nprof = 1;
    profit_resetparams(&pprofit);
    if (profit->paramlist[PARAM_X] && profit->paramlist[PARAM_Y])
      {
      pprofit.paraminit[pprofit.paramindex[PARAM_X]] = *profit->paramlist[PARAM_X];
      pprofit.paraminit[pprofit.paramindex[PARAM_Y]] = *profit->paramlist[PARAM_Y];
      }
    pprofit.paraminit[pprofit.paramindex[PARAM_DISK_FLUX]] = profit->flux;
    pprofit.niter = profit_minimize(&pprofit, PROFIT_MAXITER);
    profit_residuals(&pprofit,field,wfield, 0.0, pprofit.param,pprofit.resi);
    dchi2 = 0.5*(pprofit.chi2 - profit->chi2);
    obj2->prof_class_star = dchi2 < 50.0?
	(dchi2 > -50.0? 2.0/(1.0+exp(dchi2)) : 2.0) : 0.0;
    if (profit->flux > 0.0 && pprofit.flux > 0.0)
      obj2->prof_concentration = -2.5*log10(pprofit.flux / profit->flux);
    else  if (profit->flux > 0.0)
      obj2->prof_concentration = 99.0;
    else  if (pprofit.flux > 0.0)
      obj2->prof_concentration = -99.0;
    prof_end(pprofit.prof[0]);
    free(pprofit.prof);
    free(pprofit.covar);
    }

/* clean up. */
  free(profit->modpix);
  free(profit->psfpix);
  free(profit->pmodpix);
  free(profit->lmodpix);
  free(profit->objpix);
  free(profit->objweight);
  free(profit->resi);
/*
  free(oldparaminit);
*/
  return;
  }


/****** profit_psf ************************************************************
PROTO	void	profit_psf(profitstruct *profit)
PURPOSE	Build the local PSF at a given resolution.
INPUT	Profile-fitting structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	22/10/2008
 ***/
void	profit_psf(profitstruct *profit)
  {
   double	posin[2], posout[2], dnaxisn[2],
		*pixout,
		xcout,ycout, xcin,ycin, invpixstep, flux, norm;
   int		d,i;

  psf = profit->psf;
  psf_build(psf);

  xcout = (double)(profit->modnaxisn[0]/2) + 1.0;	/* FITS convention */
  ycout = (double)(profit->modnaxisn[1]/2) + 1.0;	/* FITS convention */
  xcin = (psf->masksize[0]/2) + 1.0;			/* FITS convention */
  ycin = (psf->masksize[1]/2) + 1.0;			/* FITS convention */
  invpixstep = profit->pixstep / psf->pixstep;

/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 1.0;					/* FITS convention */
    dnaxisn[d] = profit->modnaxisn[d]+0.5;
    }

/* Remap each pixel */
  pixout = profit->psfpix;
  flux = 0.0;
  for (i=profit->modnaxisn[0]*profit->modnaxisn[1]; i--;)
    {
    posin[0] = (posout[0] - xcout)*invpixstep + xcin;
    posin[1] = (posout[1] - ycout)*invpixstep + ycin;
    flux += ((*(pixout++) = interpolate_pix(posin, psf->maskloc,
		psf->masksize, INTERP_LANCZOS3)));
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 1.0;
    }

/* Normalize PSF flux (just in case...) */
  flux *= psf->pixstep*psf->pixstep;
  if (fabs(flux) > 0.0)
    {
    norm = 1.0/flux;
    pixout = profit->psfpix;
    for (i=profit->modnaxisn[0]*profit->modnaxisn[1]; i--;)
      *(pixout++) *= norm;
    }

  return;
  }


/****** profit_findinit *******************************************************
PROTO	void profit_findinit(profitstruct *profit)
PURPOSE	Find a suitable set of initialisation parameters
INPUT	Pointer to the profit structure involved in the fit.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2008
 ***/
void	profit_findinit(profitstruct *profit)
  {
   int	p;

  for (p=0; p<profit->nprof; p++)
    switch (profit->prof[p]->code)
      {
      default:
      break;
      }

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
VERSION	23/05/2008
 ***/
int	profit_minimize(profitstruct *profit, int niter)
  {
   double		lm_opts[5], info[LM_INFO_SZ];
   int			m,n;

/* Allocate work space */
  n = profit->nparam;
  m = profit->nresi;

  memset(profit->resi, 0, profit->nresi*sizeof(double));
  memset(profit->covar, 0, profit->nparam*profit->nparam*sizeof(double));
  profit_boundtounbound(profit, profit->paraminit);

/* Perform fit */
  lm_opts[0] = 1.0e-3;
  lm_opts[1] = 1.0e-17;
  lm_opts[2] = 1.0e-17;
  lm_opts[3] = 1.0e-17;
  lm_opts[4] = 1.0e-6;

  niter = dlevmar_dif(profit_evaluate, profit->paraminit, profit->resi,
	n, m, niter, lm_opts, info, NULL, profit->covar, profit);

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
VERSION	17/09/2008
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
		profit->ix,profit->iy, 1.0);

    reendcheck(the_field, check);
    endcheck(check);
    }

  return;
  }


/****** profit_evaluate ******************************************************
PROTO	void profit_evaluate(double *par, double *fvec, int m, int n,
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
VERSION	18/09/2008
 ***/
void	profit_evaluate(double *par, double *fvec, int m, int n,
			void *adata)
  {
   profitstruct	*profit;

  profit = (profitstruct *)adata;
  profit_unboundtobound(profit, par);
  profit_residuals(profit, the_field, the_wfield, PROFIT_DYNPARAM, par, fvec);
  profit_boundtounbound(profit, par);
  profit_printout(m, par, n, fvec, adata, 0, -1, 0 );
  return;
  }


/****** profit_residuals ******************************************************
PROTO	double *prof_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, double dynparam, double *param, double *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	pointer to the field,
	pointer to the field weight,
	dynamic compression parameter (0=no compression),
	pointer to the model parameters (output),
	pointer to the computed residuals (output).
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/03/2009
 ***/
double	*profit_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, double dynparam, double *param, double *resi)
  {
   int		p;

  memset(profit->modpix, 0,
	profit->modnaxisn[0]*profit->modnaxisn[1]*sizeof(double));
  for (p=0; p<profit->nparam; p++)
    profit->param[p] = param[p];
/* Simple PSF shortcut */
  if (profit->nprof == 1 && profit->prof[0]->code == PROF_DIRAC)
    profit_resample(profit, profit->psfpix, profit->lmodpix,
		*profit->prof[0]->flux);
  else
    {
    for (p=0; p<profit->nprof; p++)
      prof_add(profit->prof[p], profit);
    profit_convolve(profit, profit->modpix);
    profit_resample(profit, profit->modpix, profit->lmodpix, 1.0);
    }
  profit_compresi(profit, dynparam, resi);

  return resi;
  }


/****** profit_compresi ******************************************************
PROTO	double *prof_compresi(profitstruct *profit, double dynparam,
			double *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	dynamic-compression parameter (0=no compression),
	vector of residuals (output).
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/03/2009
 ***/
double	*profit_compresi(profitstruct *profit, double dynparam, double *resi)
  {
   double	*resit,
		error, x1c,x1,x2,rmin;
   PIXTYPE	*objpix, *objweight, *lmodpix,
		val,val2,wval, invsig;
   int		npix, i;
  
/* Compute vector of residuals */
  resit = resi;
  objpix = profit->objpix;
  objweight = profit->objweight;
  lmodpix = profit->lmodpix;
  error = 0.0;
  x1c = (double)(profit->objnaxisn[0]/2);
  rmin = profit->obj2->hl_radius / 2.0;
  x2 = -(double)(profit->objnaxisn[1]/2);
  npix = profit->objnaxisn[0]*profit->objnaxisn[1];
  if (dynparam > 0.0)
    {
    invsig = (PIXTYPE)(1.0/dynparam);
    for (i=npix; i--; lmodpix++)
      {
      val = *(objpix++);
      if ((wval=*(objweight++))>0.0)
        {
        val2 = (val - *lmodpix)*wval*invsig;
        val2 = val2>0.0? LOGF(1.0+val2) : -LOGF(1.0-val2);
        *(resit++) = val2*dynparam;
        error += val2*val2;
        }
      }
    profit->chi2 = dynparam*dynparam*error;
    }
  else
    {
    for (i=npix; i--; lmodpix++)
      {
      val = *(objpix++);
      if ((wval=*(objweight++))>0.0)
        {
        val2 = (val - *lmodpix)*wval;
        *(resit++) = val2;
        error += val2*val2;
        }
      }
    profit->chi2 = error;
    }

  return resi;
  }


/****** profit_resample ******************************************************
PROTO	void	prof_resample(profitstruct *profit, double *inpix,
		PIXTYPE *outpix)
PURPOSE	Resample the current full resolution model to image resolution.
INPUT	Profile-fitting structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/03/2009
 ***/
void	profit_resample(profitstruct *profit, double *inpix, PIXTYPE *outpix,
	double factor)
  {
   double	posin[2], posout[2], dnaxisn[2],
		*dx,*dy,
		xcout,ycout, xcin,ycin, invpixstep, flux;
   int		d,i;

  xcout = (double)(profit->objnaxisn[0]/2) + 1.0;	/* FITS convention */
  if ((dx=(profit->paramlist[PARAM_X])))
    xcout += *dx;
  ycout = (double)(profit->objnaxisn[1]/2) + 1.0;	/* FITS convention */
  if ((dy=(profit->paramlist[PARAM_Y])))
    ycout += *dy;
  xcin = (profit->modnaxisn[0]/2) + 1.0;		/* FITS convention */
  ycin = (profit->modnaxisn[1]/2) + 1.0;		/* FITS convention */
  invpixstep = 1.0/profit->pixstep;

/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 1.0;					/* FITS convention */
    dnaxisn[d] = profit->objnaxisn[d]+0.5;
    }

/* Remap each pixel */
  flux = 0.0;
  for (i=profit->objnaxisn[0]*profit->objnaxisn[1]; i--;)
    {
    posin[0] = (posout[0] - xcout)*invpixstep + xcin;
    posin[1] = (posout[1] - ycout)*invpixstep + ycin;
    flux += ((*(outpix++) = (PIXTYPE)(factor*interpolate_pix(posin, inpix,
		profit->modnaxisn, INTERP_LANCZOS3))));
    for (d=0; d<2; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 1.0;
    }

  profit->flux = flux;

  return;
  }


/****** profit_convolve *******************************************************
PROTO	void profit_convolve(profitstruct *profit, double *modpix)
PURPOSE	Convolve a model image with the local PSF.
INPUT	Pointer to the profit structure,
	Pointer to the image raster.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/09/2008
 ***/
void	profit_convolve(profitstruct *profit, double *modpix)
  {
  if (!profit->psfdft)
    profit_makedft(profit);

  fft_conv(modpix, profit->psfdft, profit->modnaxisn);

  return;
  }


/****** profit_makedft *******************************************************
PROTO	void profit_makedft(profitstruct *profit)
PURPOSE	Create the Fourier transform of the descrambled PSF component.
INPUT	Pointer to the profit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	22/04/2008
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

  psfwidth = profit->modnaxisn[0];
  psfheight = profit->modnaxisn[1];
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
  ppix = profit->psfpix + (psfheight/2)*psfwidth + psfwidth/2;
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

  ppix = profit->psfpix + ((psfheight/2)-hcpheight)*psfwidth + psfwidth/2;
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
  rsig = psf->fwhm/profit->pixstep;
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
			picstruct *wfield)
PURPOSE	Copy a piece of the input field image to a profit structure.
INPUT	Pointer to the profit structure,
	Pointer to the field structure,
	Pointer to the field weight structure.
OUTPUT	The number of valid pixels copied.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	18/09/2008
 ***/
int	profit_copyobjpix(profitstruct *profit, picstruct *field,
			picstruct *wfield)
  {
   double	dx2, dy2, dr2, rad2;
   PIXTYPE	*pixin,*pixout, *wpixin,*wpixout,
		backnoise2, invgain, satlevel, wthresh, pix, wpix;
   int		i,x,y, xmin,xmax,ymin,ymax, w,h,w2,dw, npix, off, gainflag,
		ix,iy;

/* First put the image background to -BIG */
  pixout = profit->objpix;
  wpixout = profit->objweight;
  w = profit->objnaxisn[0];
  h = profit->objnaxisn[1];
  for (i=w*h; i--;)
    {
    *(pixout++) = -BIG;
    *(wpixout++) = 0.0;
    }

/* Don't go further if out of frame!! */
  ix = profit->ix;
  iy = profit->iy;
  if (ix<0 || ix>=field->width || iy<field->ymin || iy>=field->ymax)
    return 0;

  backnoise2 = field->backsig*field->backsig;
  invgain = (field->gain > 0.0) ? 1.0/field->gain : 0.0;
  satlevel = field->satur_level - profit->obj->bkg;
  rad2 = h/2.0;
  if (rad2 > w/2.0)
    rad2 = w/2.0;
  rad2 *= rad2;


/* Set the image boundaries */
  pixout = profit->objpix;
  wpixout = profit->objweight;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<field->ymin)
    {
    off = (field->ymin-ymin)*w;
    pixout += off;
    wpixout += off;
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
    wpixout -= xmin;
    w2 += xmin;
    xmin = 0;
    }

/* Copy the right pixels to the destination */
  dw = w - w2;
  npix = 0;
  if (wfield)
    {
    wthresh = wfield->weight_thresh;
    gainflag = prefs.weightgain_flag;
/*-- Do the same for the weights */
    npix = 0;

    for (y=ymin; y<ymax; y++, pixout+=dw,wpixout+=dw)
      {
      dy2 = y-iy;
      dy2 *= dy2;
      pixin = &PIX(field, xmin, y);
      wpixin = &PIX(wfield, xmin, y);
      for (x=xmin; x<xmax; x++)
        {
        dx2 = (x-ix);
        dr2 = dy2 + dx2*dx2;
        pix = *(pixout++) = *(pixin++);
        if (dr2<rad2 && pix>-BIG && pix<satlevel && (wpix=*(wpixin++))<wthresh)
          {
          *(wpixout++) = 1.0 / sqrt(wpix+(pix>0.0?
		(gainflag? pix*wpix/backnoise2:pix)*invgain : 0.0));
          npix++;
          }
        else
          *(wpixout++) = 0.0;
        }
      }
    }
  else
    for (y=ymin; y<ymax; y++, pixout+=dw,wpixout+=dw)
      {
      dy2 = y-iy;
      dy2 *= dy2;
      pixin = &PIX(field, xmin, y);
      for (x=xmin; x<xmax; x++)
        {
        dx2 = (x-ix);
        dr2 = dy2 + dx2*dx2;
        pix = *(pixout++) = *(pixin++);
        if (dr2<rad2 && pix>-BIG)
          {
          *(wpixout++) = 1.0 / sqrt(backnoise2 + (pix>0.0?pix*invgain : 0.0));
          npix++;
          }
        else
          *(wpixout++) = 0.0;
        }
      }
 
  return npix;
  }


/****** profit_spiralindex ****************************************************
PROTO	double profit_spiralindex(profitstruct *profit)
PURPOSE	Compute the spiral index of a galaxy image (positive for arms
	extending counter-clockwise and negative for arms extending CW, 0 for
	no spiral pattern).
INPUT	Profile-fitting structure.
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/09/2008
 ***/
double profit_spiralindex(profitstruct *profit)
  {
   objstruct	*obj;
   obj2struct	*obj2;
   double	*dx,*dy, *fdx,*fdy, *gdx,*gdy, *gdxt,*gdyt, *pix,
		fwhm, invtwosigma2, hw,hh, ohw,ohh, x,y,xstart, tx,ty,txstart,
		gx,gy, r2, spirindex, invsig, val, sep;
   PIXTYPE	*fpix;
   int		i,j, npix;

  npix = profit->objnaxisn[0]*profit->objnaxisn[1];

  obj = profit->obj;
  obj2 = profit->obj2;
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
  gdy = NULL;			/* to avoid gcc -Wall warnings */
  QMEMCPY(gdx, gdy, double, npix);
  fdx = fft_rtf(dx, profit->objnaxisn);
  fft_conv(gdx, fdx, profit->objnaxisn);
  fdy = fft_rtf(dy, profit->objnaxisn);
  fft_conv(gdy, fdy, profit->objnaxisn);

/* Compute estimator */
  invtwosigma2 = -1.18*1.18 / (2.0*obj2->hl_radius*obj2->hl_radius);
  xstart = -hw - obj->mx + (int)(obj->mx+0.49999);
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


/****** profit_moments ****************************************************
PROTO	void profit_moments(profitstruct *profit)
PURPOSE	Compute the 2nd order moments from the unconvolved object model.
INPUT	Profile-fitting structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/09/2008
 ***/
void	 profit_moments(profitstruct *profit)
  {
   objstruct	*obj;
   obj2struct	*obj2;
   double	*pix,
		hw,hh, x,y, xstart, val,
		mx,my, sum, mx2,my2,mxy, den;
   int		ix,iy;

  obj = profit->obj;
  obj2 = profit->obj2;
  hw = (double)(profit->modnaxisn[0]/2);
  hh = (double)(profit->modnaxisn[1]/2);
  xstart = -hw;
  y = -hh;
  pix = profit->modpix;
  mx2 = my2 = mxy = mx = my = sum = 0.0;
  for (iy=profit->modnaxisn[1]; iy--; y+=1.0)
    {
    x = xstart;
    for (ix=profit->modnaxisn[0]; ix--; x+=1.0)
      {
      val = *(pix++);
      sum += val;
      mx  += val*x;
      my  += val*y;
      mx2 += val*x*x;
      mxy += val*x*y;
      my2 += val*y*y;
      }
    }

  if (sum <= 1.0/BIG)
    sum = 1.0;
  mx /= sum;
  my /= sum;
  obj2->prof_mx2 = mx2 = mx2/sum - mx*mx;
  obj2->prof_my2 = my2 = my2/sum - my*my;
  obj2->prof_mxy = mxy = mxy/sum - mx*my;
  if (mx2+my2 > 1.0/BIG)
    {
    obj2->prof_eps1 = (mx2 - my2) / (mx2+my2);
    obj2->prof_eps2 = 2.0*mxy / (mx2 + my2);
    den = mx2+my2-mxy*mxy;
    if (den>=0.0)
      den = mx2+my2+2.0*sqrt(den);
    else
      den = mx2+my2;
    obj2->prof_e1 = (mx2 - my2) / den;
    obj2->prof_e2 = 2.0*mxy / den;
    }
  else
    obj2->prof_eps1 = obj2->prof_eps2 = obj2->prof_e1 = obj2->prof_e2 = 0.0;

  return;
  }


/****** profit_addparam *******************************************************
PROTO	void profit_addparam(profitstruct *profit, paramenum paramindex,
		double **param)
PURPOSE	Add a profile parameter to the list of fitted items.
INPUT	Pointer to the profit structure,
	Parameter index,
	Pointer to the parameter pointer.
OUTPUT	-.
NOTES	-.
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
    {
    *param = profit->paramlist[paramindex] = &profit->param[profit->nparam];
    profit->paramindex[paramindex] = profit->nparam++;
    }

  return;
  }


/****** profit_resetparam ****************************************************
PROTO	void profit_resetparam(profitstruct *profit, paramenum paramtype)
PURPOSE	Set the initial, lower and upper boundary values of a profile parameter.
INPUT	Pointer to the profit structure,
	Parameter index.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2008
 ***/
void	profit_resetparam(profitstruct *profit, paramenum paramtype)
  {
   objstruct	*obj;
   obj2struct	*obj2;
   double	param, parammin,parammax;

  obj = profit->obj;
  obj2 = profit->obj2;
  param = parammin = parammax = 0.0;	/* Avoid gcc -Wall warnings*/
  switch(paramtype)
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
    case PARAM_SPHEROID_FLUX:
      param = obj2->flux_auto/2.0;
      parammin = -obj2->flux_auto/1000.0;
      parammax = 2*obj2->flux_auto;
      break;
    case PARAM_SPHEROID_REFF:
      param = obj2->hl_radius;
      parammin = 0.1;
      parammax = param * 4.0;
      break;
    case PARAM_SPHEROID_ASPECT:
      param = FLAG(obj2.prof_disk_flux)? 1.0 : obj->b/obj->a;
      parammin = FLAG(obj2.prof_disk_flux)? 0.5 : 0.01;
      parammax = 1.0;
      break;
    case PARAM_SPHEROID_POSANG:
      param = obj->theta;
      parammin = 0.0;
      parammax =  0.0;
      break;
    case PARAM_SPHEROID_SERSICN:
      param = 4.0;
      parammin = 1.0;
      parammax = 10.0;
      break;
    case PARAM_DISK_FLUX:
      param = obj2->flux_auto/2.0;
      parammin = -obj2->flux_auto/1000.0;
      parammax = 2*obj2->flux_auto;
      break;
    case PARAM_DISK_SCALE:	/* From scalelength to Re */
      param = obj2->hl_radius/1.67835*sqrt(obj->a/obj->b);
      parammin = param / 4.0;
      parammax = param * 4.0;
      break;
    case PARAM_DISK_ASPECT:
      param = obj->b/obj->a;
      parammin = 0.01;
      parammax = 1.0;
      break;
    case PARAM_DISK_POSANG:
      param = obj->theta;
      parammin = 0.0;
      parammax =  0.0;
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
      param = 0.5;
      parammin = 0.0;
      parammax = 3.0;
      break;
    case PARAM_ARMS_PITCH:
      param = 20.0;
      parammin = 5.0;
      parammax = 50.0;
      break;
    case PARAM_ARMS_PITCHVAR:
      param = 0.0;
      parammin = -1.0;
      parammax = 1.0;
      break;
//      if ((profit->spirindex=profit_spiralindex(profit, obj, obj2)) > 0.0)
//        {
//        param = -param;
//        parammin = -parammax;
//        parammax = -parammin;
//        }
//      printf("spiral index: %g  \n", profit->spirindex);
//      break;
    case PARAM_ARMS_POSANG:
      param = 0.0;
      parammin = 0.0;
      parammax = 0.0;
      break;
    case PARAM_ARMS_WIDTH:
      param = 3.0;
      parammin = 1.5;
      parammax = 11.0;
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
      parammin = 0.0;
      parammax = 0.0;
      break;
    case PARAM_INRING_FLUX:
      param = obj2->flux_auto/10.0;
      parammin = 0.0;
      parammax = 2.0*obj2->flux_auto;
      break;
    case PARAM_INRING_WIDTH:
      param = 0.3;
      parammin = 0.0;
      parammax = 0.5;
      break;
    case PARAM_INRING_ASPECT:
      param = 0.8;
      parammin = 0.4;
      parammax = 1.0;
      break;
    case PARAM_OUTRING_FLUX:
      param = obj2->flux_auto/10.0;
      parammin = 0.0;
      parammax = 2.0*obj2->flux_auto;
      break;
    case PARAM_OUTRING_START:
      param = 4.0;
      parammin = 3.5;
      parammax = 6.0;
      break;
    case PARAM_OUTRING_WIDTH:
      param = 0.3;
      parammin = 0.0;
      parammax = 0.5;
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown profile parameter in ",
		"profit_resetparam()");
      break;
   }

  if (parammin!=parammax && (param<=parammin || param>=parammax))
    param = (parammin+parammax)/2.0;
  profit_setparam(profit, paramtype, param, parammin, parammax);

  return;
  }


/****** profit_resetparams ****************************************************
PROTO	void profit_resetparams(profitstruct *profit)
PURPOSE	Set the initial, lower and upper boundary values of profile parameters.
INPUT	Pointer to the profit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/09/2008
 ***/
void	profit_resetparams(profitstruct *profit)
  {
   int		p;


  for (p=0; p<PARAM_NPARAM; p++)
    profit_resetparam(profit, (paramenum)p);

  return;
  }


/****** profit_setparam ****************************************************
PROTO	void profit_setparam(profitstruct *profit, paramenum paramtype,
		double param, double parammin, double parammax)
PURPOSE	Set the actual, lower and upper boundary values of a profile parameter.
INPUT	Pointer to the profit structure,
	Parameter index,
	Actual value,
	Lower boundary to the parameter,
	Upper boundary to the parameter.
OUTPUT	RETURN_OK if the parameter is registered, RETURN_ERROR otherwise.
AUTHOR	E. Bertin (IAP)
VERSION	15/03/2009
 ***/
int	profit_setparam(profitstruct *profit, paramenum paramtype,
		double param, double parammin, double parammax)
  {
   double	*paramptr;
   int		index;

/* Check whether the parameter has already be registered */
  if ((paramptr=profit->paramlist[(int)paramtype]))
    {
    index = profit->paramindex[(int)paramtype];
    profit->paraminit[index] = param;
    profit->parammin[index] = parammin;
    profit->parammax[index] = parammax;
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }

  
/****** profit_boundtounbound *************************************************
PROTO	void profit_boundtounbound(profitstruct *profit, double *param)
PURPOSE	Convert parameters from bounded to unbounded space.
INPUT	Pointer to the profit structure.
OUTPUT	-.
NOTES	-.
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
PROTO	void profit_unboundtobound(profitstruct *profit, double *param)
PURPOSE	Convert parameters from unbounded to bounded space.
INPUT	Pointer to the profit structure.
OUTPUT	-.
NOTES	-.
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


/****** profit_covarunboundtobound ********************************************
PROTO	void profit_covarunboundtobound(profitstruct *profit)
PURPOSE	Convert covariance matrix from unbounded to bounded space.
INPUT	Pointer to the profit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	30/04/2008
 ***/
void	profit_covarunboundtobound(profitstruct *profit)
  {
   double	*covar, *dxdy,*x,*xmin,*xmax,
		dxmin, dxmax;
   int		p,p1,p2, nparam;

  nparam = profit->nparam;
  QMALLOC(dxdy, double, nparam);
  x = profit->paraminit;
  xmin = profit->parammin;
  xmax = profit->parammax;
  for (p=0; p<profit->nparam; p++)
    if (xmin[p]!=xmax[p])
      {
      dxmin = x[p] - xmin[p];
      dxmax= xmax[p] - x[p];
      dxdy[p] = (fabs(dxmin) < 1.0/BIG && fabs(dxmax) < 1.0/BIG) ?
		0.0 : dxmin*dxmax/(dxmin+dxmax);
      }
    else
      dxdy[p] = 1.0;

  covar = profit->covar;
  for (p2=0; p2<nparam; p2++)
    for (p1=0; p1<nparam; p1++)
      *(covar++) *= dxdy[p1]*dxdy[p2];

  free(dxdy);

  return;
  }


/****** prof_init *************************************************************
PROTO	profstruct prof_init(profitstruct *profit, proftypenum profcode)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	Pointer to the profile-fitting structure,
	profile type.
OUTPUT	A pointer to an allocated prof structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	22/04/2008
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
      profit_addparam(profit, PARAM_SPHEROID_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_SPHEROID_REFF, &prof->scale);
      profit_addparam(profit, PARAM_SPHEROID_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SPHEROID_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_SPHEROID_SERSICN, &prof->extra[0]);
      break;
    case PROF_DEVAUCOULEURS:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_SPHEROID_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_SPHEROID_REFF, &prof->scale);
      profit_addparam(profit, PARAM_SPHEROID_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SPHEROID_POSANG, &prof->posangle);
      break;
    case PROF_EXPONENTIAL:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      break;
    case PROF_ARMS:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_ARMS_QUADFRAC, &prof->featfrac);
//      profit_addparam(profit, PARAM_ARMS_SCALE, &prof->featscale);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_ARMS_PITCH, &prof->featpitch);
//      profit_addparam(profit, PARAM_ARMS_PITCHVAR, &prof->featpitchvar);
      profit_addparam(profit, PARAM_ARMS_POSANG, &prof->featposang);
//      profit_addparam(profit, PARAM_ARMS_WIDTH, &prof->featwidth);
      break;
    case PROF_BAR:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
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
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_INRING_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_INRING_WIDTH, &prof->featwidth);
      profit_addparam(profit, PARAM_INRING_ASPECT, &prof->feataspect);
      break;
    case PROF_OUTRING:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_OUTRING_START, &prof->featstart);
      profit_addparam(profit, PARAM_OUTRING_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_OUTRING_WIDTH, &prof->featwidth);
      break;
    case PROF_DIRAC:
      prof->naxis = 2;
      prof->pix = NULL;
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_FLUX, &prof->flux);
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
      profit_addparam(profit, PARAM_SPHEROID_FLUX, &prof->flux);
      profit_addparam(profit, PARAM_SPHEROID_REFF, &prof->scale);
      profit_addparam(profit, PARAM_SPHEROID_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SPHEROID_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_SPHEROID_SERSICN, &prof->extra[0]);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown profile in ",
		"prof_init()");
      break;
    }

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
NOTES	-.
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
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/10/2008
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
   double	posin[PROFIT_MAXEXTRA], posout[2], dnaxisn[2],
		*pixout,
		flux,fluxfac, scaling;
   float	*pixin,
		amp,ctheta,stheta,cd11,cd12,cd21,cd22, dcd11,dcd21, dx1,dx2,
		x1,x10,x2, x1cin,x2cin, x1cout,x2cout, xscale,yscale, saspect,
		x1in,x2in, odx, ostep,
		n,k, hinvn, x1t,x2t, ca,sa, u,umin,
		armamp,arm2amp, armrdphidr, armrdphidrvar, posang,
		width, invwidth2,
		r,r2,rmin,r2min, r2minxin,r2minxout, rmax, r2max, invr2xdif,
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

  scaling = profit->pixstep / prof->typscale;

  if (prof->code!=PROF_DIRAC)
    {
/*-- Compute Profile CD matrix */
    ctheta = cos(*prof->posangle*DEG);
    stheta = sin(*prof->posangle*DEG);
    saspect = fabs(*prof->aspect);
    xscale = (*prof->scale==0.0)?
		0.0 : fabs(scaling / (*prof->scale*prof->typscale));
    yscale = (*prof->scale*saspect == 0.0)?
		0.0 : fabs(scaling / (*prof->scale*prof->typscale*saspect));
    cd11 = xscale*ctheta;
    cd12 = xscale*stheta;
    cd21 = -yscale*stheta;
    cd22 = yscale*ctheta;
    }

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
      x10 = -x1cout - dx1;
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
            dcd11 = cd11*ostep;
            dcd21 = cd21*ostep;
            odx = 0.5*(ostep-1.0);
            x1t = x1+odx;
            val = 0.0;
            for (idx2=noversamp; idx2--; odx+=ostep)
              {
              x1in = cd12*(x2+odx) + cd11*x1t;
              x2in = cd22*(x2+odx) + cd21*x1t;
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
      x10 = -x1cout - dx1;
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
          val = expf(-7.6692f*PROFIT_POWF(ra,0.125));
          noversamp  = (int)(sqrt(val)*PROFIT_OVERSAMP+0.1);
          if (noversamp < 2)
            *(pixin++) = val;
          else
            {
            ostep = 1.0/noversamp;
            dcd11 = cd11*ostep;
            dcd21 = cd21*ostep;
            odx = 0.5*(ostep-1.0);
            x1t = x1+odx;
            val = 0.0;
            for (idx2=noversamp; idx2--; odx+=ostep)
              {
              x1in = cd12*(x2+odx) + cd11*x1t;
              x2in = cd22*(x2+odx) + cd21*x1t;
              for (idx1=noversamp; idx1--;)
                {
                ra = x1in*x1in+x2in*x2in;
                val += expf(-7.6692f*PROFIT_POWF(ra,0.125));
                x1in += dcd11;
                x2in += dcd21;
                }
              }
            *(pixin++) = val*ostep*ostep;
            }
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
      umin = 0.5*logf(r2minxin + 0.00001);
      arm2amp = *prof->featfrac;
      armamp = 1.0 - arm2amp;
      armrdphidr = 1.0/tan(*prof->featpitch*DEG);
      armrdphidrvar = 0.0 /**prof->featpitchvar*/;
      posang = *prof->featposang*DEG;
      width = fabs(*prof->featwidth);
width = 3.0;
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
            u = 0.5*logf(r2 + 0.00001);
            theta = (armrdphidr+armrdphidrvar*(u-umin))*u+posang;
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
      cd22 /= *prof->feataspect;
      cd21 /= *prof->feataspect;
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
    case PROF_OUTRING:
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
    case PROF_DIRAC:
      memset(profit->pmodpix, 0,
	profit->modnaxisn[0]*profit->modnaxisn[1]*sizeof(float));
      profit->pmodpix[profit->modnaxisn[0]/2
		+ (profit->modnaxisn[1]/2)*profit->modnaxisn[0]] = 1.0;
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
  rmax = x1cout;
  if (rmax > (r = x2cout))
    rmax = r;
  rmax += 0.01;
  if (rmax<1.0)
    rmax = 1.0;
  r2max = rmax*rmax;
  rmin = rmax - 1.0;
  r2min = rmin*rmin;

/* Find best threshold (the max around the circle with radius rmax */
  dx2 = -x2cout;
  pixin = profit->pmodpix;
  thresh = -BIG;
  for (ix2=profit->modnaxisn[1]; ix2--; dx2 += 1.0)
    {
    dx1 = -x1cout;
    for (ix1=profit->modnaxisn[0]; ix1--; dx1 += 1.0)
      if ((val=*(pixin++))>thresh && (r2=dx1*dx1+dx2*dx2)>r2min && r2<r2max)
        thresh = val;
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
  prof->fluxfac = fluxfac;
  pixin = profit->pmodpix;
  pixout = profit->modpix;
  for (n=npix; n--;)
    *(pixout++) += fluxfac * *(pixin++);

  return;
  }


/****** prof_interpolate ******************************************************
PROTO	double	prof_interpolate(profstruct *prof, double *posin)
PURPOSE	Interpolate a multidimensional model profile at a given position.
INPUT	Profile structure,
	input position vector.
OUTPUT	-.
NOTES	-.
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
NOTES	-.
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
PROTO	void make_kernel(double pos, double *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/04/2008
 ***/
void	make_kernel(double pos, double *kernel, interpenum interptype)
  {
   double	x, val, sinx1,sinx2,sinx3,cosx1;

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
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0;
      val += (*kernel = cosx1/(x*x));
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
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx2=-0.5*sinx1-0.866025403785*cosx1)
				/ (x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx3=-0.5*sinx1+0.866025403785*cosx1)
				/(x*x));
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
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/4.0;
      val +=(*(kernel++) = -(sinx2=0.707106781186*(sinx1+cosx1))
				/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = cosx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -(sinx3=0.707106781186*(cosx1-sinx1))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/4.0;
      val += (*kernel = sinx3/(x*x));
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

