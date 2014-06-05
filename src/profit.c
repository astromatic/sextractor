/*
*				profit.c
*
* Fit a range of galaxy models to an image.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2006-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		04/12/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifndef HAVE_MATHIMF_H
#define _GNU_SOURCE
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"levmar/levmar.h"
#include	"fft.h"
#include	"fitswcs.h"
#include	"check.h"
#include	"graph.h"
#include	"image.h"
#include	"misc.h"
#include	"pattern.h"
#include	"psf.h"
#include	"profit.h"
#include	"subimage.h"

static double	prof_gammainc(double x, double a),
		prof_gamma(double x);
static float	prof_interpolate(profstruct *prof, float *posin);
static float	interpolate_pix(float *posin, float *pix, int *naxisn,
		interpenum interptype);

static void	make_kernel(float pos, float *kernel, interpenum interptype);

/*------------------------------- variables ---------------------------------*/

const int	interp_kernwidth[5]={1,2,4,6,8};

const int	flux_flag[PARAM_NPARAM] = {0,
					1,2,3,4,5,0,0,
					1,2,3,4,5,0,0,0,0,
					1,2,3,4,5,0,0,0,
					1,0,0,0,0,0,0,0,
					1,0,0,
					1,0,0,
					1,0,0
					};

int	the_gal;

/****** profit_init *******************************************************//**
Allocate and initialize a new model-fitting structure.
@param[in] obj		Pointer to an object
@param[in] subimage	Pointer to an array of subimages (or a single subimage)
@param[in] nsubimage	Number of subimages
@param[in] modeltype	Binary mask of model types
@param[in] conv_flag	Convolution flag (0 = no convolution)

@author 		E. Bertin (IAP)
@date			28/03/2014
 ***/
profitstruct	*profit_init(objstruct *obj, subimagestruct *subimage,
			int nsubimage, unsigned int modeltype, int conv_flag)
  {
   profitstruct		*profit;
   subprofitstruct	*subprofit;
   obj2struct		*obj2;
   double		psf_fwhm;
   int			s,t, nmodels, npix, nsub;

  QCALLOC(profit, profitstruct, 1);

  obj2 = obj->obj2;
  profit->conv_flag = conv_flag;

  QCALLOC(profit->subprofit, subprofitstruct, nsubimage);
  subprofit = profit->subprofit;
  nsub = 0;
  for (s=0; s<nsubimage; s++, subimage++)
    {
    subprofit->field = subimage->field;
    subprofit->wfield = subimage->wfield;
    if (conv_flag)
      {
      subprofit->psf = subimage->field->psf;
      subprofit->pixstep = subprofit->psf->pixstep;
      }
    else
      {
      subprofit->psf = NULL;
      subprofit->pixstep = 1.0;
      }
    subprofit->sigma = obj2? obj2->sigbkg[s] : obj->sigbkg;
    subprofit->fluxfac = 1.0;	/* Default */

/*-- Set initial guesses and boundaries */
    subprofit->guesssigbkg = subprofit->sigma;
    subprofit->guessdx = obj->mx - (int)(obj->mx+0.49999);
    subprofit->guessdy = obj->my - (int)(obj->my+0.49999);
    subprofit->guessradius = obj2? obj2->hl_radius : 1.2 * sqrtf(obj->a*obj->b);
					// 1.2 is ~2.35/2
    if (conv_flag && subprofit->guessradius < 0.5*subprofit->psf->fwhm)
      subprofit->guessradius = 0.5*subprofit->psf->fwhm;
    if ((subprofit->guessflux = obj2? obj2->flux_auto[s] : obj->dflux) <= 0.0)
      subprofit->guessflux = 0.0;
    if ((subprofit->guessfluxmax = 10.0
		* (obj2? obj2->fluxerr_auto[s] : obj->dfluxerr))
	<= subprofit->guessflux)
      subprofit->guessfluxmax = subprofit->guessflux;
    if (subprofit->guessfluxmax <= 0.0)
      subprofit->guessfluxmax = 1.0;
    subprofit->guessaspect = obj->b/obj->a;
    subprofit->guessposang = obj->theta;

/*-- Create pixmaps at image resolution */
    subprofit->ix = (int)(obj->mx + 0.49999); /* 1st pix=0 */
    subprofit->iy = (int)(obj->my + 0.49999); /* 1st pix=0 */
    psf_fwhm = conv_flag? subprofit->psf->masksize[0]*subprofit->psf->pixstep
			: 0.0;
    subprofit->objnaxisn[0] = ((int)((subimage->size[0] + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
    subprofit->objnaxisn[1] = ((int)((subimage->size[1] + psf_fwhm + 0.499)
		*1.2)/2)*2 + 1;
    if (subprofit->objnaxisn[1] < subprofit->objnaxisn[0])
      subprofit->objnaxisn[1] = subprofit->objnaxisn[0];
    else
      subprofit->objnaxisn[0] = subprofit->objnaxisn[1];

    if (subprofit->objnaxisn[0]>PROFIT_MAXOBJSIZE)
      {
      subprofit->subsamp
		= ceil((float)subprofit->objnaxisn[0]/PROFIT_MAXOBJSIZE);
      subprofit->objnaxisn[1]
		= (subprofit->objnaxisn[0] /= (int)subprofit->subsamp);
      profit->flag |= PROFLAG_OBJSUB;	/* !CHECK */
      }
    else
      subprofit->subsamp = 1.0;
    subprofit->nobjpix = subprofit->objnaxisn[0]*subprofit->objnaxisn[1];

/* Allocate memory for the data and the resampled model */
    QMALLOC16(subprofit->objpix, PIXTYPE, subprofit->nobjpix);
    QMALLOC16(subprofit->objweight, PIXTYPE, subprofit->nobjpix);
    QMALLOC16(subprofit->lmodpix, PIXTYPE, subprofit->nobjpix);
    QMALLOC16(subprofit->lmodpix2, PIXTYPE, subprofit->nobjpix);

/*-- Create pixmap at model resolution */
    subprofit->modnaxisn[0] =
	((int)(subprofit->objnaxisn[0]*subprofit->subsamp
		/subprofit->pixstep+0.4999)/2+1)*2;
    subprofit->modnaxisn[1] =
	((int)(subprofit->objnaxisn[1]*subprofit->subsamp
		/subprofit->pixstep+0.4999)/2+1)*2;
    if (subprofit->modnaxisn[1] < subprofit->modnaxisn[0])
      subprofit->modnaxisn[1] = subprofit->modnaxisn[0];
    else
      subprofit->modnaxisn[0] = subprofit->modnaxisn[1];
    if (subprofit->modnaxisn[0]>PROFIT_MAXMODSIZE)
      {
      subprofit->pixstep = (double)subprofit->modnaxisn[0] / PROFIT_MAXMODSIZE;
      subprofit->modnaxisn[0] = subprofit->modnaxisn[1] = PROFIT_MAXMODSIZE;
      profit->flag |= PROFLAG_MODSUB;	/* !CHECK */
      }
    subprofit->nmodpix = subprofit->modnaxisn[0]*subprofit->modnaxisn[1];

/* Allocate memory for the complete model */
    QFFTWF_MALLOC(subprofit->modpix, float, subprofit->nmodpix);
    QMALLOC16(subprofit->modpix2, float, subprofit->nmodpix);
    QMALLOC16(subprofit->cmodpix, float, subprofit->nmodpix);
    if (conv_flag)
      QMALLOC16(subprofit->psfpix, float, subprofit->nmodpix);
    if ((npix = subprofit_copyobjpix(subprofit, subimage)))
      {
      profit->nresi += npix;
/*---- Compute the local PSF */
      if (conv_flag)
        subprofit_psf(subprofit, obj2);
      subprofit->index = nsub++;
      subprofit++;
      }
    else
      subprofit_end(subprofit);
    }

  profit->nsubprofit = nsub;

  QMALLOC16(profit->resi, float, profit->nresi);

  QMALLOC(profit->prof, profstruct *, MODEL_NMAX);
  nmodels = 0;
  for (t=1; t<(1<<MODEL_NMAX); t<<=1)
    if (modeltype&t)
      profit->prof[nmodels++] = prof_init(profit, t);
  profit->nprof = nmodels;

  QMALLOC16(profit->covar, float, profit->nparam*profit->nparam);

  profit_resetparams(profit);

  return profit;
  }  


/****** profit_end ************************************************************
PROTO	void profit_end(profitstruct *profit)
PURPOSE	End (deallocate) a profile-fitting structure.
INPUT	Prof structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
void	profit_end(profitstruct *profit)
  {
   subprofitstruct	*subprofit;
   int			p,s;

  for (p=0; p<profit->nprof; p++)
    prof_end(profit->prof[p]);
  free(profit->prof);
  subprofit = profit->subprofit;
  for (s=profit->nsubprofit; s--;)
    subprofit_end(subprofit++);

  free(profit->subprofit);
  free(profit->resi);
  free(profit->covar);

  free(profit);

  return;
  }


/****** subprofit_end *********************************************************
PROTO	void subprofit_end(subprofitstruct *profit)
PURPOSE	End (free content) a subprofile-fitting structure.
INPUT	SubProfit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/08/2012
 ***/
void	subprofit_end(subprofitstruct *subprofit)
  {
  free(subprofit->objpix);
  free(subprofit->objweight);
  free(subprofit->lmodpix);
  free(subprofit->lmodpix2);
  QFFTWF_FREE(subprofit->modpix);
  free(subprofit->modpix2);
  free(subprofit->cmodpix);
  free(subprofit->psfpix);
  QFFTWF_FREE(subprofit->psfdft);
  fft_scratchend(subprofit->fftscratch);
  subprofit->fftscratch = NULL;

  return;
  }


/****** profit_fit ************************************************************
PROTO	int profit_fit(profitstruct *profit, objstruct *obj)
PURPOSE	Fit profile(s) convolved with the PSF to a detected object.
INPUT	Pointer to the model structure,
	pointer to the obj.
OUTPUT	Number of minimization iterations (0 in case of error).
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2012
 ***/
int	profit_fit(profitstruct *profit, objstruct *obj)
  {
   int			p, nparam, nparam2;
    checkstruct		*check;

  nparam = profit->nparam;
  nparam2 = nparam*nparam;

/* reset all flags except those set in profit_init() */
  profit->flag &= PROFLAG_OBJSUB|PROFLAG_MODSUB;

/* Check if the number of constraints exceeds the number of free parameters */
  if (profit->nresi < nparam)
    {
    if (obj->obj2) {
      if (FLAG(obj2.prof_vector))
        for (p=0; p<nparam; p++)
          obj2->prof_vector[p] = 0.0;
      if (FLAG(obj2.prof_errvector))
        for (p=0; p<nparam; p++)
          obj2->prof_errvector[p] = 0.0;
      if (FLAG(obj2.prof_errmatrix))
        for (p=0; p<nparam2; p++)
          obj2->prof_errmatrix[p] = 0.0;
      obj2->prof_niter = 0;
    }
    profit->flag |= PROFLAG_NOTCONST;
    return 0;
    }

/* Actual minimisation */
  profit->niter = profit_minimize(profit, PROFIT_MAXITER);

  return profit->niter;
  }


/****** profit_measure ********************************************************
PROTO	void profit_measure(profitstruct *profit, obj2struct *obj2)
PURPOSE	Perform measurements on the fitted model.
INPUT	Pointer to the model structure,
	pointer to the obj2.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/12/2013
 ***/
void	profit_measure(profitstruct *profit, obj2struct *obj2)
  {
    profitstruct	*pprofit, *qprofit;
    subprofitstruct	*subprofit;
    patternstruct	*pattern;
    checkstruct		*check;
    double		emx2,emy2,emxy, a , cp,sp, cn, bn, n;
    float		param0[PARAM_NPARAM], param1[PARAM_NPARAM],
			param[PARAM_NPARAM],
			**list,
			*cov,
			err, aspect, chi2, flux;
    int			*index,
			i,j,p,s, nparam, nparam2, ncomp, nsub;

  nparam = profit->nparam;
  nparam2 = nparam*nparam;
  nsub = profit->nsubprofit;
  subprofit = profit->subprofit;

  if (profit->nlimmin)
    profit->flag |= PROFLAG_MINLIM;
  if (profit->nlimmax)
    profit->flag |= PROFLAG_MAXLIM;

  for (p=0; p<nparam; p++)
    profit->paramerr[p]= sqrt(profit->covar[p*(nparam+1)]);

/* CHECK-Images */
  if ((check = prefs.check[CHECK_PROFILES]))
    {
    profit_residuals(profit, 0.0, profit->paraminit, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0);
    }

  if ((check = prefs.check[CHECK_SUBPROFILES]))
    {
    profit_residuals(profit, 0.0, profit->paraminit, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		-1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, -1.0);
    }
return;
  if ((check = prefs.check[CHECK_SPHEROIDS]))
    {
/*-- Set to 0 flux components that do not belong to spheroids */
    for (p=0; p<profit->nparam; p++)
      param[p] = profit->paraminit[p];
    list = profit->paramlist;
    index = profit->paramindex;
    for (i=0; i<PARAM_NPARAM; i++)
      if (list[i] && flux_flag[i] && i!= PARAM_SPHEROID_FLUX)
        param[index[i]] = 0.0;
    profit_residuals(profit, 0.0, param, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0);
    }
  if ((check = prefs.check[CHECK_SUBSPHEROIDS]))
    {
/*-- Set to 0 flux components that do not belong to spheroids */
    for (p=0; p<profit->nparam; p++)
      param[p] = profit->paraminit[p];
    list = profit->paramlist;
    index = profit->paramindex;
    for (i=0; i<PARAM_NPARAM; i++)
      if (list[i] && flux_flag[i] && i!= PARAM_SPHEROID_FLUX)
        param[index[i]] = 0.0;
    profit_residuals(profit, 0.0, param, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		-1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, -1.0);
    }
  if ((check = prefs.check[CHECK_DISKS]))
    {
/*-- Set to 0 flux components that do not belong to disks */
    for (p=0; p<profit->nparam; p++)
      param[p] = profit->paraminit[p];
    list = profit->paramlist;
    index = profit->paramindex;
    for (i=0; i<PARAM_NPARAM; i++)
      if (list[i] && flux_flag[i] && i!= PARAM_DISK_FLUX)
        param[index[i]] = 0.0;
    profit_residuals(profit, 0.0, param, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0);
    }
  if ((check = prefs.check[CHECK_SUBDISKS]))
    {
/*-- Set to 0 flux components that do not belong to disks */
    for (p=0; p<profit->nparam; p++)
      param[p] = profit->paraminit[p];
    list = profit->paramlist;
    index = profit->paramindex;
    for (i=0; i<PARAM_NPARAM; i++)
      if (list[i] && flux_flag[i] && i!= PARAM_DISK_FLUX)
        param[index[i]] = 0.0;
    profit_residuals(profit, 0.0, param, NULL);
    if (subprofit->subsamp>1.0)
      check_addresample(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, 1.0/subprofit->subsamp,
		-1.0/(subprofit->subsamp*subprofit->subsamp));
    else
      check_add(check, subprofit->lmodpix,
		subprofit->objnaxisn[0],subprofit->objnaxisn[1],
		subprofit->ix,subprofit->iy, -1.0);
    }
 
/* Compute compressed residuals */
  profit_residuals(profit, 10.0, profit->paraminit,profit->resi);

/* Fill measurement parameters */
  if (FLAG(obj2.prof_vector))
    {
    for (p=0; p<nparam; p++)
      obj2->prof_vector[p]= profit->paraminit[p];
    }
  if (FLAG(obj2.prof_errvector))
    {
    for (p=0; p<nparam; p++)
      obj2->prof_errvector[p]= profit->paramerr[p];
    }
  if (FLAG(obj2.prof_errmatrix))
    {
    for (p=0; p<nparam2; p++)
      obj2->prof_errmatrix[p]= profit->covar[p];
    }

  obj2->prof_niter = profit->niter;
  subprofit = profit->subprofit;
  for (s=0; s<nsub; s++, subprofit++)
    {
    if (FLAG(obj2.flux_prof))
      obj2->flux_prof[s] = subprofit->flux;
    if (FLAG(obj2.mag_prof))
      obj2->mag_prof[s] = subprofit->flux>0.0?
		-FDMAG * log(subprofit->flux) + prefs.mag_zeropoint[s] : 99.0f;
    if (FLAG(obj2.fluxerr_prof) || FLAG(obj2.magerr_prof))
      {
      err = 0.0;
      cov = profit->covar;
      index = profit->paramindex;
      list = profit->paramlist;
      for (i=0; i<PARAM_NPARAM; i++)
        if (flux_flag[i]==s+1 && list[i])
          {
          cov = profit->covar + nparam*index[i];
          for (j=0; j<PARAM_NPARAM; j++)
            if (flux_flag[j]==s+1 && list[j])
              err += cov[index[j]];
          }
      if (err<0.0)
        err = 0.0;
      if (FLAG(obj2.fluxerr_prof))
        obj2->fluxerr_prof[s] = sqrt(err);
      if (FLAG(obj2.magerr_prof))
        obj2->magerr_prof[s] = subprofit->flux>0.0?
		FDMAG * sqrt(err) / subprofit->flux : 99.0f;
      }
    }

  obj2->prof_chi2 = (profit->nresi > profit->nparam)?
		profit->chi2 / (profit->nresi - profit->nparam) : 0.0;

/* Position */
  if (FLAG(obj2.x_prof))
    {
    i = profit->paramindex[PARAM_X];
    j = profit->paramindex[PARAM_Y];
/*-- Model coordinates follow the FITS convention (first pixel at 1,1) */
    if (profit->paramlist[PARAM_X])
      {
      obj2->x_prof = (double)profit->subprofit->ix + *profit->paramlist[PARAM_X]
			+ 1.0;	/* !CHECK */	/* FITS convention */
      obj2->poserrmx2_prof = emx2 = profit->covar[i*(nparam+1)];
      }
    else
      emx2 = 0.0;
    if (profit->paramlist[PARAM_Y])
      {
      obj2->y_prof = (double)profit->subprofit->iy + *profit->paramlist[PARAM_Y]
			+ 1.0;	/* !CHECK */	/* FITS convention */
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
      obj2->poserrtheta_prof = (float)(theta/DEG);
      }

    if (FLAG(obj2.poserrcxx_prof))
      {
       double	temp;

      obj2->poserrcxx_prof = (float)(emy2/(temp=emx2*emy2-emxy*emxy));
      obj2->poserrcyy_prof = (float)(emx2/temp);
      obj2->poserrcxy_prof = (float)(-2*emxy/temp);
      }
    }

/* Equivalent noise area */
  if (FLAG(obj2.prof_noisearea))
    obj2->prof_noisearea = subprofit_noisearea(profit->subprofit); /* !CHECK */

/* Second order moments and ellipticities */
  if (FLAG(obj2.prof_mx2))
    profit_moments(profit, obj2);

/* Second order moments of the convolved model (used by other parameters) */
  if (FLAG(obj2.prof_convmx2))
    subprofit_convmoments(profit->subprofit, obj2); /* !CHECK */

/* "Hybrid" magnitudes */
  if (FLAG(obj2.fluxcor_prof))
    {
    profit_residuals(profit, 0.0, profit->paraminit, NULL);
    subprofit_fluxcor(profit->subprofit, obj2); /* !CHECK */
    }

/* Do measurements on the rasterised model (surface brightnesses) */
  if (FLAG(obj2.fluxeff_prof))
    subprofit_surface(profit, profit->subprofit, obj2);  /* !CHECK */

/* Background offset */
  if (FLAG(obj2.prof_offset_flux))
    {
    obj2->prof_offset_flux = *profit->paramlist[PARAM_BACK];
    obj2->prof_offset_fluxerr=profit->paramerr[profit->paramindex[PARAM_BACK]];
    }

/* Point source */
  if (FLAG(obj2.prof_dirac_flux))
    {
    for (s=0; s<nsub; s++)
      {
      flux = obj2->prof_dirac_flux[s] = *profit->paramlist[PARAM_DIRAC_FLUX+s];
      if (FLAG(obj2.prof_dirac_mag))
        obj2->prof_dirac_mag[s] = flux>0.0f?
		 -FDMAG * log(flux) + prefs.mag_zeropoint[s] : 99.0f;
      if (FLAG(obj2.prof_dirac_fluxerr))
        obj2->prof_dirac_fluxerr[s] =
		profit->paramerr[profit->paramindex[PARAM_DIRAC_FLUX+s]];
      if (FLAG(obj2.prof_dirac_magerr))
        obj2->prof_dirac_magerr[s] = flux>0.0f?
		 FDMAG*profit->paramerr[profit->paramindex[PARAM_DIRAC_FLUX+s]]
		/ flux : 99.0f;
      }
    }

/* Spheroid */
  if (FLAG(obj2.prof_spheroid_flux))
    {
    if ((aspect = *profit->paramlist[PARAM_SPHEROID_ASPECT]) > 1.0)
      {
      *profit->paramlist[PARAM_SPHEROID_REFF] *= aspect;
      profit->paramerr[profit->paramindex[PARAM_SPHEROID_REFF]] *= aspect;
      profit->paramerr[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			/= (aspect*aspect);
      *profit->paramlist[PARAM_SPHEROID_ASPECT] = 1.0 / aspect;
      *profit->paramlist[PARAM_SPHEROID_POSANG] += 90.0;
      }
    for (s=0; s<nsub; s++)
      {
      flux = obj2->prof_spheroid_flux[s]
		= *profit->paramlist[PARAM_SPHEROID_FLUX+s];
      if (FLAG(obj2.prof_spheroid_mag))
        obj2->prof_spheroid_mag[s] = flux>0.0f?
		 -FDMAG * log(flux) + prefs.mag_zeropoint[s] : 99.0f;
      if (FLAG(obj2.prof_spheroid_fluxerr))
        obj2->prof_spheroid_fluxerr[s] =
		profit->paramerr[profit->paramindex[PARAM_SPHEROID_FLUX+s]];
      if (FLAG(obj2.prof_spheroid_magerr))
        obj2->prof_spheroid_magerr[s] = flux>0.0f?
	      FDMAG*profit->paramerr[profit->paramindex[PARAM_SPHEROID_FLUX+s]]
		/ flux : 99.0f;
      }
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
    else
      obj2->prof_spheroid_sersicn = 4.0;
    if (FLAG(obj2.prof_spheroid_peak))
      {
      n = obj2->prof_spheroid_sersicn;
      bn = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
		+ 131.0/(1148175*n*n*n);	/* Ciotti & Bertin 1999 */
      cn = n * prof_gamma(2.0*n) * pow(bn, -2.0*n);
      obj2->prof_spheroid_peak = obj2->prof_spheroid_reff>0.0?
	obj2->prof_spheroid_flux[0]
		/ (2.0 * PI * cn
		* obj2->prof_spheroid_reff*obj2->prof_spheroid_reff
		* obj2->prof_spheroid_aspect)
	: 0.0;
      if (FLAG(obj2.prof_spheroid_fluxeff))
        obj2->prof_spheroid_fluxeff = obj2->prof_spheroid_peak * exp(-bn);
      if (FLAG(obj2.prof_spheroid_fluxmean))
        obj2->prof_spheroid_fluxmean = obj2->prof_spheroid_peak * cn;
      }
    }

/* Disk */
  if (FLAG(obj2.prof_disk_flux))
    {
    if ((aspect = *profit->paramlist[PARAM_DISK_ASPECT]) > 1.0)
      {
      *profit->paramlist[PARAM_DISK_SCALE] *= aspect;
      profit->paramerr[profit->paramindex[PARAM_DISK_SCALE]] *= aspect;
      profit->paramerr[profit->paramindex[PARAM_DISK_ASPECT]]
			/= (aspect*aspect);
      *profit->paramlist[PARAM_DISK_ASPECT] = 1.0 / aspect;
      *profit->paramlist[PARAM_DISK_POSANG] += 90.0;
      }
    for (s=0; s<nsub; s++)
      {
      flux = obj2->prof_disk_flux[s] = *profit->paramlist[PARAM_DISK_FLUX+s];
      if (FLAG(obj2.prof_disk_mag))
        obj2->prof_disk_mag[s] = flux>0.0f?
		 -FDMAG * log(flux) + prefs.mag_zeropoint[s] : 99.0f;
      if (FLAG(obj2.prof_disk_fluxerr))
        obj2->prof_disk_fluxerr[s] =
		profit->paramerr[profit->paramindex[PARAM_DISK_FLUX+s]];
      if (FLAG(obj2.prof_disk_magerr))
        obj2->prof_disk_magerr[s] = flux>0.0f?
		 FDMAG*profit->paramerr[profit->paramindex[PARAM_DISK_FLUX+s]]
		/ flux : 99.0f;
      }
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

    if (FLAG(obj2.prof_disk_peak))
      {
      obj2->prof_disk_peak = obj2->prof_disk_scale>0.0?
	obj2->prof_disk_flux[0]
	/ (2.0 * PI * obj2->prof_disk_scale*obj2->prof_disk_scale
		* obj2->prof_disk_aspect)
	: 0.0;
      if (FLAG(obj2.prof_disk_fluxeff))
        obj2->prof_disk_fluxeff = obj2->prof_disk_peak * 0.186682; /* e^-(b_n)*/
      if (FLAG(obj2.prof_disk_fluxmean))
        obj2->prof_disk_fluxmean = obj2->prof_disk_peak * 0.355007;/* b_n^(-2)*/
      }

/* Disk pattern */
    if (prefs.pattern_flag)
      {
      profit_residuals(profit, PROFIT_DYNPARAM, profit->paraminit,profit->resi);
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

  obj2->prof_flag = profit->flag;

  return;
  }


/****** profit_spread ********************************************************
PROTO	void profit_spread(profitstruct *profit, fieldstruct *field,
		fieldstruct *wfield, obj2struct *obj2)
PURPOSE	Perform star/galaxy separation by comparing the best-fitting local PSF
	and a compact exponential profile to the actual data using linear
	discriminant analysis.
INPUT	Pointer to the profile-fitting structure,
	pointer to the field,
	pointer to the field weight,
	pointer to the obj2.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/12/2013
 ***/
void	profit_spread(profitstruct *profit,  fieldstruct *field,
		fieldstruct *wfield, obj2struct *obj2)
  {
   profitstruct		*pprofit,*qprofit;
   subprofitstruct	*subprofit, *psubprofit,*qsubprofit;
   PIXTYPE		valp,valq, sig2;
   double		sump,sumq, sumpw2,sumqw2,sumpqw, sump0,sumq0;
   float		dchi2;
   int			p,s, nsub;

  nsub = profit->nsubprofit;

  pprofit = profit_init(obj2, MODEL_DIRAC, PROFIT_CONV);
  qprofit = profit_init(obj2, MODEL_EXPONENTIAL, PROFIT_CONV);

  profit_residuals(profit, PROFIT_DYNPARAM, profit->paraminit, profit->resi);
  profit_resetparams(pprofit);
  if (profit->paramlist[PARAM_X] && profit->paramlist[PARAM_Y])
    {
    pprofit->paraminit[pprofit->paramindex[PARAM_X]]
		= *profit->paramlist[PARAM_X];
    pprofit->paraminit[pprofit->paramindex[PARAM_Y]]
		= *profit->paramlist[PARAM_Y];
    }

  subprofit = profit->subprofit;
  for (s=0; s<nsub; s++)
    pprofit->paraminit[pprofit->paramindex[PARAM_DIRAC_FLUX+s]]
	= (subprofit++)->flux;

  pprofit->niter = profit_minimize(pprofit, PROFIT_MAXITER);

  profit_residuals(pprofit, PROFIT_DYNPARAM,pprofit->paraminit,pprofit->resi);

  qprofit->paraminit[qprofit->paramindex[PARAM_X]]
		= pprofit->paraminit[pprofit->paramindex[PARAM_X]];
  qprofit->paraminit[qprofit->paramindex[PARAM_Y]]
		= pprofit->paraminit[pprofit->paramindex[PARAM_Y]];
  for (s=0; s<nsub; s++)
    qprofit->paraminit[qprofit->paramindex[PARAM_DISK_FLUX+s]]
		= pprofit->paraminit[pprofit->paramindex[PARAM_DIRAC_FLUX+s]];
  qprofit->paraminit[qprofit->paramindex[PARAM_DISK_SCALE]]
	= profit->subprofit->psf->fwhm/16.0;	/* !CHECK */
  qprofit->paraminit[qprofit->paramindex[PARAM_DISK_ASPECT]] = 1.0;
  qprofit->paraminit[qprofit->paramindex[PARAM_DISK_POSANG]] = 0.0;

  profit_residuals(qprofit,PROFIT_DYNPARAM, qprofit->paraminit,qprofit->resi);

  if (FLAG(obj2.prof_class_star))
    {
    dchi2 = 0.5*(pprofit->chi2 - profit->chi2);
    obj2->prof_class_star = dchi2 < 50.0?
	(dchi2 > -50.0? 2.0/(1.0+expf(dchi2)) : 2.0) : 0.0;
    }

  psubprofit = pprofit->subprofit;
  qsubprofit = qprofit->subprofit;
  for (s=0; s<profit->nsubprofit; s++)
    {
    sump = sumq = sumpw2 = sumqw2 = sumpqw = sump0 = sumq0 = 0.0;
    for (p=0; p<psubprofit->nobjpix; p++)
      if (psubprofit->objweight[p]>0 && psubprofit->objpix[p]>-BIG)
        {
        valp = psubprofit->lmodpix[p];
        sump += (double)(valp*psubprofit->objpix[p]);
        valq = qsubprofit->lmodpix[p];
        sumq += (double)(valq*psubprofit->objpix[p]);
        sump0 += (double)(valp*valp);
        sumq0 += (double)(valp*valq);
        sig2 = 1.0f/(psubprofit->objweight[p]*psubprofit->objweight[p]);
        sumpw2 += valp*valp*sig2;
        sumqw2 += valq*valq*sig2;
        sumpqw += valp*valq*sig2;
        }

    if (FLAG(obj2.prof_concentration))
      {
      obj2->prof_concentration[s] = sump>0.0? (sumq/sump - sumq0/sump0) : 1.0;
      if (FLAG(obj2.prof_concentrationerr))
        obj2->prof_concentrationerr[s] = sump>0.0?
		sqrt(sumqw2*sump*sump+sumpw2*sumq*sumq-2.0*sumpqw*sump*sumq)
			/ (sump*sump) : 0.0;
      }
    }

  profit_end(pprofit);
  profit_end(qprofit);

  return;
  }

/****** subprofit_noisearea ***************************************************
PROTO	float subprofit_noisearea(profitstruct *profit)
PURPOSE	Return the equivalent noise area (see King 1983) of a model.
INPUT	Sub-profile-fitting structure,
OUTPUT	Equivalent noise area, in pixels.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
float	subprofit_noisearea(subprofitstruct *subprofit)
  {
   double	dval, flux,flux2;
   PIXTYPE	*pix;
   int		p;

  flux = flux2 = 0.0;
  pix = subprofit->lmodpix;
  for (p=subprofit->nobjpix; p--;)
    {
    dval = (double)*(pix++);
    flux += dval;
    flux2 += dval*dval;
    }

  return (float)(flux2>0.0? flux*flux / flux2 : 0.0);
  }


/****** subprofit_fluxcor ****************************************************
PROTO	void subprofit_fluxcor(subprofitstruct *subprofit, obj2struct *obj2)
PURPOSE	Integrate the flux within an ellipse and complete it with the wings of
	the fitted model.
INPUT	Sub-profile-fitting structure,
	pointer to the obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/03/2012
 ***/
void	subprofit_fluxcor(subprofitstruct *subprofit, obj2struct *obj2)
  {
    checkstruct		*check;
    double		mx,my, dx,dy, cx2,cy2,cxy, klim,klim2, tvobj,sigtvobj,
			tvm,tvmin,tvmout, r1,v1;
    PIXTYPE		*objpix,*objpixt,*objweight,*objweightt, *lmodpix,
			pix, weight,var;
    int			x,y, x2,y2, pos, w,h, area, corrflag;


  corrflag = (prefs.mask_type==MASK_CORRECT);
  w = subprofit->objnaxisn[0];
  h = subprofit->objnaxisn[1];
  mx = (float)(w/2);
  my = (float)(h/2);
/*
  if (FLAG(obj2.x_prof))
    {
    if (profit->paramlist[PARAM_X])
      mx += *profit->paramlist[PARAM_X];
    if (profit->paramlist[PARAM_Y])
      my += *profit->paramlist[PARAM_Y];
    }
*/
  if (obj2->auto_kronfactor>0.0)
    {
    cx2 = obj2->cxx;
    cy2 = obj2->cyy;
    cxy = obj2->cxy;
    klim2 = 0.64*obj2->auto_kronfactor*obj2->auto_kronfactor;
    }
  else
/*-- ...if not, use the circular aperture provided by the user */
    {
    cx2 = cy2 = 1.0;
    cxy = 0.0;
    klim2 = (prefs.autoaper[1]/2.0)*(prefs.autoaper[1]/2.0);
    }
/*
  cx2 = obj2->prof_convcxx;
  cy2 = obj2->prof_convcyy;
  cxy = obj2->prof_convcxy;

  lmodpix = profit->lmodpix;
  r1 = v1 = 0.0;
  for (y=0; y<h; y++)
    {
    dy = y - my;
    for (x=0; x<w; x++)
      {
      dx = x - mx;
      pix = *(lmodpix++);
      r1 += sqrt(cx2*dx*dx + cy2*dy*dy + cxy*dx*dy)*pix;
      v1 += pix;
      }
    }

  klim = r1/v1*2.0;
  klim2 = klim*klim;

if ((check = prefs.check[CHECK_APERTURES]))
sexellipse(check->pix, check->width, check->height,
obj2->x_prof-1.0, obj2->y_prof-1.0, klim*obj2->prof_conva,klim*obj2->prof_convb,
obj2->prof_convtheta, check->overlay, 0);
*/

  area = 0;
  tvmin = tvmout = tvobj = sigtvobj = 0.0;
  lmodpix = subprofit->lmodpix;
  objpixt = objpix = subprofit->objpix;
  objweightt = objweight = subprofit->objweight;
  for (y=0; y<h; y++)
    {
    for (x=0; x<w; x++, objpixt++,objweightt++)
      {
      dx = x - mx;
      dy = y - my;
      if ((cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
        {
        area++;
/*------ Here begin tests for pixel and/or weight overflows. Things are a */
/*------ bit intricated to have it running as fast as possible in the most */
/*------ common cases */
        if ((weight=*objweightt)<=0.0)
          {
          if (corrflag
		&& (x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (weight=objweight[pos = y2*w + x2])>0.0)
            {
            pix = objpix[pos];
            var = 1.0/(weight*weight);
            }
          else
            pix = var = 0.0;
          }
        else
          {
          pix = *objpixt;
          var = 1.0/(weight*weight);
          }
        tvobj += pix;
        sigtvobj += var;
        tvmin += *(lmodpix++);
//        *(lmodpix++) = pix;
        }
      else
        tvmout += *(lmodpix++);
      }
    }

//  tv -= area*bkg;

  tvm = tvmin + tvmout;
  if (tvm != 0.0)
    {
    obj2->fluxcor_prof = tvobj+obj2->flux_prof[0]*tvmout/tvm;
    obj2->fluxcorerr_prof = sqrt(sigtvobj
		+obj2->fluxerr_prof[0]*obj2->fluxerr_prof[0]*tvmout/tvm);
    }
  else
    {
    obj2->fluxcor_prof = tvobj;
    obj2->fluxcorerr_prof = sqrt(sigtvobj);
    }

/*
  if ((check = prefs.check[CHECK_OTHER]))
    check_add(check, profit->lmodpix, w, h, profit->ix,profit->iy, 1.0);
*/
  return;
  }


/****i* prof_gammainc *********************************************************
PROTO	double prof_gammainc(double x, double a)
PURPOSE	Returns the incomplete Gamma function (based on algorithm described in
	Numerical Recipes in C, chap. 6.1).
INPUT	A double,
	upper integration limit.
OUTPUT	Incomplete Gamma function.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/10/2010
*/
static double	prof_gammainc (double x, double a)

  {
   double	b,c,d,h, xn,xp, del,sum;
   int		i;

  if (a < 0.0 || x <= 0.0)
    return 0.0;

  if (a < (x+1.0))
    {
/*-- Use the series representation */
    xp = x;
    del = sum = 1.0/x;
    for (i=100;i--;)	/* Iterate to convergence */
      {
      sum += (del *= a/(++xp));
      if (fabs(del) < fabs(sum)*3e-7)
        return sum*exp(-a+x*log(a)) / prof_gamma(x);
      }
    }
  else
    {
/*-- Use the continued fraction representation and take its complement */
    b = a + 1.0 - x;
    c = 1e30;
    h = d = 1.0/b;
    for (i=1; i<=100; i++)	/* Iterate to convergence */
      {
      xn = -i*(i-x);
      b += 2.0;
      if (fabs(d=xn*d+b) < 1e-30)
        d = 1e-30;
      if (fabs(c=b+xn/c) < 1e-30)
        c = 1e-30;
      del= c * (d = 1.0/d);
      h *= del;
      if (fabs(del-1.0) < 3e-7)
        return 1.0 - exp(-a+x*log(a))*h / prof_gamma(x);
      }
    }
  error(EXIT_FAILURE, "*Error*: out of bounds in ",
		"prof_gammainc()");
  return 0.0;
  }


/****i* prof_gamma ************************************************************
PROTO	double prof_gamma(double xx)
PURPOSE	Returns the Gamma function (based on algorithm described in Numerical
	Recipes in C, chap 6.1).
INPUT	A double.
OUTPUT	Gamma function.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/09/2009
*/
static double	prof_gamma(double xx)

  {
   double		x,tmp,ser;
   static double	cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
   int			j;

  tmp=(x=xx-1.0)+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<6;j++)
    ser += cof[j]/(x+=1.0);

  return 2.50662827465*ser*exp(-tmp);
  }


/****** profit_minradius ******************************************************
PROTO	float profit_minradius(profitstruct *profit, float refffac)
PURPOSE	Returns the minimum disk radius that guarantees that each and
	every model component fits within some margin in that disk.
INPUT	Profit structure pointer,
	margin in units of (r/r_eff)^(1/n)).
OUTPUT	Radius (in pixels).
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/10/2010
*/
float	profit_minradius(profitstruct *profit, float refffac)

  {
   double	r,reff,rmax;
   int		p;

  rmax = reff = 0.0;
  for (p=0; p<profit->nprof; p++)
    {
    switch (profit->prof[p]->code)
      {
      case MODEL_BACK:
      case MODEL_DIRAC:
        reff = 0.0;
      break;
      case MODEL_SERSIC:
        reff = *profit->paramlist[PARAM_SPHEROID_REFF];
        break;
      case MODEL_DEVAUCOULEURS:
        reff = *profit->paramlist[PARAM_SPHEROID_REFF];
        break;
      case MODEL_EXPONENTIAL:
        reff = *profit->paramlist[PARAM_DISK_SCALE]*1.67835;
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown profile parameter in ",
		"profit_minradius()");
        break;
      }
    r = reff*(double)refffac;
    if (r>rmax)
      rmax = r;
    }

  return (float)rmax;
  }


/****** subprofit_psf *********************************************************
PROTO	void	subprofit_psf(subprofitstruct *subprofit, obj2struct *obj2)
PURPOSE	Build the local PSF at a given resolution.
INPUT	Pointer to sub-profile-fitting structure,
	pointer to obj2.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2012
 ***/
void	subprofit_psf(subprofitstruct *subprofit, obj2struct *obj2)
  {
   psfstruct	*psf;
   double	flux;
   float	posin[2], posout[2], dnaxisn[2],
		*pixout,
		xcout,ycout, xcin,ycin, invpixstep, norm;
   int		d,i;

  psf = subprofit->psf;
  psf_build(psf, obj2);

  xcout = (float)(subprofit->modnaxisn[0]/2) + 1.0;	/* FITS convention */
  ycout = (float)(subprofit->modnaxisn[1]/2) + 1.0;	/* FITS convention */
  xcin = (psf->masksize[0]/2) + 1.0;			/* FITS convention */
  ycin = (psf->masksize[1]/2) + 1.0;			/* FITS convention */
  invpixstep = subprofit->pixstep / psf->pixstep;

/* Initialize multi-dimensional counters */
  for (d=0; d<2; d++)
    {
    posout[d] = 1.0;					/* FITS convention */
    dnaxisn[d] = subprofit->modnaxisn[d]+0.5;
    }

/* Remap each pixel */
  pixout = subprofit->psfpix;
  flux = 0.0;
  for (i=subprofit->modnaxisn[0]*subprofit->modnaxisn[1]; i--;)
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
  flux *= subprofit->pixstep*subprofit->pixstep
		/ (subprofit->subsamp*subprofit->subsamp);
  if (fabs(flux) <= 0.0)
    error(EXIT_FAILURE, "*Error*: PSF model is empty or negative: ", psf->name);

  norm = 1.0/flux;
  pixout = subprofit->psfpix;
  for (i=subprofit->modnaxisn[0]*subprofit->modnaxisn[1]; i--;)
    *(pixout++) *= norm;  

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
VERSION	20/05/2011
 ***/
int	profit_minimize(profitstruct *profit, int niter)
  {
   double	lm_opts[5], info[LM_INFO_SZ],
		dcovar[PARAM_NPARAM*PARAM_NPARAM], dparam[PARAM_NPARAM];
   int		nfree;

  profit->iter = 0;
  memset(dcovar, 0, profit->nparam*profit->nparam*sizeof(double));

/* Perform fit */
  lm_opts[0] = 1.0e-3;		/* Initial mu */
  lm_opts[1] = 1.0e-8;		/* ||J^T e||_inf stopping factor */
  lm_opts[2] = 1.0e-8;		/* |Dp||_2 stopping factor */
  lm_opts[3] = 1.0e-8;		/* ||e||_2 stopping factor */
  lm_opts[4] = 1.0e-4;		/* Jacobian step */

  nfree = profit_boundtounbound(profit, profit->paraminit, dparam,
				PARAM_ALLPARAMS);

  niter = dlevmar_dif(profit_evaluate, dparam, NULL, nfree, profit->nresi,
			niter, lm_opts, info, NULL, dcovar, profit);

  profit_unboundtobound(profit, dparam, profit->paraminit, PARAM_ALLPARAMS);

/* Convert covariance matrix to bounded space */
  profit_covarunboundtobound(profit, dcovar, profit->covar);

  return niter;
  }


/****** profit_printout *******************************************************
PROTO	void profit_printout(int n_par, float* par, int m_dat, float* fvec,
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
VERSION	06/05/2012
 ***/
void	profit_printout(int n_par, float* par, int m_dat, float* fvec,
		void *data, int iflag, int iter, int nfev )
  {
   fieldstruct	*field;
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
    field = NULL;	/*! Should be replaced with a global variable to work!*/
    sprintf(filename, "check_%d_%04d.fits", the_gal, itero);
    check=check_init(filename, CHECK_PROFILES, 0, 1);
    check_reinit(field, check);
    check_add(check, profit->subprofit->lmodpix,
		profit->subprofit->objnaxisn[0],profit->subprofit->objnaxisn[1],
		profit->subprofit->ix,profit->subprofit->iy, 1.0);

    check_reend(field, check);
    check_end(check);
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
	pointer to a data structure (we use it for the profit structure here).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2011
 ***/
void	profit_evaluate(double *dpar, double *fvec, int m, int n, void *adata)
  {
   profitstruct		*profit;
   profstruct		**prof;
   double		*dpar0, *dresi;
   float		*modpixt, *profpixt, *resi,
			tparam, val;
   PIXTYPE		*lmodpixt,*lmodpix2t, *objpix,*weight,
			wval;
   int			c,f,i,p,q, fd,pd, jflag,sflag, nprof;

  profit = (profitstruct *)adata;

/* Detect "Jacobian-related" calls */
  jflag = pd = fd = 0;
  dpar0 = profit->dparam;
  if (profit->iter)
    {
    f = q = 0;
    for (p=0; p<profit->nparam; p++)
      {
      if (dpar[f] - dpar0[f] != 0.0)
        {
        pd = p;
        fd = f;
        q++;
        }
      if (profit->parfittype[p]!=PARFIT_FIXED)
        f++;
      }
    if (f>0 && q==1)
      jflag = 1;
    }
jflag = 0;	/* Temporarily deactivated (until problems are fixed) */
/*
  if (jflag && !(profit->nprof==1 && profit->prof[0]->code == MODEL_DIRAC))
    {
    prof = profit->prof;
    nprof = profit->nprof;

*-- "Jacobian call" *
    tparam = profit->param[pd];
    profit_unboundtobound(profit, &dpar[fd], &profit->param[pd], pd);
    sflag = 1;
    switch(profit->paramrevindex[pd])
      {
      case PARAM_BACK:
        lmodpixt = profit->lmodpix;
        lmodpix2t = profit->lmodpix2;
        val = (profit->param[pd] - tparam);
        for (i=profit->nobjpix;i--;)
          *(lmodpix2t++) = val;
        break;
      case PARAM_X:
      case PARAM_Y:
        profit_resample(profit, profit->cmodpix, profit->lmodpix2, 1.0);
        lmodpixt = profit->lmodpix;
        lmodpix2t = profit->lmodpix2;
        for (i=profit->nobjpix;i--;)
          *(lmodpix2t++) -= *(lmodpixt++);
        break;
      case PARAM_DIRAC_FLUX:
      case PARAM_SPHEROID_FLUX:
      case PARAM_DISK_FLUX:
      case PARAM_ARMS_FLUX:
      case PARAM_BAR_FLUX:
        if (nprof==1 && tparam != 0.0)
          {
          lmodpixt = profit->lmodpix;
          lmodpix2t = profit->lmodpix2;
          val = (profit->param[pd] - tparam) / tparam;
          for (i=profit->nobjpix;i--;)
            *(lmodpix2t++) = val**(lmodpixt++);
          }
        else
          {
          for (c=0; c<nprof; c++)
            if (prof[c]->flux[0] == &profit->param[pd])
              break;
          memcpy(profit->modpix, prof[c]->pix, profit->nmodpix*sizeof(float));
          profit_convolve(profit, profit->modpix);
          profit_resample(profit, profit->modpix, profit->lmodpix2,
		profit->param[pd] - tparam);
          }
        break;
      case PARAM_SPHEROID_REFF:
      case PARAM_SPHEROID_ASPECT:
      case PARAM_SPHEROID_POSANG:
      case PARAM_SPHEROID_SERSICN:
        sflag = 0;			* We are in the same switch *
        for (c=0; c<nprof; c++)
          if (prof[c]->code == MODEL_SERSIC
		|| prof[c]->code == MODEL_DEVAUCOULEURS)
            break; 
      case PARAM_DISK_SCALE:
      case PARAM_DISK_ASPECT:
      case PARAM_DISK_POSANG:
        if (sflag)
          for (c=0; c<nprof; c++)
            if (prof[c]->code == MODEL_EXPONENTIAL)
              break; 
        sflag = 0;
      case PARAM_ARMS_QUADFRAC:
      case PARAM_ARMS_SCALE:
      case PARAM_ARMS_START:
      case PARAM_ARMS_POSANG:
      case PARAM_ARMS_PITCH:
      case PARAM_ARMS_PITCHVAR:
      case PARAM_ARMS_WIDTH:
        if (sflag)
          for (c=0; c<nprof; c++)
            if (prof[c]->code == MODEL_ARMS)
              break; 
        sflag = 0;
      case PARAM_BAR_ASPECT:
      case PARAM_BAR_POSANG:
        if (sflag)
          for (c=0; c<nprof; c++)
            if (prof[c]->code == MODEL_ARMS)
              break; 
        modpixt = profit->modpix;
        profpixt = prof[c]->pix;
        val = -*prof[c]->flux[0];
        for (i=profit->nmodpix;i--;)
          *(modpixt++) = val**(profpixt++);
        memcpy(profit->modpix2, prof[c]->pix, profit->nmodpix*sizeof(float));
        prof_add(profit, prof[c], 0);
        memcpy(prof[c]->pix, profit->modpix2, profit->nmodpix*sizeof(float));
        profit_convolve(profit, profit->modpix);
        profit_resample(profit, profit->modpix, profit->lmodpix2, 1.0);
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: ",
			"unknown parameter index in profit_jacobian()");
        break;
      }
    objpix = profit->objpix;
    weight = profit->objweight;
    lmodpixt = profit->lmodpix;
    lmodpix2t = profit->lmodpix2;
    resi = profit->resi;
    dresi = fvec;
    if (PROFIT_DYNPARAM > 0.0)
      for (i=profit->nobjpix;i--; lmodpixt++, lmodpix2t++)
        {
        val = *(objpix++);
        if ((wval=*(weight++))>0.0)
          *(dresi++) = *(resi++) + *lmodpix2t
		* wval/(1.0+wval*fabs(*lmodpixt - val)/PROFIT_DYNPARAM);
        }
    else
      for (i=profit->nobjpix;i--; lmodpix2t++)
        if ((wval=*(weight++))>0.0)
          *(dresi++) = *(resi++) + *lmodpix2t * wval;
    }
  else
*/
    {
/*-- "Regular call" */
    for (p=0; p<profit->nparam; p++)
      dpar0[p] = dpar[p];
    profit_unboundtobound(profit, dpar, profit->param, PARAM_ALLPARAMS);

    profit_residuals(profit, PROFIT_DYNPARAM, profit->param, profit->resi);

    for (p=0; p<profit->nresi; p++)
      fvec[p] = profit->resi[p];
    }

//  profit_printout(m, par, n, fvec, adata, 0, -1, 0 );
  profit->iter++;

  return;
  }


/****** profit_residuals ******************************************************
PROTO	float *profit_residuals(profitstruct *profit, float dynparam,
			float *param, float *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	dynamic compression parameter (0=no compression),
	pointer to the model parameters (output),
	pointer to the computed residuals (output).
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/12/2013
 ***/
float	*profit_residuals(profitstruct *profit, float dynparam, float *param,
		float *resi)
  {
   subprofitstruct	*subprofit;
   float		*pixin, *pixout,
			pflux;
   int			i,p,s;

  subprofit = profit->subprofit;
  for (s=0; s<profit->nsubprofit; s++, subprofit++)
    memset(subprofit->modpix, 0, subprofit->nmodpix*sizeof(float));
  for (p=0; p<profit->nparam; p++)
    profit->param[p] = param[p];
/* Simple PSF shortcut */
  if (profit->conv_flag && profit->nprof == 1
	&& profit->prof[0]->code == MODEL_DIRAC)
    {
    subprofit = profit->subprofit;
    for (s=0; s<profit->nsubprofit; s++, subprofit++)
      {
      profit_resample(profit, subprofit, subprofit->psfpix,
		subprofit->lmodpix, *profit->prof[0]->flux[s]);
      subprofit->flux = *profit->prof[0]->flux[s];
      }
    }
  else
    {
    subprofit = profit->subprofit;
    for (s=0; s<profit->nsubprofit; s++, subprofit++)
      {
      subprofit->flux = 0.0;
      for (p=0; p<profit->nprof; p++)
        subprofit->flux += prof_add(subprofit, profit->prof[p], 0);
      memcpy(subprofit->cmodpix, subprofit->modpix,
		subprofit->nmodpix*sizeof(float));
      if (profit->conv_flag)
        subprofit_convolve(subprofit, subprofit->cmodpix);
      profit_resample(profit, subprofit,
	subprofit->cmodpix, subprofit->lmodpix, 1.0);
      }
    }

  if (resi)
    profit_compresi(profit, dynparam, resi);

  return resi;
  }


/****** profit_compresi ******************************************************
PROTO	float *profit_compresi(profitstruct *profit,float dynparam,float *resi)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Profile-fitting structure,
	dynamic-compression parameter (0=no compression),
	vector of residuals (output).
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
float	*profit_compresi(profitstruct *profit, float dynparam, float *resi)
  {
   subprofitstruct	*subprofit;
   double		error;
   float		*resit;
   PIXTYPE		*objpix, *objweight, *lmodpix,
			val,val2,wval, invsig;
   int			i,s;
  
/* Compute vector of residuals */
  resit = resi;
  error = 0.0;
  subprofit = profit->subprofit;
  if (dynparam > 0.0)
    {
    invsig = (PIXTYPE)(1.0/dynparam);
    for (s=profit->nsubprofit; s--; subprofit++)
      {
      objpix = subprofit->objpix;
      objweight = subprofit->objweight;
      lmodpix = subprofit->lmodpix;
      for (i=subprofit->objnaxisn[0]*subprofit->objnaxisn[1]; i--; lmodpix++)
        {
        val = *(objpix++);
        if ((wval=*(objweight++))>0.0)
          {
          val2 = (*lmodpix - val)*wval*invsig;
          val2 = val2>0.0? logf(1.0+val2) : -logf(1.0-val2);
          *(resit++) = val2*dynparam;
          error += val2*val2;
          }
        }
      }
    profit->chi2 = dynparam*dynparam*error;
    }
  else
    {
    for (s=profit->nsubprofit; s--; subprofit++)
      {
      objpix = subprofit->objpix;
      objweight = subprofit->objweight;
      lmodpix = subprofit->lmodpix;
      for (i=subprofit->objnaxisn[0]*subprofit->objnaxisn[1]; i--; lmodpix++)
        {
        val = *(objpix++);
        if ((wval=*(objweight++))>0.0)
          {
          val2 = (*lmodpix - val)*wval;
          *(resit++) = val2;
          error += val2*val2;
          }
        }
      }
    profit->chi2 = error;
    }

  return resi;
  }


/****** profit_resample ******************************************************
PROTO	int	prof_resample(profitstruct *profit, subprofitstruct *subprofit,
		float *inpix, PIXTYPE *outpix, float factor)
PURPOSE	Resample the current full resolution model to image resolution.
INPUT	Profile-fitting structure,
	sub-profile-fitting structure,
	pointer to input raster,
	pointer to output raster,
	multiplicating factor.
OUTPUT	RETURN_ERROR if the rasters don't overlap, RETURN_OK otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2013
 ***/
int	profit_resample(profitstruct *profit, subprofitstruct *subprofit,
		float *inpix, PIXTYPE *outpix, float factor)
  {
   PIXTYPE	*pixout,*pixout0;
   float	*pixin,*pixin0, *mask,*maskt, *pixinout, *dpixin,*dpixin0,
		*dpixout,*dpixout0, *dx,*dy,
		xcin,xcout,ycin,ycout, xsin,ysin, xin,yin, x,y, dxm,dym, val,
		invpixstep, norm;
   int		*start,*startt, *nmask,*nmaskt,
		i,j,k,n,t, 
		ixsout,iysout, ixout,iyout, dixout,diyout, nxout,nyout,
		iysina, nyin, hmw,hmh, ix,iy, ixin,iyin;

  invpixstep = subprofit->subsamp/subprofit->pixstep;
  factor /= subprofit->subsamp*subprofit->subsamp;

  xcin = (float)(subprofit->modnaxisn[0]/2);
  xcout = ((int)(subprofit->subsamp*subprofit->objnaxisn[0])/2 + 0.5)
		/ subprofit->subsamp - 0.5;
  if ((dx=profit->paramlist[PARAM_X]))
    xcout += *dx / subprofit->subsamp;

  xsin = xcin - xcout*invpixstep;			/* Input start x-coord*/

  if ((int)xsin >= subprofit->modnaxisn[0] || !finitef(xsin))
    return RETURN_ERROR;
  ixsout = 0;				/* Int. part of output start x-coord */
  if (xsin<0.0)
    {
    dixout = (int)(1.0-xsin/invpixstep);
/*-- Simply leave here if the images do not overlap in x */
    if (dixout >= subprofit->objnaxisn[0])
      return RETURN_ERROR;
    ixsout += dixout;
    xsin += dixout*invpixstep;
    }
  nxout = (int)((subprofit->modnaxisn[0]-xsin)/invpixstep);/* nb of interpolated
							input pixels along x */
  if (nxout>(ixout=subprofit->objnaxisn[0]-ixsout))
    nxout = ixout;
  if (!nxout)
    return RETURN_ERROR;

  ycin = (float)(subprofit->modnaxisn[1]/2);
  ycout = ((int)(subprofit->subsamp*subprofit->objnaxisn[1])/2 + 0.5)
		/ subprofit->subsamp - 0.5;
  if ((dy=profit->paramlist[PARAM_Y]))
    ycout += *dy / subprofit->subsamp;

  ysin = ycin - ycout*invpixstep;		/* Input start y-coord*/

  if ((int)ysin >= subprofit->modnaxisn[1] || !finitef(ysin))
    return RETURN_ERROR;
  iysout = 0;				/* Int. part of output start y-coord */
  if (ysin<0.0)
    {
    diyout = (int)(1.0-ysin/invpixstep);
/*-- Simply leave here if the images do not overlap in y */
    if (diyout >= subprofit->objnaxisn[1])
      return RETURN_ERROR;
    iysout += diyout;
    ysin += diyout*invpixstep;
    }
  nyout = (int)((subprofit->modnaxisn[1]-ysin)/invpixstep);/* nb of interpolated
							input pixels along y */
  if (nyout>(iyout=subprofit->objnaxisn[1]-iysout))
    nyout = iyout;
  if (!nyout)
    return RETURN_ERROR;

/* Set the yrange for the x-resampling with some margin for interpolation */
  iysina = (int)ysin;	/* Int. part of Input start y-coord with margin */
  hmh = INTERPW/2 - 1;	/* Interpolant start */
  if (iysina<0 || ((iysina -= hmh)< 0))
    iysina = 0;
  nyin = (int)(ysin+nyout*invpixstep)+INTERPW-hmh;/* Interpolated Input y size*/
  if (nyin>subprofit->modnaxisn[1])		/* with margin */
    nyin = subprofit->modnaxisn[1];
/* Express everything relative to the effective Input start (with margin) */
  nyin -= iysina;
  ysin -= (float)iysina;

/* Allocate interpolant stuff for the x direction */
  QMALLOC(mask, float, nxout*INTERPW);	/* Interpolation masks */
  QMALLOC(nmask, int, nxout);		/* Interpolation mask sizes */
  QMALLOC(start, int, nxout);		/* Int. part of Input conv starts */
/* Compute the local interpolant and data starting points in x */
  hmw = INTERPW/2 - 1;
  xin = xsin;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=nxout; j--; xin+=invpixstep)
    {
    ix = (ixin=(int)xin) - hmw;
    dxm = ixin - xin - hmw;	/* starting point in the interpolation func */
    if (ix < 0)
      {
      n = INTERPW+ix;
      dxm -= (float)ix;
      ix = 0;
      }
    else
      n = INTERPW;
    if (n>(t=subprofit->modnaxisn[0]-ix))
      n=t;
    *(startt++) = ix;
    *(nmaskt++) = n;
    norm = 0.0;
    for (x=dxm, i=n; i--; x+=1.0)
      norm += (*(maskt++) = INTERPF(x));
    norm = norm>0.0? 1.0/norm : 1.0;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

  QCALLOC(pixinout, float, nxout*nyin);	/* Intermediary frame-buffer */

/* Make the interpolation in x (this includes transposition) */
  pixin0 = inpix + iysina*subprofit->modnaxisn[0];
  dpixout0 = pixinout;
  for (k=nyin; k--; pixin0+=subprofit->modnaxisn[0], dpixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    dpixout = dpixout0;
    for (j=nxout; j--; dpixout+=nyin)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(pixin++);
      *dpixout = val;
      }
    }

/* Reallocate interpolant stuff for the y direction */
  QREALLOC(mask, float, nyout*INTERPW);	/* Interpolation masks */
  QREALLOC(nmask, int, nyout);			/* Interpolation mask sizes */
  QREALLOC(start, int, nyout);		/* Int. part of Input conv starts */

/* Compute the local interpolant and data starting points in y */
  hmh = INTERPW/2 - 1;
  yin = ysin;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=nyout; j--; yin+=invpixstep)
    {
    iy = (iyin=(int)yin) - hmh;
    dym = iyin - yin - hmh;	/* starting point in the interpolation func */
    if (iy < 0)
      {
      n = INTERPW+iy;
      dym -= (float)iy;
      iy = 0;
      }
    else
      n = INTERPW;
    if (n>(t=nyin-iy))
      n=t;
    *(startt++) = iy;
    *(nmaskt++) = n;
    norm = 0.0;
    for (y=dym, i=n; i--; y+=1.0)
      norm += (*(maskt++) = INTERPF(y));
    norm = norm>0.0? 1.0/norm : 1.0;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

/* Initialize destination buffer to zero */
  memset(outpix, 0, (size_t)subprofit->nobjpix*sizeof(PIXTYPE));

/* Make the interpolation in y and transpose once again */
  dpixin0 = pixinout;
  pixout0 = outpix+ixsout+iysout*subprofit->objnaxisn[0];
  for (k=nxout; k--; dpixin0+=nyin, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=nyout; j--; pixout+=subprofit->objnaxisn[0])
      {
      dpixin = dpixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(dpixin++);
       *pixout = (PIXTYPE)(factor*val);
      }
    }

/* Free memory */
  free(pixinout);
  free(mask);
  free(nmask);
  free(start);

  return RETURN_OK;
  }


/****** subprofit_convolve ****************************************************
PROTO	void subprofit_convolve(subprofitstruct *subprofit, float *modpix)
PURPOSE	Convolve a model image with the local PSF.
INPUT	Pointer to the subprofit structure,
	Pointer to the image raster.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2012
 ***/
void	subprofit_convolve(subprofitstruct *subprofit, float *modpix)
  {
  if (!subprofit->psfdft)
    subprofit_makedft(subprofit);

  fft_conv(modpix, subprofit->psfdft, subprofit->modnaxisn,
		&subprofit->fftscratch);

  return;
  }


/****** subprofit_makedft ****************************************************
PROTO	void subprofit_makedft(subprofitstruct *subprofit)
PURPOSE	Create the Fourier transform of the descrambled PSF component.
INPUT	Pointer to the subprofit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
void	subprofit_makedft(subprofitstruct *subprofit)
  {
   psfstruct	*psf;
   float      *mask,*maskt, *ppix;
   float       dx,dy, r,r2,rmin,rmin2,rmax,rmax2,rsig,invrsig2;
   int          width,height,npix,offset, psfwidth,psfheight,psfnpix,
                cpwidth, cpheight,hcpwidth,hcpheight, i,j,x,y;

  if (!(psf=subprofit->psf))
    return;

  psfwidth = subprofit->modnaxisn[0];
  psfheight = subprofit->modnaxisn[1];
  psfnpix = psfwidth*psfheight;
  width = subprofit->modnaxisn[0];
  height = subprofit->modnaxisn[1];
  npix = width*height;
  QCALLOC(mask, float, npix);
  cpwidth = (width>psfwidth)?psfwidth:width;
  hcpwidth = cpwidth>>1;
  cpwidth = hcpwidth<<1;
  offset = width - cpwidth;
  cpheight = (height>psfheight)?psfheight:height;
  hcpheight = cpheight>>1;
  cpheight = hcpheight<<1;

/* Frame and descramble the PSF data */
  ppix = subprofit->psfpix + (psfheight/2)*psfwidth + psfwidth/2;
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

  ppix = subprofit->psfpix + ((psfheight/2)-hcpheight)*psfwidth + psfwidth/2;
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
  rsig = psf->fwhm/subprofit->pixstep;
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
        *maskt *= (r2>rmax2)?0.0:expf((2*rmin*sqrtf(r2)-r2-rmin2)*invrsig2);
    dx = -hcpwidth;
    maskt += offset;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:expf((2*rmin*sqrtf(r2)-r2-rmin2)*invrsig2);
    }
  dy = -hcpheight;
  maskt += width*(height-cpheight);
  for (y=hcpheight; y--; dy+=1.0)
    {
    dx = 0.0;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:expf((2*rmin*sqrtf(r2)-r2-rmin2)*invrsig2);
    dx = -hcpwidth;
    maskt += offset;
    for (x=hcpwidth; x--; dx+=1.0, maskt++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *maskt *= (r2>rmax2)?0.0:expf((2*rmin*sqrtf(r2)-r2-rmin2)*invrsig2);
    }

/* Finally move to Fourier space */
  subprofit->psfdft = fft_rtf(mask, subprofit->modnaxisn);

  free(mask);

  return;
  }


/****** subprofit_copyobjpix **************************************************
PROTO	int subprofit_copyobjpix(subprofitstruct *subprofit,
				subimagestruct *subimage)
PURPOSE	Copy a piece of the input object image to a subprofit structure.
INPUT	Pointer to the subprofit structure,
	subimagestruct *subimage.
OUTPUT	The number of valid pixels copied.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
int	subprofit_copyobjpix(subprofitstruct *subprofit,
				subimagestruct *subimage)
  {
   fieldstruct		*field, *wfield;
   float		dx, dy2, dr2, rad2;
   PIXTYPE		*pixin,*spixin, *wpixin,*swpixin, *pixout,*wpixout,
			backnoise2, invgain, satlevel, wthresh, pix,spix,
			wpix,swpix;
   int			i,x,y, xmin,xmax,ymin,ymax, w,h,dw, npix, off,
			gainflag,badflag, sflag, sx,sy,sn, ix,iy, win,hin;

  field = subimage->field;
  wfield = subimage->wfield;

/* First put the image background to -BIG */
  pixout = subprofit->objpix;
  wpixout = subprofit->objweight;
  for (i=subprofit->objnaxisn[0]*subprofit->objnaxisn[1]; i--;)
    {
    *(pixout++) = -BIG;
    *(wpixout++) = 0.0;
    }
 
  ix = subprofit->ix - subimage->xmin[0];
  iy = subprofit->iy - subimage->xmin[1];
  win = subimage->size[0];
  hin = subimage->size[1];
 
  backnoise2 = field->backsig*field->backsig;
  sn = (int)subprofit->subsamp;
  sflag = (sn>1);
  w = subprofit->objnaxisn[0]*sn;
  h = subprofit->objnaxisn[1]*sn;
  if (sflag)
    backnoise2 *= (PIXTYPE)sn;
  invgain = (field->gain > 0.0) ? 1.0/field->gain : 0.0;
  satlevel = field->satur_level - subimage->bkg;
  rad2 = h/2.0;
  if (rad2 > w/2.0)
    rad2 = w/2.0;
  rad2 *= rad2;
 
/* Set the image boundaries */
  pixout = subprofit->objpix;
  wpixout = subprofit->objweight;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<0)
    {
    off = (-ymin-1)/sn + 1;
    pixout += off*subprofit->objnaxisn[0];
    wpixout += off*subprofit->objnaxisn[0];
    ymin += off*sn;
    }
  if (ymax>hin)
    ymax -= ((ymax-hin-1)/sn + 1)*sn;
 
  xmin = ix-w/2;
  xmax = xmin + w;
  dw = 0;
  if (xmax>win)
    {
    off = (xmax-win-1)/sn + 1;
    dw += off;
    xmax -= off*sn;
    }
  if (xmin<0)
    {
    off = (-xmin-1)/sn + 1;
    pixout += off;
    wpixout += off;
    dw += off;
    xmin += off*sn;
    }
 
/* Copy the right pixels to the destination */
  npix = 0;
  if (wfield)
    {
    wthresh = wfield->weight_thresh;
    gainflag = field->weightgain_flag;
    if (sflag)
      {
/*---- Sub-sampling case */
      for (y=ymin; y<ymax; y+=sn, pixout+=dw,wpixout+=dw)
        {
        for (x=xmin; x<xmax; x+=sn)
          {
          pix = wpix = 0.0;
          badflag = 0;
          for (sy=0; sy<sn; sy++)
            {
            dy2 = (y+sy-iy);
            dy2 *= dy2;
            dx = (x-ix);
            off = x + (y+sy)*win;
            spixin = subimage->image + off;
            swpixin = subimage->weight + off;
            for (sx=sn; sx--;)
              {
              dr2 = dy2 + dx*dx;
              dx++;
              spix = *(spixin++);
              swpix = *(swpixin++);
              if (dr2<rad2 && spix>-BIG && spix<satlevel && swpix<wthresh)
                {
                pix += spix;
                wpix += swpix;
                }
              else
                badflag=1;
              }
            }
          *(pixout++) = pix;
          if (!badflag) /* A single bad pixel ruins is all (saturation, etc.)*/
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
        pixin = subimage->image + xmin + y*win;
        wpixin = subimage->weight + xmin + y*win;
        for (x=xmin; x<xmax; x++)
          {
          dx = x-ix;
          dr2 = dy2 + dx*dx;
          pix = *(pixin++);
          wpix = *(wpixin++);
          if (dr2<rad2 && pix>-BIG && pix<satlevel && wpix<wthresh)
            {
            *(pixout++) = pix;
            *(wpixout++) = 1.0 / sqrt(wpix+(pix>0.0?
		(gainflag? pix*wpix/backnoise2:pix)*invgain : 0.0));
            npix++;
            }
          else
            *(pixout++) = *(wpixout++) = 0.0;
          }
        }
    }
  else
    {
    if (sflag)
      {
/*---- Sub-sampling case */
      for (y=ymin; y<ymax; y+=sn, pixout+=dw, wpixout+=dw)
        {
        for (x=xmin; x<xmax; x+=sn)
          {
          pix = 0.0;
          badflag = 0;
          for (sy=0; sy<sn; sy++)
            {
            dy2 = y+sy-iy;
            dy2 *= dy2;
            dx = x-ix;
            spixin = subimage->image + x + (y+sy)*win;
            for (sx=sn; sx--;)
              {
              dr2 = dy2 + dx*dx;
              dx++;
              spix = *(spixin++);
              if (dr2<rad2 && spix>-BIG && spix<satlevel)
                pix += spix;
              else
                badflag=1;
              }
            }
          *(pixout++) = pix;
          if (!badflag) /* A single bad pixel ruins is all (saturation, etc.)*/
            {
            *(wpixout++) = 1.0 / sqrt(backnoise2 + (pix>0.0?pix*invgain:0.0));
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
        pixin = subimage->image + xmin + y*win;
        for (x=xmin; x<xmax; x++)
          {
          dx = x-ix;
          dr2 = dy2 + dx*dx;
          pix = *(pixin++);
          if (dr2<rad2 && pix>-BIG && pix<satlevel)
            {
            *(pixout++) = pix;
            *(wpixout++) = 1.0 / sqrt(backnoise2 + (pix>0.0?pix*invgain : 0.0));
            npix++;
            }
          else
            *(pixout++) = *(wpixout++) = 0.0;
          }
        }
    }
 
  return npix;
  }
 

/****** subprofit_submodpix ***************************************************
PROTO	void	subprofit_submodpix(subprofitstruct *subprofitmod,
		PIXTYPE *pixout, int ix, int iy, int width, int height,
		float oversamp, float fac)
PURPOSE	Subtract a rasterized model from an image.
INPUT	Pointer to the sub-profile structure containing the model raster,
	image pointer,
	destination image center x coordinate (integer),
	destination image center y coordinate (integer),
	destination image width,
	destination image height,
	oversampling factor,
	multiplicative factor to apply to the model pixels.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2014
 ***/
void	subprofit_submodpix(subprofitstruct *subprofitmod,
		PIXTYPE *pixout, int ix, int iy, int width, int height,
		float oversamp, float fac)
  {
   float	dx, dy2, dr2, rad2;
   PIXTYPE	*pixin,*spixin,
		pix,spix;
   int		i,x,y, xmin,xmax,ymin,ymax, w,h,dw, npix, off, gainflag,
		badflag, sflag, sx,sy,sn, win,hin;

/* Don't go further if out of frame!! */
  win = subprofitmod->objnaxisn[0];
  hin = subprofitmod->objnaxisn[1];
  ix -= subprofitmod->ix - win/2;
  iy -= subprofitmod->iy - hin/2;

  sn = (int)oversamp;
  sflag = (sn>1);
  w = width*sn;
  h = height*sn;
  rad2 = h/2.0;
  if (rad2 > w/2.0)
    rad2 = w/2.0;
  rad2 *= rad2;

/* Set the image boundaries */
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<0)
    {
    off = (-ymin-1)/sn + 1;
    pixout += off*width;
    ymin += off*sn;
    }
  if (ymax>hin)
    ymax -= ((ymax-hin-1)/sn + 1)*sn;

  xmin = ix-w/2;
  xmax = xmin + w;
  dw = 0;
  if (xmax>win)
    {
    off = (xmax-win-1)/sn + 1;
    dw += off;
    xmax -= off*sn;
    }
  if (xmin<0)
    {
    off = (-xmin-1)/sn + 1;
    pixout += off;
    dw += off;
    xmin += off*sn;
    }

/* Copy the right pixels to the destination */
  npix = 0;
  if (sflag)
    {
/*-- Sub-sampling case */
    for (y=ymin; y<ymax; y+=sn, pixout+=dw)
      {
      for (x=xmin; x<xmax; x+=sn)
        {
        pix = 0.0;
        badflag = 0;
        for (sy=0; sy<sn; sy++)
          {
          dy2 = y+sy-iy;
          dy2 *= dy2;
          dx = x-ix;
          spixin = subprofitmod->lmodpix + x + (y+sy)*win;
          for (sx=sn; sx--;)
            {
            dr2 = dy2 + dx*dx;
            dx++;
            spix = *(spixin++);
            pix += spix;
            }
          }
        *(pixout++) -= fac*pix;
        }
      }
    }
  else
    for (y=ymin; y<ymax; y++, pixout+=dw)
      {
      dy2 = y-iy;
      dy2 *= dy2;
      pixin = subprofitmod->lmodpix + xmin + y*win;
      for (x=xmin; x<xmax; x++)
        {
        dx = x-ix;
        dr2 = dy2 + dx*dx;
        *(pixout++) -= fac**(pixin++);
        }
      }

  return;
  }


/****** subprofit_spiralindex *************************************************
PROTO	float subprofit_spiralindex(subprofitstruct *subprofit)
PURPOSE	Compute the spiral index of a galaxy image (positive for arms
	extending counter-clockwise and negative for arms extending CW, 0 for
	no spiral pattern).
INPUT	Profile-fitting structure.
OUTPUT	Vector of residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/08/2012
 ***/
float subprofit_spiralindex(subprofitstruct *subprofit, obj2struct *obj2)
  {
   float	*dx,*dy, *fdx,*fdy, *gdx,*gdy, *gdxt,*gdyt, *pix,
		fwhm, invtwosigma2, hw,hh, ohw,ohh, x,y,xstart, tx,ty,txstart,
		gx,gy, r2, spirindex, invsig, val, sep;
   PIXTYPE	*fpix;
   int		i,j, npix;

  npix = subprofit->objnaxisn[0]*subprofit->objnaxisn[1];

/* Compute simple derivative vectors at a fraction of the object scale */
  fwhm = subprofit->guessradius * 2.0 / 4.0;
  if (fwhm < 2.0)
    fwhm = 2.0;
  sep = 2.0;

  invtwosigma2 = -(2.35*2.35/(2.0*fwhm*fwhm));
  hw = (float)(subprofit->objnaxisn[0]/2);
  ohw = subprofit->objnaxisn[0] - hw;
  hh = (float)(subprofit->objnaxisn[1]/2);
  ohh = subprofit->objnaxisn[1] - hh;
  txstart = -hw;
  ty = -hh;
  QFFTWF_MALLOC(dx, float, npix);
  pix = dx;
  for (j=subprofit->objnaxisn[1]; j--; ty+=1.0)
    {
    tx = txstart;
    y = ty < -0.5? ty + hh : ty - ohh;
    for (i=subprofit->objnaxisn[0]; i--; tx+=1.0)
      {
      x = tx < -0.5? tx + hw : tx - ohw;
      *(pix++) = exp(invtwosigma2*((x+sep)*(x+sep)+y*y))
		- exp(invtwosigma2*((x-sep)*(x-sep)+y*y));
      }
    }
  QFFTWF_MALLOC(dy, float, npix);
  pix = dy;
  ty = -hh;
  for (j=subprofit->objnaxisn[1]; j--; ty+=1.0)
    {
    tx = txstart;
    y = ty < -0.5? ty + hh : ty - ohh;
    for (i=subprofit->objnaxisn[0]; i--; tx+=1.0)
      {
      x = tx < -0.5? tx + hw : tx - ohw;
      *(pix++) = exp(invtwosigma2*(x*x+(y+sep)*(y+sep)))
		- exp(invtwosigma2*(x*x+(y-sep)*(y-sep)));
      }
    }

  QFFTWF_MALLOC(gdx, float, npix);
  gdxt = gdx;
  fpix = subprofit->objpix;
  invsig = npix/subprofit->sigma;
  for (i=npix; i--; fpix++)
    {
    val = *fpix > -1e29? *fpix*invsig : 0.0;
    *(gdxt++) = (val>0.0? log(1.0+val) : -log(1.0-val));
    }

  QFFTWF_MALLOC(gdy, float, npix);
  memcpy(gdy, gdx, npix*sizeof(float));
  fdx = fft_rtf(dx, subprofit->objnaxisn);
  fft_conv(gdx, fdx, subprofit->objnaxisn, &subprofit->fftscratch);
  fdy = fft_rtf(dy, subprofit->objnaxisn);
  fft_conv(gdy, fdy, subprofit->objnaxisn, &subprofit->fftscratch);

/* Compute estimator */
  invtwosigma2 = -1.18*1.18/(2.0*subprofit->guessradius*subprofit->guessradius);
  xstart = -hw - obj2->mx + (int)(obj2->mx+0.49999);
  y = -hh -  obj2->my + (int)(obj2->my+0.49999);;
  spirindex = 0.0;
  gdxt = gdx;
  gdyt = gdy;
  for (j=subprofit->objnaxisn[1]; j--; y+=1.0)
    {
    x = xstart;
    for (i=subprofit->objnaxisn[0]; i--; x+=1.0)
      {
      gx = *(gdxt++);
      gy = *(gdyt++);
      if ((r2=x*x+y*y)>0.0)
        spirindex += (x*y*(gx*gx-gy*gy)+gx*gy*(y*y-x*x))/r2
			* exp(invtwosigma2*r2);
      }
    }

  QFFTWF_FREE(dx);
  QFFTWF_FREE(dy);
  QFFTWF_FREE(fdx);
  QFFTWF_FREE(fdy);
  QFFTWF_FREE(gdx);
  QFFTWF_FREE(gdy);

  return spirindex;
  }


/****** profit_moments *******************************************************
PROTO	void profit_moments(profitstruct *profit, obj2struct *obj2)
PURPOSE	Compute the 2nd order moments from the unconvolved object model.
INPUT	Profile-fitting structure,
	Pointer to obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	22/04/2011
 ***/
void	 profit_moments(profitstruct *profit, obj2struct *obj2)
  {
   profstruct	*prof;
   double	dpdmx2[6], cov[4],
		*jac,*jact, *pjac,*pjact, *dcovar,*dcovart,
		*dmx2,*dmy2,*dmxy,
		m0,invm0, mx2,my2,mxy, den,invden,
		temp, temp2,invtemp2,invstemp2,
		pmx2,theta, flux, dval;
   float	 *covart;
   int		findex[MODEL_NMAX],
		i,j,p, nparam;

/*  hw = (float)(profit->modnaxisn[0]/2);*/
/*  hh = (float)(profit->modnaxisn[1]/2);*/
/*  r2max = hw<hh? hw*hw : hh*hh;*/
/*  xstart = -hw;*/
/*  y = -hh;*/
/*  pix = profit->modpix;*/
/*  mx2 = my2 = mxy = mx = my = sum = 0.0;*/
/*  for (iy=profit->modnaxisn[1]; iy--; y+=1.0)*/
/*    {*/
/*    x = xstart;*/
/*    for (ix=profit->modnaxisn[0]; ix--; x+=1.0)*/
/*      if (y*y+x*x <= r2max)*/
/*        {*/
/*        val = *(pix++);*/
/*        sum += val;*/
/*        mx  += val*x;*/
/*        my  += val*y;*/
/*        mx2 += val*x*x;*/
/*        mxy += val*x*y;*/
/*        my2 += val*y*y;*/
/*        }*/
/*      else*/
/*        pix++;*/
/*    }*/

/*  if (sum <= 1.0/BIG)*/
/*    sum = 1.0;*/
/*  mx /= sum;*/
/*  my /= sum;*/
/*  obj2->prof_mx2 = mx2 = mx2/sum - mx*mx;*/
/*  obj2->prof_my2 = my2 = my2/sum - my*my;*/
/*  obj2->prof_mxy = mxy = mxy/sum - mx*my;*/

  nparam = profit->nparam;
  if (FLAG(obj2.prof_e1err) || FLAG(obj2.prof_pol1err))
    {
/*-- Set up Jacobian matrices */
    QCALLOC(jac, double, nparam*3);
    QMALLOC(pjac, double, (nparam<2? 6 : nparam*3));
    QMALLOC(dcovar, double, nparam*nparam);
    dcovart = dcovar;
    covart = profit->covar;
    for (i=nparam*nparam; i--;)
      *(dcovart++) = (double)(*(covart++));
    dmx2 = jac;
    dmy2 = jac+nparam;
    dmxy = jac+2*nparam;
    }
  else
    jac = pjac = dcovar = dmx2 = dmy2 = dmxy = NULL;

  m0 = mx2 = my2 = mxy = 0.0;
  for (p=0; p<profit->nprof; p++)
    {
    prof = profit->prof[p];
    findex[p] = prof_moments(profit, prof, pjac);
    flux = *prof->flux[0];
    m0 += flux;
    mx2 += prof->mx2*flux;
    my2 += prof->my2*flux;
    mxy += prof->mxy*flux;
    if (jac)
      {
      jact = jac;
      pjact = pjac;
      for (j=nparam*3; j--;)
        *(jact++) += flux * *(pjact++);
      }
    }
  invm0 = 1.0 / m0;
  obj2->prof_mx2 = (mx2 *= invm0);
  obj2->prof_my2 = (my2 *= invm0);
  obj2->prof_mxy = (mxy *= invm0);
/* Complete the flux derivative of moments */
  if (jac)
    {
    for (p=0; p<profit->nprof; p++)
      {
      prof = profit->prof[p];
      dmx2[findex[p]] = prof->mx2 - mx2;
      dmy2[findex[p]] = prof->my2 - my2;
      dmxy[findex[p]] = prof->mxy - mxy;
      }
    jact = jac;
    for (j=nparam*3; j--;)
      *(jact++) *= invm0;
    }

/* Handle fully correlated profiles (which cause a singularity...) */
  if ((temp2=mx2*my2-mxy*mxy)<0.00694)
    {
    mx2 += 0.0833333;
    my2 += 0.0833333;
    temp2 = mx2*my2-mxy*mxy;
    }

/* Use the Jacobians to compute the moment covariance matrix */
  if (jac)
    propagate_covar(dcovar, jac, obj2->prof_mx2cov, nparam, 3,
						pjac);	/* We re-use pjac */

  if (FLAG(obj2.prof_pol1))
    {
/*--- "Polarisation", i.e. module = (a^2-b^2)/(a^2+b^2) */
    if (mx2+my2 > 1.0/BIG)
      {
      obj2->prof_pol1 = (mx2 - my2) / (mx2+my2);
      obj2->prof_pol2 = 2.0*mxy / (mx2 + my2);
      if (FLAG(obj2.prof_pol1err))
        {
/*------ Compute the Jacobian of polarisation */
        invden = 1.0/(mx2+my2);
        dpdmx2[0] =  2.0*my2*invden*invden;
        dpdmx2[1] = -2.0*mx2*invden*invden;
        dpdmx2[2] =  0.0;
        dpdmx2[3] = -2.0*mxy*invden*invden;
        dpdmx2[4] = -2.0*mxy*invden*invden;
        dpdmx2[5] =  2.0*invden;

/*------ Use the Jacobian to compute the polarisation covariance matrix */
        propagate_covar(obj2->prof_mx2cov, dpdmx2, cov, 3, 2,
						pjac);	/* We re-use pjac */
        obj2->prof_pol1err = (float)sqrt(cov[0]<0.0? 0.0: cov[0]);
        obj2->prof_pol2err = (float)sqrt(cov[3]<0.0? 0.0: cov[3]);
        obj2->prof_pol12corr = (dval=cov[0]*cov[3]) > 0.0?
					(float)(cov[1]/sqrt(dval)) : 0.0;
        }
      }
    else
      obj2->prof_pol1 = obj2->prof_pol2
	= obj2->prof_pol1err = obj2->prof_pol2err = obj2->prof_pol12corr = 0.0;
    }

  if (FLAG(obj2.prof_e1))
    {
/*--- "Ellipticity", i.e. module = (a-b)/(a+b) */
    if (mx2+my2 > 1.0/BIG)
      {
      den = (temp2>=0.0) ? mx2+my2+2.0*sqrt(temp2) : mx2+my2;
      invden = 1.0/den;
      obj2->prof_e1 = (float)(invden * (mx2 - my2));
      obj2->prof_e2 = (float)(2.0 * invden * mxy);
      if (FLAG(obj2.prof_e1err))
        {
/*------ Compute the Jacobian of ellipticity */
        invstemp2 = (temp2>=0.0) ? 1.0/sqrt(temp2) : 0.0;
        dpdmx2[0] = ( den - (1.0+my2*invstemp2)*(mx2-my2))*invden*invden;
        dpdmx2[1] = (-den - (1.0+mx2*invstemp2)*(mx2-my2))*invden*invden;
        dpdmx2[2] = 2.0*mxy*invstemp2*(mx2-my2)*invden*invden;
        dpdmx2[3] = -2.0*mxy*(1.0+my2*invstemp2)*invden*invden;
        dpdmx2[4] = -2.0*mxy*(1.0+mx2*invstemp2)*invden*invden;
        dpdmx2[5] =  (2.0*den+4.0*mxy*mxy*invstemp2)*invden*invden;

/*------ Use the Jacobian to compute the ellipticity covariance matrix */
        propagate_covar(obj2->prof_mx2cov, dpdmx2, cov, 3, 2,
					pjac);	/* We re-use pjac */
        obj2->prof_e1err = (float)sqrt(cov[0]<0.0? 0.0: cov[0]);
        obj2->prof_e2err = (float)sqrt(cov[3]<0.0? 0.0: cov[3]);
        obj2->prof_e12corr = (dval=cov[0]*cov[3]) > 0.0?
					(float)(cov[1]/sqrt(dval)) : 0.0;
        }
      }
    else
      obj2->prof_e1 = obj2->prof_e2
	= obj2->prof_e1err = obj2->prof_e2err = obj2->prof_e12corr = 0.0;
    }

  if (FLAG(obj2.prof_cxx))
    {
    invtemp2 = (temp2>=0.0) ? 1.0/temp2 : 0.0;
    obj2->prof_cxx = (float)(my2*invtemp2);
    obj2->prof_cyy = (float)(mx2*invtemp2);
    obj2->prof_cxy = (float)(-2*mxy*invtemp2);
    }

  if (FLAG(obj2.prof_a))
    {
    if ((fabs(temp=mx2-my2)) > 0.0)
      theta = atan2(2.0 * mxy,temp) / 2.0;
    else
      theta = PI/4.0;

    temp = sqrt(0.25*temp*temp+mxy*mxy);
    pmx2 = 0.5*(mx2+my2);
    obj2->prof_a = (float)sqrt(pmx2 + temp);
    obj2->prof_b = (float)sqrt(pmx2 - temp);
    obj2->prof_theta = theta*180.0/PI;
    }

/* Free memory used by Jacobians */
  free(jac);
  free(pjac);
  free(dcovar);

  return;
  }


/****** subprofit_convmoments *************************************************
PROTO	void profit_convmoments(subprofitstruct *subprofit, obj2struct *obj2)
PURPOSE	Compute the 2nd order moments of the convolved object model.
INPUT	Sub-profile-fitting structure,
	Pointer to obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/05/2012
 ***/
void	 subprofit_convmoments(subprofitstruct *subprofit, obj2struct *obj2)
  {
   double	hw,hh, r2max, x,xstart,y, mx2,my2,mxy,mx,my,sum, dval,
		temp,temp2,invtemp2, pmx2, theta;
   PIXTYPE	*pix;
   int		ix,iy, w,h;

  w = subprofit->modnaxisn[0];
  h = subprofit->modnaxisn[1];
  hw = (double)(w/2);
  hh = (double)(h/2);

  r2max = hw<hh? hw*hw : hh*hh;
  xstart = -hw;
  y = -hh;
  pix = subprofit->cmodpix;
  mx2 = my2 = mxy = mx = my = sum = 0.0;
  for (iy=h; iy--; y+=1.0)
    {
    x = xstart;
    for (ix=w; ix--; x+=1.0)
      if (y*y+x*x <= r2max)
        {
        dval = *(pix++);
        sum += dval;
        mx  += dval*x;
        my  += dval*y;
        mx2 += dval*x*x;
        mxy += dval*x*y;
        my2 += dval*y*y;
        }
      else
        pix++;
    }

  if (sum <= 1.0/BIG)
    sum = 1.0;
  mx /= sum;
  my /= sum;
  obj2->prof_convmx2 = (mx2 = mx2/sum - mx*mx)
			*subprofit->pixstep*subprofit->pixstep;
  obj2->prof_convmy2 = (my2 = my2/sum - my*my)
			*subprofit->pixstep*subprofit->pixstep;
  obj2->prof_convmxy = (mxy = mxy/sum - mx*my)
			*subprofit->pixstep*subprofit->pixstep;

/* Handle fully correlated profiles (which cause a singularity...) */
  if ((temp2=mx2*my2-mxy*mxy)<0.00694)
    {
    mx2 += 0.0833333;
    my2 += 0.0833333;
    temp2 = mx2*my2-mxy*mxy;
    }

  temp2 *= subprofit->pixstep*subprofit->pixstep;

  if (FLAG(obj2.prof_convcxx))
    {
    invtemp2 = (temp2>=0.0) ? 1.0/temp2 : 0.0;
    obj2->prof_convcxx = (float)(my2*invtemp2);
    obj2->prof_convcyy = (float)(mx2*invtemp2);
    obj2->prof_convcxy = (float)(-2*mxy*invtemp2);
    }

  if (1 /*FLAG(obj2.prof_conva)*/)
    {
    if ((fabs(temp=mx2-my2)) > 0.0)
      theta = atan2(2.0 * mxy,temp) / 2.0;
    else
      theta = PI/4.0;

    temp = sqrt(0.25*temp*temp+mxy*mxy);
    pmx2 = 0.5*(mx2+my2);
    obj2->prof_conva = (float)sqrt(pmx2 + temp)*subprofit->pixstep;
    obj2->prof_convb = (float)sqrt(pmx2 - temp)*subprofit->pixstep;
    obj2->prof_convtheta = theta/DEG;
    }

  return;
  }


/****** subprofit_surface ****************************************************
PROTO	void subprofit_surface(profitstruct *profit, subprofitstruct *subprofit,
			obj2struct *obj2)
PURPOSE	Compute surface brightnesses from the unconvolved object model.
INPUT	Pointer to the profile-fitting structure,
	Pointer to the sub-profile-fitting structure,
	Pointer to obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
void	 subprofit_surface(profitstruct *profit, subprofitstruct *subprofit,
			obj2struct *obj2)
  {
   subprofitstruct	hdsubprofit;
   double		dsum,dhsum,dsumoff, dhval, frac, seff;
   float		*spix, *spixt,
			val,vmax,
			scalefac, imsizefac, flux, lost, sum, lostfluxfrac;
   int			i,p, imax, npix, neff;

/* Allocate "high-definition" raster only to make measurements */
  hdsubprofit.modnaxisn[0] = hdsubprofit.modnaxisn[1] = PROFIT_HIDEFRES;
  npix = hdsubprofit.nmodpix
	= hdsubprofit.modnaxisn[0]*hdsubprofit.modnaxisn[1];
/* Find best image size factor from fitting results */
  imsizefac = 2.0*profit_minradius(profit, PROFIT_REFFFAC)/subprofit->pixstep
	/ (float)subprofit->modnaxisn[0];
  if (imsizefac<0.01)
    imsizefac = 0.01;
  else if (imsizefac>100.0)
    imsizefac = 100.0;
  scalefac = (float)hdsubprofit.modnaxisn[0] / (float)subprofit->modnaxisn[0]
	/ imsizefac;
  hdsubprofit.pixstep = subprofit->pixstep / scalefac;
  hdsubprofit.fluxfac = 1.0/(hdsubprofit.pixstep*hdsubprofit.pixstep);
  QCALLOC(hdsubprofit.modpix, float,npix*sizeof(float));

  for (p=0; p<profit->nparam; p++)
    profit->param[p] = profit->paraminit[p];
  lost = sum = 0.0;

  for (p=0; p<profit->nprof; p++)
    {
    sum += (flux = prof_add(&hdsubprofit, profit->prof[p],0));
    lost += flux*profit->prof[p]->lostfluxfrac;
    }
  lostfluxfrac = sum > 0.0? lost / sum : 0.0;
/*
char filename[256];
checkstruct *check;
sprintf(filename, "raster_%02d.fits", the_gal);
check=check_init(filename, CHECK_OTHER, 0);
check->width = hdprofit.modnaxisn[0];
check->height = hdprofit.modnaxisn[1];
check_reinit(the_field, check);
memcpy(check->pix,hdprofit.modpix,check->npix*sizeof(float));
check_reend(the_field, check);
check_end(check);
*/
  if (FLAG(obj2.fluxeff_prof))
    {
/*-- Sort model pixel values */
    spix = NULL;			/* to avoid gcc -Wall warnings */
    QMEMCPY(hdsubprofit.modpix, spix, float, npix);
    fqmedian(spix, npix);
/*-- Build a cumulative distribution */
    dsum = 0.0;
    spixt = spix;
    for (i=npix; i--;)
      dsum += (double)*(spixt++);
/*-- Find matching surface brightness */
    if (lostfluxfrac > 1.0)
      lostfluxfrac = 0.0;
    dhsum = 0.5 * dsum / (1.0-lostfluxfrac);
    dsum = lostfluxfrac * dsum / (1.0-lostfluxfrac);
    neff = 0;
    spixt = spix;
    for (i=npix; i--;)
      if ((dsum += (double)*(spixt++)) >= dhsum)
        {
        neff = i;
        break;
        }
    dhval = (double)*(spixt-1);
    seff = neff;
    dsumoff = 0.0;
    if (spixt>=spix+2)
      if (dhval > 0.0 && (frac = (dsum - dhsum) / dhval) < 1.0)
        {
        seff += frac;
        dsumoff = frac*dhval;
        dhval = dsumoff + (1.0 - frac)*(double)*(spixt-2);
        }
    obj2->fluxeff_prof = dhval;
    if (FLAG(obj2.fluxmean_prof))
      {
      dsum = dsumoff;
      for (i=neff; i--;)
        dsum += (double)*(spixt++);
      obj2->fluxmean_prof = seff > 0.0? dsum / seff : 0.0;
      }
    free(spix);
    }

/* Compute model peak */
  if (FLAG(obj2.peak_prof))
    {
/*-- Find position of maximum pixel in current hi-def raster */
    imax = 0;
    vmax = -BIG;
    spixt = hdsubprofit.modpix;
    for (i=npix; i--;)
      if ((val=*(spixt++))>vmax)
        {
        vmax = val;
        imax = i;
        }
    imax = npix-1 - imax;
    obj2->peak_prof = hdsubprofit.modpix[imax];
    }

/* Free hi-def model raster */
  free(hdsubprofit.modpix);

  return;
  }


/****** profit_addparam *******************************************************
PROTO	void profit_addparam(profitstruct *profit, paramenum paramindex,
		float **param)
PURPOSE	Add a profile parameter to the list of fitted items.
INPUT	Pointer to the profit structure,
	Parameter index,
	Pointer to the parameter pointer.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/03/2010
 ***/
void	profit_addparam(profitstruct *profit, paramenum paramindex,
		float **param)
  {
/* Check whether the parameter has already be registered */
  if (profit->paramlist[paramindex])
/*-- Yes */
    *param = profit->paramlist[paramindex];
  else
/*-- No */
    {
    *param = profit->paramlist[paramindex] = &profit->param[profit->nparam];
    profit->paramindex[paramindex] = profit->nparam;
    profit->paramrevindex[profit->nparam++] = paramindex;
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
VERSION	24/04/2013
 ***/
void	profit_resetparam(profitstruct *profit, paramenum paramtype)
  {
   subprofitstruct	*subprofit;
   float		param, parammin,parammax, range;
   parfitenum		fittype;

  subprofit = profit->subprofit;	/* Take the first subprofit !CHECK */
  param = parammin = parammax = 0.0;	/* Avoid gcc -Wall warnings*/

  if (!profit->paramlist[(int)paramtype])
    return;

  switch(paramtype)
    {
    case PARAM_BACK:
      fittype = PARFIT_LINBOUND;
      param = 0.0;
      parammin = -6.0*subprofit->guesssigbkg;
      parammax =  6.0*subprofit->guesssigbkg;
      break;
    case PARAM_X:
      fittype = PARFIT_LINBOUND;
      param = subprofit->guessdx;
      range = subprofit->guessradius*4.0;
      if (range>subprofit->objnaxisn[0]*subprofit->subsamp*2.0)
        range = subprofit->objnaxisn[0]*subprofit->subsamp*2.0;
      parammin = -range;
      parammax =  range;
      break;
    case PARAM_Y:
      fittype = PARFIT_LINBOUND;
      param = subprofit->guessdy;
      range = subprofit->guessradius*4.0;
      if (range>subprofit->objnaxisn[1]*subprofit->subsamp*2.0)
        range = subprofit->objnaxisn[1]*subprofit->subsamp*2.0;
      parammin = -range;
      parammax =  range;
      break;
    case PARAM_DIRAC_FLUX:
    case PARAM_DIRAC_FLUX2:
    case PARAM_DIRAC_FLUX3:
    case PARAM_DIRAC_FLUX4:
    case PARAM_DIRAC_FLUX5:
      fittype = PARFIT_LOGBOUND;
      param = subprofit[paramtype-PARAM_DIRAC_FLUX].guessflux/profit->nprof;
      parammin = 0.00001*subprofit[paramtype-PARAM_DIRAC_FLUX].guessfluxmax;
      parammax = 10.0*subprofit[paramtype-PARAM_DIRAC_FLUX].guessfluxmax;
      break;
    case PARAM_SPHEROID_FLUX:
    case PARAM_SPHEROID_FLUX2:
    case PARAM_SPHEROID_FLUX3:
    case PARAM_SPHEROID_FLUX4:
    case PARAM_SPHEROID_FLUX5:
      fittype = PARFIT_LOGBOUND;
      param = subprofit[paramtype-PARAM_SPHEROID_FLUX].guessflux/profit->nprof;
      parammin = 0.00001*subprofit[paramtype-PARAM_SPHEROID_FLUX].guessfluxmax;
      parammax = 10.0*subprofit[paramtype-PARAM_SPHEROID_FLUX].guessfluxmax;
      break;
    case PARAM_SPHEROID_REFF:
      fittype = PARFIT_LOGBOUND;
      param = FLAG(obj2.prof_disk_flux)? subprofit->guessradius
			: subprofit->guessradius/sqrtf(subprofit->guessaspect);
      parammin = 0.01;
      parammax = param * 10.0;
      break;
    case PARAM_SPHEROID_ASPECT:
      fittype = PARFIT_LOGBOUND;
      param = FLAG(obj2.prof_disk_flux)? 1.0 : subprofit->guessaspect;
      parammin = FLAG(obj2.prof_disk_flux)? 0.5 : 0.01;
      parammax = FLAG(obj2.prof_disk_flux)? 2.0 : 100.0;
      break;
    case PARAM_SPHEROID_POSANG:
      fittype = PARFIT_UNBOUND;
      param = subprofit->guessposang;
      parammin = 90.0;
      parammax =  90.0;
      break;
    case PARAM_SPHEROID_SERSICN:
      fittype = PARFIT_LINBOUND;
      param = 4.0;
      parammin = FLAG(obj2.prof_disk_flux)? 1.0 : 0.3;
      parammax = 10.0;
      break;
    case PARAM_DISK_FLUX:
    case PARAM_DISK_FLUX2:
    case PARAM_DISK_FLUX3:
    case PARAM_DISK_FLUX4:
    case PARAM_DISK_FLUX5:
      fittype = PARFIT_LOGBOUND;
      param = subprofit[paramtype-PARAM_DISK_FLUX].guessflux/profit->nprof;
      parammin = 0.00001*subprofit[paramtype-PARAM_DISK_FLUX].guessfluxmax;
      parammax = 10.0*subprofit[paramtype-PARAM_DISK_FLUX].guessfluxmax;
      break;
    case PARAM_DISK_SCALE:	/* From scalelength to Re */
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessradius/(1.67835*sqrtf(subprofit->guessaspect));
      parammin = 0.01/1.67835;
      parammax = param * 10.0;
      break;
    case PARAM_DISK_ASPECT:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessaspect;
      parammin = 0.01;
      parammax = 100.0;
      break;
    case PARAM_DISK_POSANG:
      fittype = PARFIT_UNBOUND;
      param = subprofit->guessposang;
      parammin = 90.0;
      parammax =  90.0;
      break;
    case PARAM_ARMS_FLUX:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessflux/profit->nprof;
      parammin = 0.00001*subprofit->guessfluxmax;
      parammax = 10.0*subprofit->guessfluxmax;
      break;
    case PARAM_ARMS_QUADFRAC:
      fittype = PARFIT_LINBOUND;
      param = 0.5;
      parammin = 0.0;
      parammax = 1.0;
      break;
    case PARAM_ARMS_SCALE:
      fittype = PARFIT_LINBOUND;
      param = 1.0;
      parammin = 0.5;
      parammax = 10.0;
      break;
    case PARAM_ARMS_START:
      fittype = PARFIT_LINBOUND;
      param = 0.5;
      parammin = 0.0;
      parammax = 3.0;
      break;
    case PARAM_ARMS_PITCH:
      fittype = PARFIT_LINBOUND;
      param = 20.0;
      parammin = 5.0;
      parammax = 50.0;
      break;
    case PARAM_ARMS_PITCHVAR:
      fittype = PARFIT_LINBOUND;
      param = 0.0;
      parammin = -1.0;
      parammax = 1.0;
      break;
//      if ((profit->spirindex=profit_spiralindex(profit, obj2)) > 0.0)
//        {
//        param = -param;
//        parammin = -parammax;
//        parammax = -parammin;
//        }
//      printf("spiral index: %g  \n", profit->spirindex);
//      break;
    case PARAM_ARMS_POSANG:
      fittype = PARFIT_UNBOUND;
      param = 0.0;
      parammin = 0.0;
      parammax = 0.0;
      break;
    case PARAM_ARMS_WIDTH:
      fittype = PARFIT_LINBOUND;
      param = 3.0;
      parammin = 1.5;
      parammax = 11.0;
      break;
    case PARAM_BAR_FLUX:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessflux/profit->nprof;
      parammin = 0.00001*subprofit->guessfluxmax;
      parammax = 4.0*subprofit->guessfluxmax;
      break;
    case PARAM_BAR_ASPECT:
      fittype = PARFIT_LOGBOUND;
      param = 0.3;
      parammin = 0.2;
      parammax = 0.5;
      break;
    case PARAM_BAR_POSANG:
      fittype = PARFIT_UNBOUND;
      param = 0.0;
      parammin = 90.0;
      parammax = 90.0;
      break;
    case PARAM_INRING_FLUX:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessflux/profit->nprof;
      parammin = 0.00001*subprofit->guessfluxmax;
      parammax = 4.0*subprofit->guessfluxmax;
      break;
    case PARAM_INRING_WIDTH:
      fittype = PARFIT_LINBOUND;
      param = 0.3;
      parammin = 0.0;
      parammax = 0.5;
      break;
    case PARAM_INRING_ASPECT:
      fittype = PARFIT_LOGBOUND;
      param = 0.8;
      parammin = 0.4;
      parammax = 1.0;
      break;
    case PARAM_OUTRING_FLUX:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessflux/profit->nprof;
      parammin = 0.00001*subprofit->guessfluxmax;
      parammax = 4.0*subprofit->guessfluxmax;
      break;
    case PARAM_OUTRING_START:
      fittype = PARFIT_LINBOUND;
      param = 4.0;
      parammin = 3.5;
      parammax = 6.0;
      break;
    case PARAM_OUTRING_WIDTH:
      fittype = PARFIT_LINBOUND;
      param = 0.3;
      parammin = 0.0;
      parammax = 0.5;
      break;
    case PARAM_MOFFAT_FLUX:
      fittype = PARFIT_LOGBOUND;
      param = subprofit[paramtype-PARAM_MOFFAT_FLUX].guessflux/profit->nprof;
      parammin = 0.00001*subprofit[paramtype-PARAM_MOFFAT_FLUX].guessfluxmax;
      parammax = 100.0*subprofit[paramtype-PARAM_MOFFAT_FLUX].guessfluxmax;
      break;
    case PARAM_MOFFAT_ALPHA:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessradius/sqrtf(subprofit->guessaspect) * 2.0;
		/* Approximatively 1/sqrt(2^(1/beta)-1) with beta = 4 */
      parammin = 0.1;
      parammax = param * 10.0;
      break;
    case PARAM_MOFFAT_ASPECT:
      fittype = PARFIT_LOGBOUND;
      param = subprofit->guessaspect;
      parammin = 0.001;
      parammax = 1000.0;
      break;
    case PARAM_MOFFAT_POSANG:
      fittype = PARFIT_UNBOUND;
      param = subprofit->guessposang;
      parammin = 90.0;
      parammax =  90.0;
      break;
    case PARAM_MOFFAT_BETA:
      fittype = PARFIT_LINBOUND;
      param = 4.0;
      parammin = 0.5;
      parammax = 8.0;
      break;
    case PARAM_MOFFAT_MINKP:
      fittype = PARFIT_LINBOUND;
      param = 2.0;
      parammin = 1.9;
      parammax = 8.0;
      break;
    case PARAM_MOFFAT_OFFSET:
      fittype = PARFIT_LOGBOUND;
      param = 0.001;
      parammin = 0.000001;
      parammax = 5.0;
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown profile parameter in ",
		"profit_resetparam()");
      break;
   }

  if (parammin!=parammax && (param<=parammin || param>=parammax))
    param = (parammin+parammax)/2.0;
  else if (parammin==0.0 && parammax==0.0)
    parammax = 1.0;
  profit_setparam(profit, paramtype, param, parammin, parammax, fittype);

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
		float param, float parammin, float parammax,
		parfitenum parfittype)
PURPOSE	Set the actual, lower and upper boundary values of a profile parameter.
INPUT	Pointer to the profit structure,
	parameter index,
	actual value,
	lower boundary to the parameter,
	upper boundary to the parameter,
	parameter fitting type.
OUTPUT	RETURN_OK if the parameter is registered, RETURN_ERROR otherwise.
AUTHOR	E. Bertin (IAP)
VERSION	20/05/2011
 ***/
int	profit_setparam(profitstruct *profit, paramenum paramtype,
		float param, float parammin,float parammax,
		parfitenum parfittype)
  {
   float	*paramptr;
   int		index;

/* Check whether the parameter has already be registered */
  if ((paramptr=profit->paramlist[(int)paramtype]))
    {
    index = profit->paramindex[(int)paramtype];
    profit->paraminit[index] = param;
    profit->parammin[index] = parammin;
    profit->parammax[index] = parammax;
    profit->parfittype[index] = parfittype;
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }

  
/****** profit_boundtounbound *************************************************
PROTO	int profit_boundtounbound(profitstruct *profit,
				float *param, double *dparam, int index)
PURPOSE	Convert parameters from bounded to unbounded space.
INPUT	Pointer to the profit structure,
	input array of single precision parameters,
	output (incomplete) array of double precision parameters,
	parameter selection index (<0 = all parameters)
OUTPUT	Number of free parameters.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/05/2011
 ***/
int	profit_boundtounbound(profitstruct *profit,
				float *param, double *dparam, int index)
  {
   double	num,den, tparam;
   int		f,p, pstart,np;

  if (index<0)
    {
    pstart = 0;
    np = profit->nparam;
    }
  else
    {
    pstart = index;
    np = pstart+1;
    }

  f = 0;
  for (p=pstart ; p<np; p++)
    switch(profit->parfittype[p])
      {
      case PARFIT_FIXED:
        break;
      case PARFIT_UNBOUND:
        if (profit->parammax[p] > 0.0 || profit->parammax[p] < 0.0)
          dparam[f] = param[p-pstart] / profit->parammax[p];
        f++;
        break;
      case PARFIT_LINBOUND:
        tparam = param[p-pstart];
        num = tparam - profit->parammin[p];
        den = profit->parammax[p] - tparam;
        dparam[f] = num>1e-50? (den>1e-50? log(num/den): 50.0) : -50.0;
        f++;
        break;
      case PARFIT_LOGBOUND:
        tparam = log(param[p-pstart]);
        num = tparam - log(profit->parammin[p]);
        den = log(profit->parammax[p]) - tparam;
        dparam[f] = num>1e-50? (den>1e-50? log(num/den): 50.0) : -50.0;
        f++;
        break;
      default:
        error(EXIT_FAILURE,
		"*Internal Error*: Unknown parameter fitting type in ",
		"profit_boundtounbound()");
      }

  return f;
  }


/****** profit_unboundtobound *************************************************
PROTO	int profit_unboundtobound(profitstruct *profit,
				double *dparam, float *param, int index)
PURPOSE	Convert parameters from unbounded to bounded space.
INPUT	Pointer to the profit structure,
	input (incomplete) array of double precision parameters,
	output array of single precision parameters.
	parameter selection index (<0 = all parameters)
OUTPUT	Number of free parameters.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/05/2011
 ***/
int	profit_unboundtobound(profitstruct *profit,
				double *dparam, float *param, int index)
  {
   int		f,p, pstart,np;

  if (index<0)
    {
    pstart = 0;
    np = profit->nparam;
    }
  else
    {
    pstart = index;
    np = pstart+1;
    }

  f = 0;
  for (p=pstart; p<np; p++)
    switch(profit->parfittype[p])
      {
      case PARFIT_FIXED:
        param[p-pstart] = profit->paraminit[p];
        break;
      case PARFIT_UNBOUND:
        param[p-pstart] = dparam[f]*profit->parammax[p];
        f++;
        break;
      case PARFIT_LINBOUND:
        param[p-pstart] = (profit->parammax[p] - profit->parammin[p])
		/ (1.0 + exp(-(dparam[f]>50.0? 50.0
				: (dparam[f]<-50.0? -50.0: dparam[f]))))
		+ profit->parammin[p];
        f++;
        break;
      case PARFIT_LOGBOUND:
        param[p-pstart] = exp(log(profit->parammax[p]/profit->parammin[p])
		/ (1.0 + exp(-(dparam[f]>50.0? 50.0
				: (dparam[f]<-50.0? -50.0: dparam[f])))))
		*profit->parammin[p];
        f++;
        break;
      default:
        error(EXIT_FAILURE,
		"*Internal Error*: Unknown parameter fitting type in ",
		"profit_unboundtobound()");
      }
 
  return f;
  }


/****** profit_covarunboundtobound ********************************************
PROTO	int profit_covarunboundtobound(profitstruct *profit,
				double *dcovar, float *covar)
PURPOSE	Convert covariance matrix from unbounded to bounded space.
INPUT	Pointer to the profit structure,
	input (incomplete) matrix of double precision covariances,
	output matrix of single precision covariances.
OUTPUT	Number of parameters that hit a boundary.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/05/2011
 ***/
int	profit_covarunboundtobound(profitstruct *profit,
				double *dcovar, float *covar)
  {
   double	*dxdy,
		xmmin,maxmx, maxmmin;
   float	*x,*xmin,*xmax;
   parfitenum	*fittype;
   int		*fflag,
		f,f1,f2, p,p1,p2, nfree, nparam, nmin,nmax;

  nparam = profit->nparam;
  fittype = profit->parfittype;
  QMALLOC16(dxdy, double, PARAM_NPARAM);
  x = profit->paraminit;
  xmin = profit->parammin;
  xmax = profit->parammax;
  nmin = nmax = 0;
  f = 0;
  for (p=0; p<profit->nparam; p++)
    switch(fittype[p])
      {
      case PARFIT_FIXED:
        break;
      case PARFIT_UNBOUND:
        dxdy[f++] = xmax[p];
        break;
      case PARFIT_LINBOUND:
        xmmin   = x[p] - xmin[p];
        maxmx   = xmax[p] - x[p];
        maxmmin = xmax[p] - xmin[p];
        dxdy[f++] = xmmin * maxmx / maxmmin;
        if (xmmin<0.001*maxmmin)
          nmin++;
        else if (maxmx<0.001*maxmmin)
          nmax++;
        break;
      case PARFIT_LOGBOUND:
        xmmin   = log(x[p]/xmin[p]);
        maxmx   = log(xmax[p]/x[p]);
        maxmmin = log(xmax[p]/xmin[p]);
        dxdy[f++] = x[p] * xmmin * maxmx / maxmmin;
        if (xmmin<0.001*maxmmin)
          nmin++;
        else if (maxmx<0.001*maxmmin)
          nmax++;
        break;
      default:
        error(EXIT_FAILURE,
		"*Internal Error*: Unknown parameter fitting type in ",
		"profit_boundtounbound()");
      }

  nfree = f;

  memset(profit->covar, 0, nparam*nparam*sizeof(float));
  f2 = 0;
  for (p2=0; p2<nparam; p2++)
    {
    if (fittype[p2]!=PARFIT_FIXED)
      {
      f1 = 0;
      for (p1=0; p1<nparam; p1++)
        if (fittype[p1]!=PARFIT_FIXED)
          {
          covar[p2*nparam+p1] = (float)(dcovar[f2*nfree+f1]*dxdy[f1]*dxdy[f2]);
          f1++;
          }
      f2++;
      }
    }

  free(dxdy);

  profit->nlimmin = nmin;
  profit->nlimmax = nmax;

  return nmin+nmax;
  }


/****** prof_init *************************************************************
PROTO	profstruct prof_init(profitstruct *profit, unsigned int modeltype)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	Pointer to the profile-fitting structure,
	model type.
OUTPUT	A pointer to an allocated prof structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	30/12/2013
 ***/
profstruct	*prof_init(profitstruct *profit, unsigned int modeltype)
  {
   profstruct	*prof;
   float	*pix,
		rmax2, re2, dy2,r2, scale, zero, k,n, hinvn;
   int		width,height, ixc,iyc, ix,iy, nsub,
		d,s;

  QCALLOC(prof, profstruct, 1);
  prof->code = modeltype;
  switch(modeltype)
    {
    case MODEL_BACK:
      prof->name = "background offset";
      prof->naxis = 2;
      prof->pix = NULL;
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, PARAM_BACK, &prof->flux[s]);
      prof->typscale = 1.0;
      break;
    case MODEL_DIRAC:
      prof->name = "point source";
      prof->naxis = 2;
/*
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
*/
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_DIRAC_FLUX+s), &prof->flux[s]);
      break;
    case MODEL_SERSIC:
      prof->name = "Sersic spheroid";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_SPHEROID_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_SPHEROID_REFF, &prof->scale);
      profit_addparam(profit, PARAM_SPHEROID_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SPHEROID_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_SPHEROID_SERSICN, &prof->extra[0]);
      break;
    case MODEL_DEVAUCOULEURS:
      prof->name = "de Vaucouleurs spheroid";
      prof->naxis = 2;
/*
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
*/
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_SPHEROID_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_SPHEROID_REFF, &prof->scale);
      profit_addparam(profit, PARAM_SPHEROID_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_SPHEROID_POSANG, &prof->posangle);
      break;
    case MODEL_EXPONENTIAL:
      prof->name = "exponential disk";
      prof->naxis = 2;
/*
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
*/
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_DISK_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      break;
    case MODEL_ARMS:
      prof->name = "spiral arms";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_ARMS_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_ARMS_QUADFRAC, &prof->featfrac);
//      profit_addparam(profit, PARAM_ARMS_SCALE, &prof->featscale);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      profit_addparam(profit, PARAM_ARMS_PITCH, &prof->featpitch);
//      profit_addparam(profit, PARAM_ARMS_PITCHVAR, &prof->featpitchvar);
      profit_addparam(profit, PARAM_ARMS_POSANG, &prof->featposang);
//      profit_addparam(profit, PARAM_ARMS_WIDTH, &prof->featwidth);
      break;
    case MODEL_BAR:
      prof->name = "bar";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_BAR_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_BAR_ASPECT, &prof->feataspect);
      profit_addparam(profit, PARAM_ARMS_POSANG, &prof->featposang);
      break;
    case MODEL_INRING:
      prof->name = "inner ring";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_ARMS_START, &prof->featstart);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_INRING_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_INRING_WIDTH, &prof->featwidth);
      profit_addparam(profit, PARAM_INRING_ASPECT, &prof->feataspect);
      break;
    case MODEL_OUTRING:
      prof->name = "outer ring";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      profit_addparam(profit, PARAM_DISK_SCALE, &prof->scale);
      profit_addparam(profit, PARAM_DISK_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_DISK_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_OUTRING_START, &prof->featstart);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_OUTRING_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_OUTRING_WIDTH, &prof->featwidth);
      break;
    case MODEL_MOFFAT:
      prof->name = "Generalized Moffat profile";
      prof->naxis = 2;
      prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      prof->naxisn[2] = 1;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      prof->typscale = 1.0;
      profit_addparam(profit, PARAM_X, &prof->x[0]);
      profit_addparam(profit, PARAM_Y, &prof->x[1]);
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_MOFFAT_FLUX+s),
		&prof->flux[s]);
      profit_addparam(profit, PARAM_MOFFAT_ALPHA, &prof->scale);
      profit_addparam(profit, PARAM_MOFFAT_ASPECT, &prof->aspect);
      profit_addparam(profit, PARAM_MOFFAT_POSANG, &prof->posangle);
      profit_addparam(profit, PARAM_MOFFAT_BETA, &prof->extra[0]);
      profit_addparam(profit, PARAM_MOFFAT_MINKP, &prof->extra[1]);
      profit_addparam(profit, PARAM_MOFFAT_OFFSET, &prof->extra[2]);
      break;
    case MODEL_TABULATED:	/* An example of tabulated profile */
      prof->name = "tabulated model";
      prof->naxis = 3;
      width =  prof->naxisn[0] = profit->subprofit->modnaxisn[0];
      height = prof->naxisn[1] = profit->subprofit->modnaxisn[1];
      nsub = prof->naxisn[2] = PROFIT_MAXSMODSIZE;
      prof->npix = prof->naxisn[0]*prof->naxisn[1]*prof->naxisn[2];
      QMALLOC(prof->pix, float, prof->npix);
      ixc = width/2;
      iyc = height/2;
      rmax2 = (ixc - 1.0)*(ixc - 1.0);
      re2 = width/64.0;
      prof->typscale = re2;
      re2 *= re2;
      zero = prof->extrazero[0] = 0.2;
      scale = prof->extrascale[0]= 8.0/PROFIT_MAXSMODSIZE;
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
      for (s=0; s<profit->nsubprofit; s++)
        profit_addparam(profit, (paramenum)((int)PARAM_SPHEROID_FLUX+s),
		&prof->flux[s]);
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

  if (prof->naxis>2)
    {
    prof->kernelnlines = 1;
    for (d=0; d<prof->naxis; d++)
      {
      prof->interptype[d] = INTERP_BILINEAR;
      prof->kernelnlines *= 
	(prof->kernelwidth[d] = interp_kernwidth[prof->interptype[d]]);
      }
    prof->kernelnlines /= prof->kernelwidth[0];
    QMALLOC16(prof->kernelbuf, float, prof->kernelnlines);
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
PROTO	float prof_add(subprofitstruct *subprofit, profstruct *prof,
		int extfluxfac_flag)
PURPOSE	Add a model profile to an image.
INPUT	Sub-profile-fitting structure,
	profile structure,
	flag (0 if flux correction factor is to be computed internally) 
OUTPUT	Total (asymptotic) flux contribution.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	30/12/2013
 ***/
float	prof_add(subprofitstruct *subprofit, profstruct *prof,
		int extfluxfac_flag)
  {
   double	xscale, yscale, saspect, ctheta,stheta, flux, scaling, bn, n,
		dx1cout,dx2cout, ddx1[36],ddx2[36];
   float	posin[PROFIT_MAXEXTRA], posout[2], dnaxisn[2],
		*pixin, *pixin2, *pixout,
		pflux, fluxfac, amp,cd11,cd12,cd21,cd22, dx1,dx2,
		x1,x10,x2, x1cin,x2cin, x1cout,x2cout, x1max,x2max, x1in,x2in,
		k, hinvn, x1t,x2t, ca,sa, u,umin,
		armamp,arm2amp, armrdphidr, armrdphidrvar, posang,
		width, invwidth2,
		r,r2,rmin, r2minxin,r2minxout, rmax, r2max,
		r2max1, r2max2, r2min, invr2xdif,
		val, theta, thresh, ra,rb, num,num2,den, ang,angstep,
		invn, dr, krpinvn,dkrpinvn, rs,rs2,
		a11,a12,a21,a22, invdet, dca,dsa, a0,a2,a3, p1,p2, off,
		krspinvn, ekrspinvn, selem, mofbeta,mofp,invmofalpha2,invmofp2;
   int		npix, threshflag,
		a,d,e,i,s, ix1,ix2, ix1max,ix2max, nang, nx2,
		npix2;

  npix = subprofit->nmodpix;
  pflux = *prof->flux[subprofit->index];

  if (prof->code==MODEL_BACK)
    {
    pixout = subprofit->modpix;
    for (i=npix; i--;)
      *(pixout++) += pflux;
    prof->lostfluxfrac = 0.0;
    return 0.0;
    }

  scaling = subprofit->pixstep / prof->typscale;
  QMALLOC(prof->pix, float, subprofit->nmodpix);

  if (prof->code!=MODEL_DIRAC)
    {
/*-- Compute Profile CD matrix */
    ctheta = cos(*prof->posangle*DEG);
    stheta = sin(*prof->posangle*DEG);
    saspect = fabs(*prof->aspect);
    xscale = (*prof->scale==0.0)?
		0.0 : fabs(scaling / (*prof->scale*prof->typscale));
    yscale = (*prof->scale*saspect == 0.0)?
		0.0 : fabs(scaling / (*prof->scale*prof->typscale*saspect));
    cd11 = (float)(xscale*ctheta);
    cd12 = (float)(xscale*stheta);
    cd21 = (float)(-yscale*stheta);
    cd22 = (float)(yscale*ctheta);
    dx1 = 0.0;	/* Shifting operations have been moved to profit_resample() */
    dx2 = 0.0;	/* Shifting operations have been moved to profit_resample() */

    x1cout = (float)(subprofit->modnaxisn[0]/2);
    x2cout = (float)(subprofit->modnaxisn[1]/2);
    nx2 = subprofit->modnaxisn[1]/2 + 1;

/*-- Compute the largest r^2 that fits in the frame */
    num = cd11*cd22-cd12*cd21;
    num *= num;
    x1max = x1cout - 1.0;
    x2max = x2cout - 1.0;
    den = fabs(cd12*cd12+cd22*cd22);
    num2 = x1max*x1max*num;
    r2max1 = num2<PROFIT_MAXR2MAX*den? num2 / den : PROFIT_MAXR2MAX;
    den = fabs(cd11*cd11+cd21*cd21);
    num2 = x2max*x2max*num;
    r2max2 = num2<PROFIT_MAXR2MAX*den? num2 / den : PROFIT_MAXR2MAX;
    r2max = (r2max1 < r2max2? r2max1 : r2max2);
    rmax = sqrtf(r2max);
    }

  switch(prof->code)
    {
    case MODEL_DIRAC:
      prof->pix[subprofit->modnaxisn[0]/2
		+ (subprofit->modnaxisn[1]/2)*subprofit->modnaxisn[0]] = 1.0;
      prof->lostfluxfrac = 0.0;
      threshflag = 0;
      break;
    case MODEL_SERSIC:
    case MODEL_DEVAUCOULEURS:
    case MODEL_EXPONENTIAL:
/*---- Compute sharp/smooth transition radius */
      rs = PROFIT_SMOOTHR*(xscale>yscale?xscale:yscale);
      if (rs<=0)
        rs = 1.0;
      rs2 = rs*rs;
/*---- The consequence of sampling on flux is compensated by PSF normalisation*/
      if (prof->code==MODEL_EXPONENTIAL)
        bn = n = 1.0;
      else if (prof->code==MODEL_DEVAUCOULEURS)
        {
        n = 4.0;
        bn = 7.66924944;
        }
      else
        {
        n = fabs(*prof->extra[0]);
        bn = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
		+ 131.0/(1148175*n*n*n);	/* Ciotti & Bertin 1999 */
        }
      invn = 1.0/n;
      hinvn = 0.5/n;
      k = -bn;
/*---- Compute central polynomial terms */
      krspinvn = prof->code==MODEL_EXPONENTIAL? -rs : k*expf(logf(rs)*invn);
      ekrspinvn = expf(krspinvn);
      p2 = krspinvn*invn*invn;
      p1 = krspinvn*p2;
      a0 = (1+(1.0/6.0)*(p1+(1.0-5.0*n)*p2))*ekrspinvn;
      a2 = (-1.0/2.0)*(p1+(1.0-3.0*n)*p2)/rs2*ekrspinvn;
      a3 = (1.0/3.0)*(p1+(1.0-2.0*n)*p2)/(rs2*rs)*ekrspinvn;
/*---- Compute the smooth part of the profile */
      x10 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = prof->pix;
      if (prof->code==MODEL_EXPONENTIAL)
        for (ix2=nx2; ix2--; x2+=1.0)
          {
          x1 = x10;
          for (ix1=subprofit->modnaxisn[0]; ix1--; x1+=1.0)
            {
            x1in = cd12*x2 + cd11*x1;
            x2in = cd22*x2 + cd21*x1;
            ra = x1in*x1in+x2in*x2in;
            if (ra>r2max)
              {
              *(pixin++) = 0.0;
              continue;
              }
            val = ra<rs2? a0+ra*(a2+a3*sqrtf(ra)) : expf(-sqrtf(ra));
            *(pixin++) = val;
            }
          }
      else
        for (ix2=nx2; ix2--; x2+=1.0)
          {
          x1 = x10;
          for (ix1=subprofit->modnaxisn[0]; ix1--; x1+=1.0)
            {
            x1in = cd12*x2 + cd11*x1;
            x2in = cd22*x2 + cd21*x1;
            ra = x1in*x1in+x2in*x2in;
            if (ra>r2max)
              {
              *(pixin++) = 0.0;
              continue;
              }
            val = ra<rs2? a0+ra*(a2+a3*sqrtf(ra)) : expf(k*expf(logf(ra)*hinvn));
            *(pixin++) = val;
            }
          }
/*---- Copy the symmetric part */
      if ((npix2=(subprofit->modnaxisn[1]-nx2)*subprofit->modnaxisn[0]) > 0)
        {
        pixin2 = pixin - subprofit->modnaxisn[0] - 1;
        if (!(subprofit->modnaxisn[0]&1))
          {
          *(pixin++) = 0.0;
          npix2--;
          }
        for (i=npix2; i--;)
          *(pixin++) = *(pixin2--);
        }

/*---- Compute the sharp part of the profile */
      ix1max = subprofit->modnaxisn[0];
      ix2max = subprofit->modnaxisn[1];
      dx1cout = x1cout + 0.4999999;
      dx2cout = x2cout + 0.4999999;
      invdet = 1.0/fabsf(cd11*cd22 - cd12*cd21);
      a11 = cd22*invdet;
      a12 = -cd12*invdet;
      a21 = -cd21*invdet;
      a22 = cd11*invdet;
      nang = 72 / 2;		/* 36 angles; only half of them are computed*/
      angstep = PI/nang;
      ang = 0.0;
      for (a=0; a<nang; a++)
        {
#ifdef HAVE_SINCOSF
        sincosf(ang, &dsa, &dca);
#else
        dsa = sinf(ang);
        dca = cosf(ang);
#endif
        ddx1[a] = a11*dsa+a12*dca;
        ddx2[a] = a21*dsa+a22*dca;
        ang += angstep;
        }
      r = DEXPF(-4.0);
      dr = DEXPF(0.05);
      selem = 0.5*angstep*(dr - 1.0/dr)/(xscale*yscale);
      krpinvn = k*DEXPF(-4.0*invn);
      dkrpinvn = DEXPF(0.05*invn);
      pixin = prof->pix;
      for (; r<rs; r *= dr)
        {
        r2 = r*r;
        val = (expf(krpinvn) - (a0 + r2*(a2+a3*r)))*r2*selem;
        for (a=0; a<nang; a++)
          {
          ix1 = (int)(dx1cout + r*ddx1[a]);
          ix2 = (int)(dx2cout + r*ddx2[a]);
          if (ix1>=0 && ix1<ix1max && ix2>=0 && ix2<ix2max)
            pixin[ix2*ix1max+ix1] += val;
          ix1 = (int)(dx1cout - r*ddx1[a]);
          ix2 = (int)(dx2cout - r*ddx2[a]);
          if (ix1>=0 && ix1<ix1max && ix2>=0 && ix2<ix2max)
            pixin[ix2*ix1max+ix1] += val;
          }
        krpinvn *= dkrpinvn;
        }

      prof->lostfluxfrac = prof->code==MODEL_EXPONENTIAL?
		(1.0 + rmax)*exp(-rmax)
		:1.0 - prof_gammainc(2.0*n, bn*pow(r2max, hinvn));
      threshflag = 0;
      break;
    case MODEL_ARMS:
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
      pixin = prof->pix;
      for (ix2=subprofit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=subprofit->modnaxisn[0]; ix1--;)
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
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
      break;
    case MODEL_BAR:
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
      pixin = prof->pix;
      for (ix2=subprofit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=subprofit->modnaxisn[0]; ix1--;)
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
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
      break;
    case MODEL_INRING:
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
      pixin = prof->pix;
      for (ix2=subprofit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=subprofit->modnaxisn[0]; ix1--;)
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
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
      break;
    case MODEL_OUTRING:
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
      pixin = prof->pix;
      for (ix2=subprofit->modnaxisn[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=subprofit->modnaxisn[0]; ix1--;)
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
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
      break;
    case MODEL_MOFFAT:
      mofbeta = *prof->extra[0];
      mofp = *prof->extra[1];
      off = *prof->extra[2];
      invmofp2 = 1.0 / mofp;
      x10 = -x1cout - dx1;
      x2 = -x2cout - dx2;
      pixin = prof->pix;
      for (ix2=nx2; ix2--; x2+=1.0)
        {
        x1 = x10;
        for (ix1=subprofit->modnaxisn[0]; ix1--; x1+=1.0)
          {
          x1in = fabs(cd12*x2 + cd11*x1);
          x2in = fabs(cd22*x2 + cd21*x1);
          ra = powf(powf(x1in,mofp) + powf(x2in,mofp),invmofp2) - off;
          if (ra<0.0)
            ra = 0.0;
          ra *= ra;
          if (ra>r2max)
            {
            *(pixin++) = 0.0;
            continue;
            }
          val = expf(-logf(1.0 + ra)*mofbeta);
          *(pixin++) = val;
          }
        }
/*---- Copy the symmetric part */
      if ((npix2=(subprofit->modnaxisn[1]-nx2)*subprofit->modnaxisn[0]) > 0)
        {
        pixin2 = pixin - subprofit->modnaxisn[0] - 1;
        if (!(subprofit->modnaxisn[0]&1))
          {
          *(pixin++) = 0.0;
          npix2--;
          }
        for (i=npix2; i--;)
          *(pixin++) = *(pixin2--);
        }
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
      break;
    default:
/*---- Tabulated profile: remap each pixel */
/*---- Initialize multi-dimensional counters */
     for (d=0; d<2; d++)
        {
        posout[d] = 1.0;	/* FITS convention */
        dnaxisn[d] = subprofit->modnaxisn[d] + 0.99999;
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
            posin[d] += (float)prof->naxisn[d];
          else
            posin[d] = 1.0;
          }
        else if (posin[d] > (float)prof->naxisn[d])
          {
          if (prof->extracycleflag[e])
          posin[d] = (prof->extracycleflag[e])?
		  fmod(posin[d], (float)prof->naxisn[d])
		: (float)prof->naxisn[d];
          }
        }
      x1cin = (float)(prof->naxisn[0]/2);
      x2cin = (float)(prof->naxisn[1]/2);
      pixin = prof->pix;
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
      prof->lostfluxfrac = 0.0;
      threshflag = 1;
    break;
    }

/* For complex profiles, threshold to the brightest pixel value at border */
  if (threshflag)
    {
/*-- Now find truncation threshold */
/*-- Find the shortest distance to a vignet border */
    rmax = x1cout;
    if (rmax > (r = x2cout))
      rmax = r;
    rmax += 0.01;
    if (rmax<1.0)
      rmax = 1.0;
    r2max = rmax*rmax;
    rmin = rmax - 1.0;
    r2min = rmin*rmin;

/*-- Find best threshold (the max around the circle with radius rmax */
    dx2 = -x2cout;
    pixin = prof->pix;
    thresh = -BIG;
    for (ix2=subprofit->modnaxisn[1]; ix2--; dx2 += 1.0)
      {
      dx1 = -x1cout;
      for (ix1=subprofit->modnaxisn[0]; ix1--; dx1 += 1.0)
        if ((val=*(pixin++))>thresh && (r2=dx1*dx1+dx2*dx2)>r2min && r2<r2max)
          thresh = val;
      }

/*-- Threshold and measure the flux */
    flux = 0.0;
    pixin = prof->pix;
    for (i=npix; i--; pixin++)
      if (*pixin >= thresh)
        flux += (double)*pixin;
      else
        *pixin = 0.0;
    }
  else
    {
    flux = 0.0;
    pixin = prof->pix;
    for (i=npix; i--;)
      flux += (double)*(pixin++);
    }

  if (extfluxfac_flag)
    fluxfac = prof->fluxfac;
  else
    {
    if (prof->lostfluxfrac < 1.0)
      flux /= (1.0 - prof->lostfluxfrac);

    prof->fluxfac = fluxfac = fabs(flux)>0.0? subprofit->fluxfac/fabs(flux):0.0;
    }

  pixin = prof->pix;
  for (i=npix; i--;)
    *(pixin++) *= fluxfac;

/* Correct final flux */
  pixin = prof->pix;
  pixout = subprofit->modpix;
  for (i=npix; i--;)
    *(pixout++) += pflux**(pixin++);

  QFREE(prof->pix);
  return pflux;
  }


/****** prof_moments **********************************************************
PROTO	int	prof_moments(profitstruct *profit, profstruct *prof)
PURPOSE	Computes (analytically or numerically) the 2nd moments of a profile.
INPUT	Profile-fitting structure,
	profile structure,
	optional pointer to 3xnparam Jacobian matrix.
OUTPUT	Index to the profile flux for further processing.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	20/08/2010
 ***/
int	prof_moments(profitstruct *profit, profstruct *prof, double *jac)
  {
   double	*dmx2,*dmy2,*dmxy,
		m20, a2, ct,st, mx2fac, my2fac,mxyfac, dc2,ds2,dcs,
		bn,bn2, n,n2, nfac,nfac2, hscale2, dmdn;
   int		nparam, index;

  if (jac)
/*-- Clear output Jacobian */
    {
    nparam = profit->nparam;
    memset(jac, 0, nparam*3*sizeof(double));
    dmx2 = jac;
    dmy2 = jac + nparam;
    dmxy = jac + 2*nparam;
    }
  else
    dmx2 = dmy2 = dmxy = NULL;		/* To avoid gcc -Wall warnings */

  m20 = 0.0;				/* to avoid gcc -Wall warnings */
  index = 0;				/* to avoid gcc -Wall warnings */


  if (prof->posangle)
    {
    a2 = *prof->aspect**prof->aspect;
    ct = cos(*prof->posangle*DEG);
    st = sin(*prof->posangle*DEG);
    mx2fac = ct*ct + st*st*a2;
    my2fac = st*st + ct*ct*a2;
    mxyfac = ct*st * (1.0 - a2);
    if (jac)
      {
      dc2 = -2.0*ct*st*DEG;
      ds2 =  2.0*ct*st*DEG;
      dcs = (ct*ct - st*st)*DEG;
      }
    else
      dc2 = ds2 = dcs = 0.0;		/* To avoid gcc -Wall warnings */
    switch(prof->code)
      {
      case MODEL_SERSIC:
        n = fabs(*prof->extra[0]);
        bn = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
		+ 131.0/(1148175*n*n*n);	/* Ciotti & Bertin 1999 */
        nfac  = prof_gamma(4.0*n) / (prof_gamma(2.0*n)*pow(bn, 2.0*n));
        hscale2 = 0.5 * *prof->scale**prof->scale;
        m20 = hscale2 * nfac;
        if (jac)
          {
          dmx2[profit->paramindex[PARAM_SPHEROID_REFF]]
			= *prof->scale * nfac * mx2fac;
          dmy2[profit->paramindex[PARAM_SPHEROID_REFF]]
			= *prof->scale * nfac * my2fac;
          dmxy[profit->paramindex[PARAM_SPHEROID_REFF]]
			= *prof->scale * nfac * mxyfac;
          n2 = n+0.01;
          bn2 = 2.0*n2 - 1.0/3.0 + 4.0/(405.0*n2) + 46.0/(25515.0*n2*n2)
		+ 131.0/(1148175*n2*n2*n2);	/* Ciotti & Bertin 1999 */
          nfac2 = prof_gamma(4.0*n2) / (prof_gamma(2.0*n2)*pow(bn2, 2.0*n2));
          dmdn = 100.0 * hscale2 * (nfac2-nfac);
          dmx2[profit->paramindex[PARAM_SPHEROID_SERSICN]] = dmdn * mx2fac;
          dmy2[profit->paramindex[PARAM_SPHEROID_SERSICN]] = dmdn * my2fac;
          dmxy[profit->paramindex[PARAM_SPHEROID_SERSICN]] = dmdn * mxyfac;
          dmx2[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= 2.0 * m20 * st*st * *prof->aspect;
          dmy2[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= 2.0 * m20 * ct*ct * *prof->aspect;
          dmxy[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= -2.0 * m20 * ct*st * *prof->aspect;
          dmx2[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (dc2+ds2*a2);
          dmy2[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (ds2+dc2*a2);
          dmxy[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (1.0-a2)*dcs;
          }
        index = profit->paramindex[PARAM_SPHEROID_FLUX];
        break; 
      case MODEL_DEVAUCOULEURS:
        m20 = 10.83995 * *prof->scale**prof->scale;
        if (jac)
          {
          dmx2[profit->paramindex[PARAM_SPHEROID_REFF]]
			= 21.680 * *prof->scale * mx2fac;
          dmy2[profit->paramindex[PARAM_SPHEROID_REFF]]
			= 21.680 * *prof->scale * my2fac;
          dmxy[profit->paramindex[PARAM_SPHEROID_REFF]]
			= 21.680 * *prof->scale * mxyfac;
          dmx2[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= 2.0 * m20 * st*st * *prof->aspect;
          dmy2[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= 2.0 * m20 * ct*ct * *prof->aspect;
          dmxy[profit->paramindex[PARAM_SPHEROID_ASPECT]]
			= -2.0 * m20 * ct*st * *prof->aspect;
          dmx2[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (dc2+ds2*a2);
          dmy2[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (ds2+dc2*a2);
          dmxy[profit->paramindex[PARAM_SPHEROID_POSANG]] = m20 * (1.0-a2)*dcs;
          }
        index = profit->paramindex[PARAM_SPHEROID_FLUX];
        break;
      case MODEL_EXPONENTIAL:
        m20 = 3.0 * *prof->scale**prof->scale;
        if (jac)
          {
          dmx2[profit->paramindex[PARAM_DISK_SCALE]]
			= 6.0 * *prof->scale * mx2fac;
          dmy2[profit->paramindex[PARAM_DISK_SCALE]]
			= 6.0 * *prof->scale * my2fac;
          dmxy[profit->paramindex[PARAM_DISK_SCALE]]
			= 6.0 * *prof->scale * mxyfac;
          dmx2[profit->paramindex[PARAM_DISK_ASPECT]]
			= 2.0 * m20 * st*st * *prof->aspect;
          dmy2[profit->paramindex[PARAM_DISK_ASPECT]]
			= 2.0 * m20 * ct*ct * *prof->aspect;
          dmxy[profit->paramindex[PARAM_DISK_ASPECT]]
			= -2.0 * m20 * ct*st * *prof->aspect;
          dmx2[profit->paramindex[PARAM_DISK_POSANG]] = m20 * (dc2 + ds2*a2);
          dmy2[profit->paramindex[PARAM_DISK_POSANG]] = m20 * (ds2 + dc2*a2);
          dmxy[profit->paramindex[PARAM_DISK_POSANG]] = m20 * (1.0 - a2) * dcs;
          }
        index = profit->paramindex[PARAM_DISK_FLUX];
        break;
      case MODEL_ARMS:
        m20 = 1.0;
        index = profit->paramindex[PARAM_ARMS_FLUX];
        break;
      case MODEL_BAR:
        m20 = 1.0;
        index = profit->paramindex[PARAM_BAR_FLUX];
        break;
      case MODEL_INRING:
        m20 = 1.0;
        index = profit->paramindex[PARAM_INRING_FLUX];
        break;
      case MODEL_OUTRING:
        m20 = 1.0;
        index = profit->paramindex[PARAM_OUTRING_FLUX];
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown oriented model in ",
		"prof_moments()");
      break;
      }

    prof->mx2 = m20*mx2fac;
    prof->my2 = m20*my2fac;
    prof->mxy = m20*mxyfac;
    
    }
  else
    prof->mx2 = prof->my2 = prof->mxy = 0.0;

  return index;
  }


/****** prof_interpolate ******************************************************
PROTO	float	prof_interpolate(profstruct *prof, float *posin)
PURPOSE	Interpolate a multidimensional model profile at a given position.
INPUT	Profile structure,
	input position vector.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/12/2006
 ***/
static float	prof_interpolate(profstruct *prof, float *posin)
  {
   float		dpos[2+PROFIT_MAXEXTRA],
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
PROTO	void interpolate_pix(float *posin, float *pix, int naxisn,
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
static float	interpolate_pix(float *posin, float *pix, int *naxisn,
			interpenum interptype)
  {
   float	buffer[INTERP_MAXKERNELWIDTH],
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
PROTO	void make_kernel(float pos, float *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/07/2011
 ***/
void	make_kernel(float pos, float *kernel, interpenum interptype)
  {
   float	x, val, sinx1,sinx2,sinx3,cosx1;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR)
    {
    *(kernel++) = 1.0-pos;
    *kernel = pos;
    }
  else if (interptype == INTERP_LANCZOS2)
    {
    if (pos<1e-5 && pos>-1e-5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/2.0*(pos+1.0);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
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
    if (pos<1e-5 && pos>-1e-5)
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
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
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
    if (pos<1e-5 && pos>-1e-5)
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
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
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

