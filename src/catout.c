/*
*				catout.c
*
* Functions related to catalogue output.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2023 CFHT/IAP/CNRS/SorbonneU
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
*	Last modified:		26/02/2023
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
#include	"param.h"
#include	"sexhead.h"
#include	"sexhead1.h"
#include	"sexheadsc.h"
#include	"xml.h"

objstruct	outobj, flagobj;
obj2struct	outobj2, flagobj2;

catstruct	*fitscat;
tabstruct	*objtab = NULL;
FILE		*ascfile;
char		*buf;
int		catopen_flag = 0;

double		ddummy;
int		idummy;

/******************************* readcatparams *******************************/
/*
Read the catalog config file
*/
void	readcatparams(char *filename)
  {
   keystruct	*key;
   FILE		*infile;
   char		str[MAXCHAR], *keyword, *sstr;
   int		i, size;

/* Prepare the OBJECTS tables*/
  objtab = new_tab("OBJECTS");

  if ((infile = fopen(filename,"r")) == NULL)
    error(EXIT_FAILURE, "*ERROR*: can't read ", filename);

/* Scan the catalog config file*/
  thecat.nparam = 0;
  while (fgets(str, MAXCHAR, infile))
    {
    sstr = str + strspn(str," \t");
    if (*sstr!=(char)'#' && *sstr!=(char)'\n')
      {
      keyword = strtok(sstr, " \t{[(\n\r");
      if (keyword &&
	(i = findkey(keyword,(char *)objkey,sizeof(keystruct)))!=RETURN_ERROR)
        {
        key = objkey+i;
        if (add_key(key, objtab, 0) != RETURN_ERROR)
          thecat.nparam++;
        *((char *)key->ptr) = (char)'\1';
        if (key->naxis)
          {
          for (i=0; i<key->naxis; i++)
            key->naxisn[i] = 1;
          size=t_size[key->ttype];
          for (i=0; (sstr = strtok(NULL, " \t,;.)]}\r")) && *sstr!=(char)'#'
		&& *sstr!=(char)'\n'; i++)
            {
            if (i>=key->naxis)
              error(EXIT_FAILURE, "*Error*: too many dimensions for keyword ",
		keyword);
            if (!(size*=(key->naxisn[i]=atoi(sstr))))
              error(EXIT_FAILURE, "*Error*: wrong array syntax for keyword ",
		keyword);
            }
          key->nbytes = size;
          }
        }
      else
        warning(keyword, " catalog parameter unknown");
      }
    }

  fclose(infile);

/* Now we copy the flags to the proper structures */

  flagobj = outobj;
  flagobj2 = outobj2;
/* Differentiate between outobj and outobj2 vectors */
  memset(&outobj2, 0, sizeof(outobj2));
  updateparamflags();

  return;
  }


/******************************* alloccatparams ******************************/
/*
Allocate memory for parameter arrays
*/
void	alloccatparams(void)
  {
   keystruct	*key;
   int		i;

/* Go back to multi-dimensional arrays for memory allocation */
  if (thecat.nparam)
    for (i=objtab->nkey, key=objtab->key; i--; key = key->nextkey) 
      if (key->naxis)
        {
/*------ Only outobj2 vectors are dynamic */
        if (!*((char **)key->ptr))
          {
          QMALLOC(*((char **)key->ptr), char, key->nbytes);
          key->ptr = *((char **)key->ptr);
          key->allocflag = 1;
          }
        }

  return;
  }

/*************************** changecatparamarrays ****************************/
/*
Change parameter array dimensions
*/
void	changecatparamarrays(char *keyword, int *axisn, int naxis)
  {
   keystruct	*key;
   int		d,i, size;

  if (thecat.nparam)
    for (i=objtab->nkey, key=objtab->key; i--; key = key->nextkey) 
      if (key->naxis && !strcmp(keyword, key->name))
        {
        size = t_size[key->ttype];
        if (key->naxis != naxis)
          key->naxis = naxis;
        for (d=0; d<naxis; d++)
          size *= (key->naxisn[d]=axisn[d]);
        key->nbytes = size;
        break;
        }

  return;
  }

/***************************** updateparamflags ******************************/
/*
Update parameter flags according to their mutual dependencies.
*/
void	updateparamflags()

  {
   int	i;

/*----------------------------- Model-fitting -----------------------------*/

  prefs.pattern_flag |= FLAG(obj2.prof_disk_patternvector)
			| FLAG(obj2.prof_disk_patternmodvector)
			| FLAG(obj2.prof_disk_patternargvector)
			| FLAG(obj2.prof_disk_patternspiral);
  FLAG(obj2.prof_disk_scale) |= prefs.pattern_flag;

  FLAG(obj2.prof_concentration) |= FLAG(obj2.prof_concentrationerr);
  FLAG(obj2.fluxcorerr_prof) |= FLAG(obj2.magcorerr_prof);
  FLAG(obj2.fluxcor_prof) |= FLAG(obj2.magcor_prof)
			| FLAG(obj2.fluxcorerr_prof);

  FLAG(obj2.poserraw_prof) |= FLAG(obj2.poserrbw_prof);
  FLAG(obj2.poserrcxxw_prof) |= FLAG(obj2.poserrcyyw_prof)
			| FLAG(obj2.poserrcxyw_prof);
  FLAG(obj2.poserrthetas_prof) |= FLAG(obj2.poserrtheta1950_prof)
				| FLAG(obj2.poserrtheta2000_prof);
  FLAG(obj2.poserrthetaw_prof) |= FLAG(obj2.poserrthetas_prof);

  FLAG(obj2.poserrmx2w_prof) |= FLAG(obj2.poserrmy2w_prof)
			| FLAG(obj2.poserrmxyw_prof)
			| FLAG(obj2.poserrthetaw_prof)|FLAG(obj2.poserraw_prof)
			| FLAG(obj2.poserrcxxw_prof);

  FLAG(obj2.poserra_prof) |= FLAG(obj2.poserrb_prof)
			| FLAG(obj2.poserrtheta_prof)
			| FLAG(obj2.poserraw_prof);
  FLAG(obj2.poserrcxx_prof) |= FLAG(obj2.poserrcyy_prof)
			| FLAG(obj2.poserrcxy_prof);
  FLAG(obj2.poserrmx2_prof) |= FLAG(obj2.poserrmy2_prof)
			| FLAG(obj2.poserrmxy_prof)
			| FLAG(obj2.poserrmx2w_prof);
  FLAG(obj2.alpha1950_prof) |= FLAG(obj2.delta1950_prof)
			| FLAG(obj2.poserrtheta1950_prof);
  FLAG(obj2.alpha2000_prof) |= FLAG(obj2.delta2000_prof)
			| FLAG(obj2.alpha1950_prof)
			| FLAG(obj2.poserrtheta2000_prof);
  FLAG(obj2.alphas_prof) |= FLAG(obj2.deltas_prof)
			| FLAG(obj2.alpha2000_prof);
  FLAG(obj2.xw_prof) |= FLAG(obj2.yw_prof)
			| FLAG(obj2.alphas_prof);
  FLAG(obj2.xf_prof) |= FLAG(obj2.yf_prof);
  FLAG(obj2.x_prof) |= FLAG(obj2.y_prof)
			| FLAG(obj2.xw_prof)
			| FLAG(obj2.xf_prof)
			| FLAG(obj2.poserra_prof)
			| FLAG(obj2.poserrcxx_prof)
			| FLAG(obj2.poserrmx2_prof)
			| FLAG(obj2.prof_concentration)
			| FLAG(obj2.prof_class_star)
			| FLAG(obj2.fluxcor_prof);

  FLAG(obj2.prof_dirac_fluxratio) |= FLAG(obj2.prof_dirac_fluxratioerr);
  FLAG(obj2.prof_dirac_mag) |= FLAG(obj2.prof_dirac_magerr);
  FLAG(obj2.prof_spheroid_fluxratio) |= FLAG(obj2.prof_spheroid_fluxratioerr);
  FLAG(obj2.prof_spheroid_mag) |= FLAG(obj2.prof_spheroid_magerr);
  FLAG(obj2.prof_spheroid_reff) |= FLAG(obj2.prof_spheroid_refferr);
  FLAG(obj2.prof_spheroid_aspect) |= FLAG(obj2.prof_spheroid_aspecterr);
  FLAG(obj2.prof_spheroid_theta) |= FLAG(obj2.prof_spheroid_thetaerr);
  FLAG(obj2.prof_spheroid_sersicn) |= FLAG(obj2.prof_spheroid_sersicnerr);
  FLAG(obj2.prof_disk_fluxratio) |= FLAG(obj2.prof_disk_fluxratioerr);
  FLAG(obj2.prof_disk_mag) |= FLAG(obj2.prof_disk_magerr);
  FLAG(obj2.prof_disk_scale) |= FLAG(obj2.prof_disk_scaleerr);
  FLAG(obj2.prof_disk_aspect) |= FLAG(obj2.prof_disk_aspecterr);
  FLAG(obj2.prof_disk_inclination) |= FLAG(obj2.prof_disk_inclinationerr);
  FLAG(obj2.prof_disk_theta) |= FLAG(obj2.prof_disk_thetaerr);
  FLAG(obj2.prof_bar_fluxratio) |= FLAG(obj2.prof_bar_fluxratioerr);
  FLAG(obj2.prof_bar_mag) |= FLAG(obj2.prof_bar_magerr);
  FLAG(obj2.prof_bar_length) |= FLAG(obj2.prof_bar_lengtherr);
  FLAG(obj2.prof_bar_aspect) |= FLAG(obj2.prof_bar_aspecterr);
  FLAG(obj2.prof_bar_theta) |= FLAG(obj2.prof_bar_thetaerr);
  FLAG(obj2.prof_arms_fluxratio) |= FLAG(obj2.prof_arms_fluxratioerr);
  FLAG(obj2.prof_arms_mag) |= FLAG(obj2.prof_arms_magerr);
  FLAG(obj2.prof_e1errw) |= FLAG(obj2.prof_e2errw) | FLAG(obj2.prof_e12corrw);
  FLAG(obj2.prof_e1w) |= FLAG(obj2.prof_e2w) | FLAG(obj2.prof_e1errw);
  FLAG(obj2.prof_pol1errw) |= FLAG(obj2.prof_pol2errw)
			| FLAG(obj2.prof_pol12corrw);
  FLAG(obj2.prof_pol1w) |= FLAG(obj2.prof_pol2w) | FLAG(obj2.prof_pol1errw);
  FLAG(obj2.prof_aw) |= FLAG(obj2.prof_bw);
  FLAG(obj2.prof_cxxw) |= FLAG(obj2.prof_cyyw) | FLAG(obj2.prof_cxyw);
  FLAG(obj2.prof_thetas) |= FLAG(obj2.prof_theta1950)
			| FLAG(obj2.prof_theta2000);
  FLAG(obj2.prof_thetaw) |= FLAG(obj2.prof_thetas);

  FLAG(obj2.prof_mx2w) |= FLAG(obj2.prof_my2w)
			| FLAG(obj2.prof_mxyw)
			| FLAG(obj2.prof_thetaw) | FLAG(obj2.prof_aw)
			| FLAG(obj2.prof_cxxw)
			| FLAG(obj2.prof_e1w) | FLAG(obj2.prof_pol1w);

  FLAG(obj2.dtheta1950) |= FLAG(obj2.prof_theta1950)
			| FLAG(obj2.prof_spheroid_theta1950)
			| FLAG(obj2.prof_disk_theta1950)
//			| FLAG(obj2.prof_arms_theta1950)
			| FLAG(obj2.prof_bar_theta1950)
			| FLAG(obj2.poserrtheta1950_psf)
			| FLAG(obj2.win_theta1950)
			| FLAG(obj2.winposerr_theta1950)
			| FLAG(obj2.poserr_theta1950)
			| FLAG(obj2.theta1950)
			| FLAG(obj2.poserrtheta1950_prof);
  FLAG(obj2.dtheta2000) |= FLAG(obj2.prof_theta2000)
			| FLAG(obj2.prof_spheroid_theta2000)
			| FLAG(obj2.prof_disk_theta2000)
//			| FLAG(obj2.prof_arms_theta2000)
			| FLAG(obj2.prof_bar_theta2000)
			| FLAG(obj2.poserrtheta2000_psf)
			| FLAG(obj2.win_theta2000)
			| FLAG(obj2.winposerr_theta2000)
			| FLAG(obj2.poserr_theta2000)
			| FLAG(obj2.theta2000)
			| FLAG(obj2.poserrtheta2000_prof);


  FLAG(obj2.prof_spheroid_thetas) |= FLAG(obj2.prof_spheroid_theta2000)
			| FLAG(obj2.prof_spheroid_theta1950);
  FLAG(obj2.prof_spheroid_thetaw) |= FLAG(obj2.prof_spheroid_thetas);
  FLAG(obj2.prof_disk_thetas) |= FLAG(obj2.prof_disk_theta2000)
			| FLAG(obj2.prof_disk_theta1950);
  FLAG(obj2.prof_disk_thetaw) |= FLAG(obj2.prof_disk_thetas);
  FLAG(obj2.prof_bar_thetas) |= FLAG(obj2.prof_bar_theta2000)
			| FLAG(obj2.prof_bar_theta1950);
  FLAG(obj2.prof_bar_thetaw) |= FLAG(obj2.prof_bar_thetas);
/*
  FLAG(obj2.prof_arms_thetaw) |= FLAG(obj2.prof_arms_thetas)
			| FLAG(obj2.prof_arms_theta2000)
			| FLAG(obj2.prof_arms_theta1950);
*/
  FLAG(obj2.prof_arms_scalew) |= FLAG(obj2.prof_arms_startw)
			| FLAG(obj2.prof_arms_scaleerrw)
			| FLAG(obj2.prof_arms_starterrw);

  FLAG(obj2.prof_bar_lengthw) |= FLAG(obj2.prof_bar_aspectw)
			| FLAG(obj2.prof_bar_thetaw)
			| FLAG(obj2.prof_bar_lengtherrw)
			| FLAG(obj2.prof_bar_aspecterrw)
			| FLAG(obj2.prof_bar_thetaerrw);

  FLAG(obj2.prof_disk_scalew) |= FLAG(obj2.prof_disk_aspectw)
			| FLAG(obj2.prof_disk_thetaw)
			| FLAG(obj2.prof_disk_scaleerrw)
			| FLAG(obj2.prof_disk_aspecterrw)
			| FLAG(obj2.prof_disk_thetaerrw)
			| FLAG(obj2.prof_arms_scalew);
  FLAG(obj2.prof_disk_fluxmean) |= FLAG(obj2.prof_disk_mumean);
  FLAG(obj2.prof_disk_fluxeff) |= FLAG(obj2.prof_disk_mueff);
  FLAG(obj2.prof_disk_peak) |= FLAG(obj2.prof_disk_mumax)
				| FLAG(obj2.prof_disk_fluxeff)
				| FLAG(obj2.prof_disk_fluxmean);

  FLAG(obj2.prof_spheroid_reffw) |= FLAG(obj2.prof_spheroid_aspectw)
			| FLAG(obj2.prof_spheroid_thetaw)
			| FLAG(obj2.prof_spheroid_refferrw)
			| FLAG(obj2.prof_spheroid_aspecterrw)
			| FLAG(obj2.prof_spheroid_thetaerrw);

  FLAG(obj2.prof_spheroid_fluxmean) |= FLAG(obj2.prof_spheroid_mumean);
  FLAG(obj2.prof_spheroid_fluxeff) |= FLAG(obj2.prof_spheroid_mueff);
  FLAG(obj2.prof_spheroid_peak) |= FLAG(obj2.prof_spheroid_mumax)
				| FLAG(obj2.prof_spheroid_fluxeff)
				| FLAG(obj2.prof_spheroid_fluxmean);

  FLAG(obj2.prof_flagw) |= FLAG(obj2.prof_spheroid_reffw)
			| FLAG(obj2.prof_disk_scalew)
			| FLAG(obj2.prof_bar_lengthw)
			| FLAG(obj2.prof_arms_scalew)
			| FLAG(obj2.prof_mx2w);

  FLAG(obj2.prof_e1err) |= FLAG(obj2.prof_e2err) | FLAG(obj2.prof_e12corr)
			| FLAG(obj2.prof_e1errw);
			
  FLAG(obj2.prof_e1) |= FLAG(obj2.prof_e2)
			| FLAG(obj2.prof_e1err) | FLAG(obj2.prof_e1w);
  FLAG(obj2.prof_pol1err) |= FLAG(obj2.prof_pol2err)
			| FLAG(obj2.prof_pol12corr) | FLAG(obj2.prof_pol1errw);
  FLAG(obj2.prof_pol1) |= FLAG(obj2.prof_pol2)
			| FLAG(obj2.prof_pol1err) | FLAG(obj2.prof_pol1w);
  FLAG(obj2.prof_a) |= FLAG(obj2.prof_b) | FLAG(obj2.prof_theta)
			| FLAG(obj2.prof_aw);
  FLAG(obj2.prof_cxx) |= FLAG(obj2.prof_cyy)
			| FLAG(obj2.prof_cxy) | FLAG(obj2.prof_cxxw);
  FLAG(obj2.prof_mx2) |= FLAG(obj2.prof_my2) | FLAG(obj2.prof_mxy)
			| FLAG(obj2.prof_e1) | FLAG(obj2.prof_pol1)
			| FLAG(obj2.prof_a) | FLAG(obj2.prof_cxx)
			| FLAG(obj2.prof_mx2w);
  FLAG(obj2.prof_conva) |= FLAG(obj2.prof_convb) | FLAG(obj2.prof_convtheta);
  FLAG(obj2.prof_convcxx) |= FLAG(obj2.prof_convcyy)
			/*| FLAG(obj2.fluxcor_prof) */;
  FLAG(obj2.prof_convmx2) |= FLAG(obj2.prof_convcxx) | FLAG(obj2.prof_conva);
  FLAG(obj2.fluxmean_prof) |= FLAG(obj2.mumean_prof);
  FLAG(obj2.peak_prof) |= FLAG(obj2.mumax_prof);
  FLAG(obj2.fluxeff_prof) |= FLAG(obj2.mueff_prof)
			| FLAG(obj2.fluxmean_prof) | FLAG(obj2.peak_prof);

  FLAG(obj2.prof_arms_flux) |= FLAG(obj2.prof_arms_fluxerr)
			| FLAG(obj2.prof_arms_mag)
			| FLAG(obj2.prof_arms_fluxratio)
			| FLAG(obj2.prof_arms_scalew)
			| FLAG(obj2.prof_arms_scale)
			| FLAG(obj2.prof_arms_posang)
			| FLAG(obj2.prof_arms_pitch)
			| FLAG(obj2.prof_arms_start)
			| FLAG(obj2.prof_arms_quadfrac);
  FLAG(obj2.prof_bar_theta) |= FLAG(obj2.prof_bar_lengthw);
  FLAG(obj2.prof_bar_flux) |= FLAG(obj2.prof_bar_fluxerr)
			| FLAG(obj2.prof_bar_mag)
			| FLAG(obj2.prof_bar_fluxratio)
			| FLAG(obj2.prof_bar_lengthw)
			| FLAG(obj2.prof_bar_length)
			| FLAG(obj2.prof_bar_aspect)
			| FLAG(obj2.prof_bar_posang)
			| FLAG(obj2.prof_bar_theta)
			| FLAG(obj2.prof_arms_flux);
  FLAG(obj2.prof_disk_flux) |= FLAG(obj2.prof_disk_fluxerr)
			| FLAG(obj2.prof_disk_mag)
			| FLAG(obj2.prof_disk_fluxratio)
			| FLAG(obj2.prof_disk_scalew)
			| FLAG(obj2.prof_disk_scale)
			| FLAG(obj2.prof_disk_aspect)
			| FLAG(obj2.prof_disk_inclination)
			| FLAG(obj2.prof_disk_theta)
			| FLAG(obj2.prof_disk_peak)
			| FLAG(obj2.prof_bar_flux);
  FLAG(obj2.prof_spheroid_flux) |= FLAG(obj2.prof_spheroid_fluxerr)
			| FLAG(obj2.prof_spheroid_mag)
			| FLAG(obj2.prof_spheroid_fluxratio)
			| FLAG(obj2.prof_spheroid_reffw)
			| FLAG(obj2.prof_spheroid_reff)
			| FLAG(obj2.prof_spheroid_aspect)
			| FLAG(obj2.prof_spheroid_theta)
			| FLAG(obj2.prof_spheroid_sersicn)
			| FLAG(obj2.prof_spheroid_peak);
  FLAG(obj2.prof_dirac_flux) |= FLAG(obj2.prof_dirac_fluxerr)
			| FLAG(obj2.prof_dirac_mag)
			| FLAG(obj2.prof_dirac_fluxratio);
  FLAG(obj2.prof_offset_flux) |= FLAG(obj2.prof_offset_fluxerr);
  FLAG(obj2.fluxerr_dprof) |= FLAG(obj2.magerr_dprof);
  FLAG(obj2.flux_dprof) |= FLAG(obj2.mag_dprof)
			| FLAG(obj2.fluxerr_dprof);
  prefs.dprof_flag |= FLAG(obj2.dprof_chi2)
			| FLAG(obj2.dprof_flag)
			| FLAG(obj2.dprof_niter)
			| FLAG(obj2.flux_dprof);

  FLAG(obj2.fluxerr_prof) |= FLAG(obj2.magerr_prof)
			| FLAG(obj2.prof_concentrationerr)
			| FLAG(obj2.fluxcorerr_prof)
			| FLAG(obj2.prof_arms_fluxratio)
			| FLAG(obj2.prof_bar_fluxratio)
			| FLAG(obj2.prof_disk_fluxratio)
			| FLAG(obj2.prof_spheroid_fluxratio)
			| FLAG(obj2.prof_dirac_fluxratio);
  FLAG(obj2.flux_prof) |= FLAG(obj2.mag_prof)
			| FLAG(obj2.fluxerr_prof)
			| FLAG(obj2.fluxcor_prof);

  prefs.prof_flag |= FLAG(obj2.prof_chi2) | FLAG(obj2.prof_niter)
			| FLAG(obj2.prof_vector) | FLAG(obj2.prof_errvector)
			| FLAG(obj2.prof_errmatrix)
			| FLAG(obj2.x_prof) | FLAG(obj2.y_prof)
			| FLAG(obj2.prof_flag)
			| FLAG(obj2.flux_prof)
			| FLAG(obj2.prof_mx2)
			| FLAG(obj2.fluxeff_prof)
			| FLAG(obj2.prof_disk_flux)
			| FLAG(obj2.prof_spheroid_flux)
			| FLAG(obj2.prof_dirac_flux)
			| FLAG(obj2.prof_offset_flux)
			| prefs.dprof_flag;


/* If only global parameters are requested, fit a Sersic model */
  if (prefs.prof_flag && !(FLAG(obj2.prof_spheroid_flux)
			| FLAG(obj2.prof_disk_flux)
			| FLAG(obj2.prof_dirac_flux)
			| FLAG(obj2.prof_offset_flux)))
    {
    FLAG(obj2.prof_spheroid_flux) |= prefs.prof_flag;
    FLAG(obj2.prof_spheroid_sersicn) |= prefs.prof_flag;
    }


/*------------------------------ Astrometry ---------------------------------*/
  FLAG(obj2.win_aw) |= FLAG(obj2.win_bw) | FLAG(obj2.win_polarw);
  FLAG(obj2.win_cxxw) |= FLAG(obj2.win_cyyw) | FLAG(obj2.win_cxyw);
  FLAG(obj2.win_thetas) |= FLAG(obj2.win_theta1950)
			| FLAG(obj2.win_theta2000);
  FLAG(obj2.win_thetaw) |= FLAG(obj2.win_thetas);

  FLAG(obj2.win_mx2w) |= FLAG(obj2.win_my2w)
			| FLAG(obj2.win_mxyw)
			| FLAG(obj2.win_thetaw) | FLAG(obj2.win_aw)
			| FLAG(obj2.win_cxxw);

  FLAG(obj2.win_a) |= FLAG(obj2.win_b) | FLAG(obj2.win_theta)
			| FLAG(obj2.win_polar) | FLAG(obj2.win_aw);
  FLAG(obj2.win_cxx) |= FLAG(obj2.win_cyy)
			| FLAG(obj2.win_cxy) | FLAG(obj2.win_cxxw);
  FLAG(obj2.win_mx2) |= FLAG(obj2.win_my2)
			| FLAG(obj2.win_mxy)
			| FLAG(obj2.win_a) | FLAG(obj2.win_cxx)
			| FLAG(obj2.win_mx2w);

  FLAG(obj2.winposerr_aw) |= FLAG(obj2.winposerr_bw);
  FLAG(obj2.winposerr_cxxw) |= FLAG(obj2.winposerr_cyyw)
			| FLAG(obj2.winposerr_cxyw);
  FLAG(obj2.winposerr_thetas) |= FLAG(obj2.winposerr_theta1950)
			| FLAG(obj2.winposerr_theta2000);
  FLAG(obj2.winposerr_thetaw) |= FLAG(obj2.winposerr_thetas);

  FLAG(obj2.winposerr_mx2w) |= FLAG(obj2.winposerr_my2w)
			| FLAG(obj2.winposerr_mxyw)
			| FLAG(obj2.winposerr_thetaw) | FLAG(obj2.winposerr_aw)
			| FLAG(obj2.winposerr_cxxw);

  FLAG(obj2.winposerr_a) |= FLAG(obj2.winposerr_b) | FLAG(obj2.winposerr_theta);
  FLAG(obj2.winposerr_cxx) |= FLAG(obj2.winposerr_cyy)
			| FLAG(obj2.winposerr_cxy);
  FLAG(obj2.winposerr_mx2) |= FLAG(obj2.winposerr_my2)
			| FLAG(obj2.winposerr_mxy)
			| FLAG(obj2.winposerr_a) | FLAG(obj2.winposerr_cxx)
			| FLAG(obj2.winposerr_mx2w);

  FLAG(obj2.winpos_alpha1950) |= FLAG(obj2.winpos_delta1950)
			| FLAG(obj2.win_theta1950)
			| FLAG(obj2.winposerr_theta1950);
  FLAG(obj2.winpos_alpha2000) |= FLAG(obj2.winpos_delta2000)
			| FLAG(obj2.winpos_alpha1950)
			| FLAG(obj2.win_theta2000)
			| FLAG(obj2.winposerr_theta2000);
  FLAG(obj2.winpos_alphas) |= FLAG(obj2.winpos_deltas)
			| FLAG(obj2.winpos_alpha2000);
  FLAG(obj2.winpos_xf) |= FLAG(obj2.winpos_yf);
  FLAG(obj2.winpos_xw) |= FLAG(obj2.winpos_yw)
			| FLAG(obj2.winpos_alphas);

  FLAG(obj2.poserr_aw) |= FLAG(obj2.poserr_bw);
  FLAG(obj2.poserr_cxxw) |= FLAG(obj2.poserr_cyyw) | FLAG(obj2.poserr_cxyw);
  FLAG(obj2.poserr_thetas) |= FLAG(obj2.poserr_theta1950)
				| FLAG(obj2.poserr_theta2000);
  FLAG(obj2.poserr_thetaw) |= FLAG(obj2.poserr_thetas);

  FLAG(obj2.poserr_mx2w) |= FLAG(obj2.poserr_my2w) | FLAG(obj2.poserr_mxyw)
			| FLAG(obj2.poserr_thetaw) | FLAG(obj2.poserr_aw)
			| FLAG(obj2.poserr_cxxw);

  FLAG(obj2.poserr_a) |= FLAG(obj2.poserr_b) | FLAG(obj2.poserr_theta)
			| FLAG(obj2.winposerr_a);
  FLAG(obj2.poserr_cxx) |= FLAG(obj2.poserr_cyy) | FLAG(obj2.poserr_cxy);
  FLAG(obj.poserr_mx2) |= FLAG(obj.poserr_my2) | FLAG(obj.poserr_mxy)
			| FLAG(obj2.poserr_a) | FLAG(obj2.poserr_cxx)
			| FLAG(obj2.poserr_mx2w) | FLAG(obj2.winposerr_mx2);

  FLAG(obj2.peakalpha1950) |= FLAG(obj2.peakdelta1950);
  FLAG(obj2.alpha1950) |= FLAG(obj2.delta1950) |  FLAG(obj2.theta1950)
			| FLAG(obj2.poserr_theta1950) | FLAG(obj2.dtheta1950);
;
  FLAG(obj2.peakalpha2000) |= FLAG(obj2.peakdelta2000)
			| FLAG(obj2.peakalpha1950);
  FLAG(obj2.alpha2000) |= FLAG(obj2.delta2000) | FLAG(obj2.alpha1950)
			| FLAG(obj2.dtheta2000);
  FLAG(obj2.peakalphas) |= FLAG(obj2.peakdeltas) | FLAG(obj2.peakalpha2000);
  FLAG(obj2.alphas) |= FLAG(obj2.deltas) | FLAG(obj2.alpha2000);
  FLAG(obj2.thetas) |= FLAG(obj2.theta1950) | FLAG(obj2.theta2000);
  FLAG(obj2.thetaw) |= FLAG(obj2.thetas);
  FLAG(obj2.aw) |= FLAG(obj2.bw) | FLAG(obj2.polarw);
  FLAG(obj2.cxxw) |= FLAG(obj2.cyyw) | FLAG(obj2.cxyw);

  FLAG(obj2.mx2w) |= FLAG(obj2.my2w) | FLAG(obj2.mxyw)
			| FLAG(obj2.thetaw) | FLAG(obj2.aw) | FLAG(obj2.cxxw);
  
  FLAG(obj2.peakxw) |= FLAG(obj2.peakyf);
  FLAG(obj2.peakxw) |= FLAG(obj2.peakyw) | FLAG(obj2.peakalphas);
  FLAG(obj.peakx) |= FLAG(obj.peaky) | FLAG(obj2.peakxw) | FLAG(obj2.peakxf);

  FLAG(obj2.mxf) |= FLAG(obj2.myf);

  FLAG(obj2.mxw) |= FLAG(obj2.myw) | FLAG(obj2.mx2w) | FLAG(obj2.alphas)
		| FLAG(obj2.poserr_mx2w)
		| ((FLAG(obj2.assoc) | FLAG(obj2.assoc_number))
			&& (prefs.assoccoord_type==ASSOCCOORD_WORLD));
  FLAG(obj2.mamaposx) |= FLAG(obj2.mamaposy);
  FLAG(obj2.fluxerr_win) |= FLAG(obj2.snr_win);
  FLAG(obj2.flux_win) |= FLAG(obj2.mag_win)|FLAG(obj2.magerr_win)
			| FLAG(obj2.flux_win) | FLAG(obj2.fluxerr_win);
  FLAG(obj2.winpos_dgeox) |= FLAG(obj2.winpos_dgeoy);
  FLAG(obj2.winpos_x) |= FLAG(obj2.winpos_y) | FLAG(obj2.winpos_dgeox)
			| FLAG(obj2.winposerr_mx2) | FLAG(obj2.win_mx2)
			| FLAG(obj2.winpos_xw) | FLAG(obj2.winpos_xf)
			| FLAG(obj2.win_flag)
			| FLAG(obj2.flux_win) |FLAG(obj2.winpos_niter);

  FLAG(obj2.area_flagw) |= FLAG(obj2.npixw) | FLAG(obj2.fdnpixw)
			| FLAG(obj2.fwhmw)
			| FLAG(obj2.maxmu) | FLAG(obj2.threshmu)
			| FLAG(obj2.prof_spheroid_mumax)
			| FLAG(obj2.prof_disk_mumax)
			| FLAG(obj2.mumax_prof);

/*------------------------------ Photometry ---------------------------------*/

  FLAG(obj2.fluxerr_best) |= FLAG(obj2.magerr_best);

  FLAG(obj2.flux_best) |= FLAG(obj2.mag_best) | FLAG(obj2.fluxerr_best);

  FLAG(obj2.hl_radius) |= FLAG(obj2.winpos_x) | prefs.prof_flag;

  FLAG(obj2.flux_auto)  |= FLAG(obj2.mag_auto) | FLAG(obj2.magerr_auto)
			| FLAG(obj2.fluxerr_auto)
			| FLAG(obj2.kronfactor)
			| FLAG(obj2.flux_best)
			| FLAG(obj2.flux_radius)
			| FLAG(obj2.hl_radius)
			| FLAG(obj2.fluxcor_prof);
  FLAG(obj2.flux_petro) |= FLAG(obj2.mag_petro) | FLAG(obj2.magerr_petro)
			| FLAG(obj2.fluxerr_petro)
			| FLAG(obj2.petrofactor);

  FLAG(obj2.fluxerr_isocor) |= FLAG(obj2.magerr_isocor)
				| FLAG(obj2.fluxerr_best);

  FLAG(obj2.flux_isocor) |= FLAG(obj2.mag_isocor) | FLAG(obj2.fluxerr_isocor)
			 | FLAG(obj2.flux_best);

  FLAG(obj2.flux_aper) |= FLAG(obj2.mag_aper)|FLAG(obj2.magerr_aper)
			    | FLAG(obj2.fluxerr_aper);

  FLAG(obj2.flux_galfit) |= FLAG(obj2.mag_galfit) | FLAG(obj2.magerr_galfit)
			    | FLAG(obj2.fluxerr_galfit);

/*---------------------------- External flags -------------------------------*/
  VECFLAG(obj.imaflag) |= VECFLAG(obj.imanflag);

/*------------------------------ PSF-fitting --------------------------------*/
  FLAG(obj2.fwhm_psf) |= FLAG(obj2.fwhmw_psf);
  FLAG(obj2.poserraw_psf) |= FLAG(obj2.poserrbw_psf);
  FLAG(obj2.poserrcxxw_psf) |= FLAG(obj2.poserrcyyw_psf)
			| FLAG(obj2.poserrcxyw_psf);
  FLAG(obj2.poserrthetas_psf) |= FLAG(obj2.poserrtheta1950_psf)
				| FLAG(obj2.poserrtheta2000_psf);
  FLAG(obj2.poserrthetaw_psf) |= FLAG(obj2.poserrthetas_psf);

  FLAG(obj2.poserrmx2w_psf) |= FLAG(obj2.poserrmy2w_psf)
			| FLAG(obj2.poserrmxyw_psf)
			| FLAG(obj2.poserrthetaw_psf) | FLAG(obj2.poserraw_psf)
			| FLAG(obj2.poserrcxxw_psf);

  FLAG(obj2.poserra_psf) |= FLAG(obj2.poserrb_psf)
			| FLAG(obj2.poserrtheta_psf);
  FLAG(obj2.poserrcxx_psf) |= FLAG(obj2.poserrcyy_psf)
			| FLAG(obj2.poserrcxy_psf);
  FLAG(obj2.poserrmx2_psf) |= FLAG(obj2.poserrmy2_psf)
			| FLAG(obj2.poserrmxy_psf)
			| FLAG(obj2.poserra_psf) | FLAG(obj2.poserrcxx_psf)
			| FLAG(obj2.poserrmx2w_psf);

  FLAG(obj2.alpha1950_psf) |= FLAG(obj2.delta1950_psf)
			| FLAG(obj2.poserrtheta1950_psf);
  FLAG(obj2.alpha2000_psf) |= FLAG(obj2.delta2000_psf)
			| FLAG(obj2.alpha1950_psf)
			| FLAG(obj2.poserrtheta2000_psf);
  FLAG(obj2.alphas_psf) |= FLAG(obj2.deltas_psf) | FLAG(obj2.alpha2000_psf);

  FLAG(obj2.xw_psf) |= FLAG(obj2.yw_psf) | FLAG(obj2.poserrmx2w_psf)
			| FLAG(obj2.alphas_psf);

  FLAG(obj2.fluxerr_psf) |= FLAG(obj2.poserrmx2_psf) | FLAG(obj2.magerr_psf);

  FLAG(obj2.mx2_pc) |= FLAG(obj2.my2_pc) | FLAG(obj2.mxy_pc)
			| FLAG(obj2.a_pc) | FLAG(obj2.b_pc)
			| FLAG(obj2.theta_pc) | FLAG(obj2.vector_pc)
			| FLAG(obj2.gdposang) | FLAG(obj2.gdscale)
			| FLAG(obj2.gdaspect) | FLAG(obj2.flux_galfit)
			| FLAG(obj2.gde1) | FLAG(obj2.gde2)
			| FLAG(obj2.gbposang) | FLAG(obj2.gbscale)
			| FLAG(obj2.gbaspect) | FLAG(obj2.gbratio);

  FLAG(obj2.flux_psf) |= FLAG(obj2.mag_psf) | FLAG(obj2.x_psf)
			| FLAG(obj2.y_psf) | FLAG(obj2.xw_psf)
			| FLAG(obj2.fluxerr_psf)
			| FLAG(obj2.niter_psf)
			| FLAG(obj2.chi2_psf)
			| FLAG(obj2.mx2_pc);

/*-------------------------------- Others -----------------------------------*/
  FLAG(obj.fwhm) |= FLAG(obj2.fwhmw);

  FLAG(obj.iso[0]) |= FLAG(obj2.sprob);
  for (i=0; i<NISO; i++)
    FLAG(obj.iso[0]) |= FLAG(obj.iso[i]);

  FLAG(obj.wflag) |= FLAG(obj.nzwpix) | FLAG(obj.nzdwpix);

  return; 
  }


/********************************** dumpparams *******************************/
/*
Initialize the catalog header
*/
void	dumpparams(void)
  {
   int		i;

  for (i=0; *objkey[i].name ; i++)
    if (*objkey[i].unit)
      printf("#%-24.24s %-57.57s [%s]\n",
		objkey[i].name, objkey[i].comment, objkey[i].unit);
    else
      printf("#%-24.24s %-57.57s\n",
		objkey[i].name, objkey[i].comment);

  return;
  }


/********************************** initcat **********************************/
/*
Initialize the catalog header
*/
void	initcat(void)
  {
   keystruct	*key;
   int		i, n;

  if (prefs.cat_type == CAT_NONE)
    return;

  update_tab(objtab);
  if (prefs.cat_type == ASCII_HEAD || prefs.cat_type == ASCII ||
	prefs.cat_type == ASCII_SKYCAT || prefs.cat_type == ASCII_VO)
    {
    if (prefs.pipe_flag)
      ascfile = stdout;
    else
      if (!(ascfile = fopen(prefs.cat_name, "w+")))
        error(EXIT_FAILURE,"*Error*: cannot open ", prefs.cat_name);
//  setlinebuf(ascfile);
    if (prefs.cat_type == ASCII_HEAD && (key = objtab->key))
      for (i=0,n=1; i++<objtab->nkey; key=key->nextkey)
        {
        if (*key->unit)
          fprintf(ascfile, "# %3d %-22.22s %-58.58s [%s]\n",
		n, key->name,key->comment, key->unit);
        else
          fprintf(ascfile, "# %3d %-22.22s %-58.58s\n",
		n, key->name,key->comment);
        n += key->nbytes/t_size[key->ttype];
        }
    else if (prefs.cat_type == ASCII_SKYCAT && (key = objtab->key))
      {
      if (objtab->nkey<3)
        error(EXIT_FAILURE,"The SkyCat format requires at least 4 parameters:",
	      " Id Ra Dec Mag");
/*--- We add a tab between rows, as required by Skycat */
      fprintf(ascfile, skycathead, 8.0);
      for (i=1,key=key->nextkey; i++<objtab->nkey; key=key->nextkey)
        {
        if (i>4)
          fprintf(ascfile, "\t%s", key->name);
        sprintf(gstr, "\t%s", key->printf);
        strcpy(key->printf, gstr);
        }
      fprintf(ascfile, "\n------------------\n");
      }
    else if (prefs.cat_type == ASCII_VO && objtab->key) 
      {
      write_xml_header(ascfile);
      write_vo_fields(ascfile);
      fprintf(ascfile, "   <DATA><TABLEDATA>\n");
      }
    }
  else
    {
    fitscat = new_cat(1);
    init_cat(fitscat);
    strcpy(fitscat->filename, prefs.cat_name);
    if (open_cat(fitscat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE,"*Error*: cannot open for writing ",prefs.cat_name);
    switch(prefs.cat_type)
      {
      case FITS_LDAC:
      case FITS_TPX:
/*------ Save a "pure" primary HDU */
        save_tab(fitscat, fitscat->tab);
        break;

      case FITS_10:
/*------ Add to the primary HDU extraction parameters */
        key = headkey1;
        while (*key->name)
          addkeyto_head(fitscat->tab, key++);
        save_tab(fitscat, fitscat->tab);
        break;
      default:
        error (EXIT_FAILURE, "*Internal Error*: Unknown FITS type in ",
		"initcat()");
      }
    }

  catopen_flag = 1;

  return;
  }


/****** write_vo_fields *******************************************************
PROTO	int	write_vo_fields(FILE *file)
PURPOSE	Write the list of columns to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
void	write_vo_fields(FILE *file)
  {
   keystruct	*key;
   char		datatype[40], arraysize[40], str[40];
   int		i, d;

  if (!objtab || !objtab->key)
    return;
  key=objtab->key;
  for (i=0; i++<objtab->nkey; key=key->nextkey)
    {
/*--- indicate datatype, arraysize, width and precision attributes */
/*--- Handle multidimensional arrays */
    arraysize[0] = '\0';
    if (key->naxis>1)
      {
      for (d=0; d<key->naxis; d++)
        {
        sprintf(str, "%s%d", d?"x":" arraysize=\"", key->naxisn[d]);
        strcat(arraysize, str);
        }
      strcat(arraysize, "\"");
      }
    switch(key->ttype)
      {
      case T_BYTE:	strcpy(datatype, "unsignedByte"); break;
      case T_SHORT:	strcpy(datatype, "short"); break;
      case T_LONG:	strcpy(datatype, "int"); break;
      case T_FLOAT:	strcpy(datatype, "float"); break;
      case T_DOUBLE:	strcpy(datatype, "double"); break;
      default:		error(EXIT_FAILURE,
			"*Internal Error*: Unknown datatype in ",
			"initcat()");
      }
    fprintf(file,
	"  <FIELD name=\"%s\" ucd=\"%s\" datatype=\"%s\" unit=\"%s\"%s>\n",
	key->name, key->voucd, datatype,key->vounit, arraysize);
    fprintf(file, "   <DESCRIPTION>%s</DESCRIPTION>\n", key->comment);
    fprintf(file, "  </FIELD>\n");
    }

  return;
  }


/********************************* reinitcat *********************************/
/*
Initialize the catalog header
*/
void	reinitcat(picstruct *field)
  {
   tabstruct	*tab, *asctab;
   keystruct	*key;

  if (prefs.cat_type == CAT_NONE)
    return;

  if (prefs.cat_type != ASCII_HEAD && prefs.cat_type != ASCII &&
	prefs.cat_type != ASCII_SKYCAT && prefs.cat_type != ASCII_VO)
    {
    update_tab(objtab);
    switch(prefs.cat_type)
      {
      case FITS_LDAC:
/*------ We create a dummy table (only used through its header) */
        QCALLOC(asctab, tabstruct, 1);
        asctab->headnblock = field->tab->headnblock;
        QMEMCPY(field->tab->headbuf, asctab->headbuf, char,
		asctab->headnblock*FBSIZE);
        key = headkey;
        while (*key->name)
          addkeyto_head(asctab, key++);
        tab = new_tab("LDAC_IMHEAD");
        add_tab(tab, fitscat, 0);
        key = new_key("Field Header Card");
        key->ptr = asctab->headbuf;
        asctab->headbuf = NULL;
        free_tab(asctab);
        key->htype = H_STRING;
        key->ttype = T_STRING;
        key->nobj = 1;
        key->nbytes = 80*(fitsfind(key->ptr, "END     ")+1);
        key->naxis = 2;
        QMALLOC(key->naxisn, int, key->naxis);
        key->naxisn[0] = 80;
        key->naxisn[1] = key->nbytes/80;
        add_key(key, tab, 0);
        save_tab(fitscat, tab);
        strcpy(objtab->extname, "LDAC_OBJECTS");
        break;

      case FITS_TPX:
/*------ We create a dummy table (only used through its header) */
        QCALLOC(asctab, tabstruct, 1);
        asctab->headnblock = field->tab->headnblock;
        QMEMCPY(field->tab->headbuf, asctab->headbuf, char,
		asctab->headnblock*FBSIZE);
        key = headkey;
        while (*key->name)
          addkeyto_head(asctab, key++);
        tab = new_tab("TPX_IMHEAD");
        add_tab(tab, fitscat, 0);
        key = new_key("Field Header Card");
        key->ptr = asctab->headbuf;
        asctab->headbuf = NULL;
        free_tab(asctab);
        key->htype = H_STRING;
        key->ttype = T_STRING;
        key->nobj = fitsfind(key->ptr, "END     ")+1;
        key->nbytes = 80;
        key->naxis = 1;
        QMALLOC(key->naxisn, int, key->naxis);
        key->naxisn[0] = 80;
        add_key(key, tab, 0);
        save_tab(fitscat, tab);
        strcpy(objtab->extname, "TPX_OBJECTS");
        break;

      case FITS_10:
/*------ Add to the primary HDU extraction parameters */
/*
        key = headkey1;
        while (*key->name)
          addkeyto_head(fitscat->tab, key++);
        save_tab(fitscat, fitscat->tab);
*/
        break;

      default:
        error (EXIT_FAILURE, "*Internal Error*: Unknown FITS type in ",
		"reinitcat()");
      }

    objtab->cat = fitscat;
    init_writeobj(fitscat, objtab, &buf);
    }

  zerocat();

  return;
  }


/********************************* writecat **********************************/
/*
Write out in the catalog each one object.
*/
void	writecat(int n, objliststruct *objlist)
  {
  outobj = objlist->obj[n];

  switch(prefs.cat_type)
    {
    case FITS_10:
    case FITS_LDAC:
    case FITS_TPX:
      write_obj(objtab, buf);
      break;

    case ASCII:
    case ASCII_HEAD:
    case ASCII_SKYCAT:
      print_obj(ascfile, objtab);
      break;
    case ASCII_VO:
      voprint_obj(ascfile, objtab);
      break;

    case CAT_NONE:
      break;

    default:
      error (EXIT_FAILURE, "*Internal Error*: Unknown catalog type", "");
    }

  return;
  }


/********************************** endcat ***********************************/
/*
Terminate the catalog output.
*/
void	endcat(char *error)
  {
   keystruct	*key;
   int		i;

  if (!catopen_flag)
    {
    if (prefs.cat_type == ASCII_VO)
      write_xmlerror(prefs.cat_name, error);
    return;
    }
  switch(prefs.cat_type)
    {
    case ASCII:
    case ASCII_HEAD:
      if (!prefs.pipe_flag)
        fclose(ascfile);
      break;

    case ASCII_SKYCAT:
      fprintf(ascfile, "%s", skycattail);
      if (!prefs.pipe_flag)
        fclose(ascfile);
      break;

    case ASCII_VO:
      fprintf(ascfile, "    </TABLEDATA></DATA>\n");
      fprintf(ascfile, "  </TABLE>\n");
/*---- Add configuration file meta-data */
      write_xml_meta(ascfile, error);
      fprintf(ascfile, "</RESOURCE>\n");
      fprintf(ascfile, "</VOTABLE>\n");

      if (!prefs.pipe_flag)
        fclose(ascfile);
      break;

    case FITS_LDAC:
    case FITS_TPX:
    case FITS_10:
      free_cat(&fitscat,1);
      break;

    case CAT_NONE:
      break;

    default:
      break;
    }

/* Free allocated memory for arrays */
  key = objtab->key;
  for (i=objtab->nkey; i--; key=key->nextkey)
    if (key->naxis && key->allocflag)
      free(key->ptr);

  objtab->key = NULL;
  objtab->nkey = 0;
  free_tab(objtab);
  objtab = NULL;

  return;
  }


/******************************** reendcat ***********************************/
/*
Terminate the catalog output.
*/
void	reendcat()
  {
   keystruct	*key;
   tabstruct	*tab;
   OFF_T2	pos;
   char		*head;

  switch(prefs.cat_type)
    {
    case ASCII:
    case ASCII_HEAD:
    case ASCII_SKYCAT:
    case ASCII_VO:
      break;

    case FITS_LDAC:
    case FITS_TPX:
      end_writeobj(fitscat, objtab, buf);
      key = NULL;
      if (!(tab=fitscat->tab->prevtab)
	|| !(key=name_to_key(tab, "Field Header Card")))
        error(EXIT_FAILURE,"*Error*: cannot update table ", "ASCFIELD");
      head = key->ptr;
      fitswrite(head, "SEXNDET ", &thecat.ndetect,H_INT,T_LONG);
      fitswrite(head, "SEXNFIN ", &thecat.ntotal, H_INT,T_LONG);
      fitswrite(head, "SEXDATE ", thecat.ext_date, H_STRING, T_STRING);
      fitswrite(head, "SEXTIME ", thecat.ext_time, H_STRING, T_STRING);
      fitswrite(head, "SEXELAPS", &thecat.ext_elapsed, H_FLOAT, T_DOUBLE);
      QFTELL(fitscat->file, pos, fitscat->filename);
      QFSEEK(fitscat->file, tab->headpos, SEEK_SET, fitscat->filename);
      save_tab(fitscat, tab);
      remove_tab(fitscat, "LDAC_IMHEAD", 0);
      QFSEEK(fitscat->file, pos, SEEK_SET, fitscat->filename);
      break;

    case FITS_10:
      end_writeobj(fitscat, objtab, buf);
      fitswrite(fitscat->tab->headbuf,"SEXNDET ",&thecat.ndetect,H_INT,T_LONG);
      fitswrite(fitscat->tab->headbuf,"SEXNFIN ",&thecat.ntotal, H_INT,T_LONG);
      QFTELL(fitscat->file, pos, fitscat->filename);
      QFSEEK(fitscat->file, fitscat->tab->headpos, SEEK_SET,fitscat->filename);
      save_tab(fitscat, fitscat->tab);
      QFSEEK(fitscat->file, pos, SEEK_SET, fitscat->filename);
      break;

    case CAT_NONE:
      break;

    default:
      break;
    }

  return;
  }


/********************************** zerocat **********************************/
/*
Reset to 0 all measurements in the current catalog.
*/
void	zerocat(void)
  {
   keystruct	*key;
   int		i;

  key = objtab->key;
  for (i=objtab->nkey; i--; key=key->nextkey)
    memset(key->ptr, 0, key->nbytes);

  return;
  }


