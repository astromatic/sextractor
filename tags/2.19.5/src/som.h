/*
*				som.h
*
* Include file for som.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1995-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------------- constants ------------------------------*/

#define	INPUT_MAXDIM		9	/* Maximum dimensionality of input */
#define	SOM_MAXDIM		6	/* Maximum dimensionality of the SOM */

/*------------------------------- SOM flags --------------------------------*/

#define		SOM_NODE	0x01	/* Compute at some exact node pos */
#define		SOM_PHOTOM	0x02	/* Do photometry */
#define		SOM_GRADIENT	0x04	/* Compute interpolated SOM gradient */
#define		SOM_LINE	0x08	/* Proceed along a specific line */

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  int		inputdim;		/* Dimensionality of input vector */
  int		*inputsize;		/* Dimensions of the input vector */
  int		ninput;			/* Total number of inputs */
  int		nextrainput;		/* Number of extra inputs */
  int		neurdim;		/* Dimensionality of the SOM */
  int		*neursize;		/* Dimensions of the SOM */
  int		nneur;			/* Total number of neurons */
  int		*neurstep;		/* Stepping through the SOM */
  float		*weight;		/* Weights */
  int		nweight;		/* Total number of weights */
  float		*input;			/* Input data */
  float		*inputw;		/* Input data weighting */
  float		*proto;			/* Current composite prototype */
  float		*dproto;		/* Current composite gradients */
  float		*vector;		/* Current SOM coordinates */
  float		*dvector;		/* Current SOM search direction */
  float		learnrate, clearnrate;	/* Starting and current learn. rates */
  float		learndecay;		/* Learning decay rate */
  float		kernw, ckernw;		/* Starting and current kernel width */
  float		kernwdecay;		/* Kernel width decay rate */
  float		xy_stiff;		/* Stiffness of the X/Y mapping */
  int		*freq;			/* Number of winning times per node */
  int		ntrain;			/* # of training examples so far */
  int		nsweep;			/* # of sweeps through the whole set */
  float		amp, sigamp;		/* Best fitting amplitude and error */
  float		stderror;		/* Global reduced error */
  }	somstruct;

somstruct	*thesom;

/*---------------------------------- protos --------------------------------*/

extern somstruct	*som_load(char *filename);

extern float		som_err(somstruct *som, float dist, int flag),
			som_linmin(somstruct *som);

extern int		som_mkweight(somstruct *som,float back,float backnoise,
				float gain);

extern void		som_conjgrad(somstruct *som, float ftol),
			som_end(somstruct *som),
			som_phot(somstruct *som, float back,float backnoise,
				float gain, float dx, float dy,
				float *vector, float clip),
			som_start(somstruct *som, float *context,
				int ncontext, float x, float y);
