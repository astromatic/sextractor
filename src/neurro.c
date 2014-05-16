/*
*				neurro.c
*
* Read-only implementation of a backprogation neural network.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		28/03/2012
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
#include	"neurro.h"

static double neur_activ(double x);

brainstruct	*brain;

/****** neur_init ***********************************************************
PROTO	void neur_init(void)
PURPOSE	Initialize the neural network ("brain").
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2012
 ***/
void	neur_init(void)

  {
  QMALLOC(brain, brainstruct, 1);

  return;
  }


/****** neur_end ***********************************************************
PROTO	void neur_end(void)
PURPOSE	End the neural network ("brain").
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2012
 ***/
void	neur_end(void)

  {
  free(brain);

  return;
  }


/****** neur_resp ***********************************************************
PROTO	void neur_resp(double *input, double *output)
PURPOSE	Compute the neural network response to an input vector.
INPUT	Pointer to the input vector,
	pointer to the output vector (filled with new values in output).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2012
 ***/
void	neur_resp(double *input, double *output)

  {
   int	i, j, l, lastlay = brain->layersnb-1;
   double	neursum;

  for (i=0; i<brain->nn[0]; i++)
    brain->n[0][i] = input[i]*brain->inscale[i] + brain->inbias[i];
  for (l=0; l<lastlay; l++)
    for (j=0; j<brain->nn[l+1]; j++)
      {
      neursum = brain->b[l][j];
      for (i=0; i<brain->nn[l]; i++)
        neursum += brain->w[l][i][j] * brain->n[l][i];
      brain->n[l+1][j] = neur_activ(neursum);
      }
  for (i=0; i<brain->nn[lastlay]; i++)
    output[i] = (brain->n[lastlay][i]-brain->outbias[i])
	/ brain->outscale[i];

  return;
  }


/*i**** neur_activ ***********************************************************
PROTO	static double neur_activ(double x)
PURPOSE	Activation function (sigmoid function) of the input scalar.
INPUT	Input value.
OUTPUT	Output activation value.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2012
 ***/
static double	neur_activ(double x)

  {
  return 1.0 / (1.0 + exp(-x));
  }


/****** neur_getnnw *********************************************************
PROTO	void neur_getnnw(char *filename)
PURPOSE	Read the ".nnw" ASCII file containing the neural network weights.
INPUT	File name.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2012
 ***/
void	neur_getnnw(char *filename)
  {

   FILE	*infile;
   int	i, j, k, step;
   char	str[MAXCHAR], *sstr, *null;

  if ((infile = fopen(filename,"r")) == NULL)
    error(EXIT_FAILURE,"*ERROR*: can't read ", filename);

  fgets(str, MAXCHAR, infile);
  if (strncmp(str,"NNW",3))
    error(EXIT_FAILURE, filename, " is NOT a NNW table!");

  step = 1;
  i=j=0;			/* To avoid gcc -Wall warnings */
  while (fgets(str, MAXCHAR, infile))
    {
    sstr = &str[(int)strspn(str," \t")];
    if (sstr[0]!=(char)'#' && sstr[0]!=(char)'\n')
      {
      null = sstr;
      switch(step)
        {
        case 1:	brain->layersnb = atoi(strtok(sstr, " \t\n"));
		for (i=0; i<brain->layersnb; i++)
		  brain->nn[i] = atoi(strtok(NULL, " \t\n"));
		step++;
		break;

        case 2:	for (i=0; i<brain->nn[0]; i++)
		  {
		  brain->inbias[i] = atof(strtok(null, " \t\n"));
		  null = NULL;
		  }
		step++;
		break;

        case 3:	for (i=0; i<brain->nn[0]; i++)
		  {
		  brain->inscale[i] = atof(strtok(null, " \t\n"));
		  null = NULL;
		  }
		i=j=0;
		step++;
		break;

        case 4:	if (j == brain->nn[i+1])
		  {
		  j = 0;
		  i++;
		  }
		if (i < brain->layersnb-1)
		  {
		  for (k=0; k<brain->nn[i]; k++)
		    {
		    brain->w[i][k][j] = atof(strtok(null, " \t\n"));
		    null = NULL;
		    }
		  brain->b[i][j] = atof(strtok(NULL, " \t\n"));
		  j++;
		  break;
		  }
		else
		  step++;

        case 5:	for (i=0; i<brain->nn[brain->layersnb-1]; i++)
		  {
		  brain->outbias[i] = atof(strtok(null, " \t\n"));
		  null = NULL;
		  }
		step++;
		break;
        case 6:	for (i=0; i<brain->nn[brain->layersnb-1]; i++)
		  {
		  brain->outscale[i] = atof(strtok(null, " \t\n"));
		  null = NULL;
		  }
		step++;
		break;
	default:error(EXIT_FAILURE, "*Error*: inconsistency in ",
		filename);
	}

      }
    }

  fclose(infile);

  return;
  }

