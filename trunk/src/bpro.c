/*
*				bpro.c
*
* Manage backpropagation neural networks.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"bpro.h"

/******************************** play_bpann *********************************/
/*
Single forward pass through the ANN.
*/
void	play_bpann(bpannstruct *bpann, NFLOAT *invec, NFLOAT *outvec)
  {
   NFLOAT	u, *neuroni,*neuronj, *weight;
   int		i,j,l,lp,ll, lflag;

  ll = bpann->nlayers-1;
  memcpy(bpann->neuron[0], invec, bpann->nn[0]*sizeof(float));
  lflag = bpann->linearoutflag;
  for (lp=0,l=1; lp<ll; l++, lp++)
    {
    neuronj = bpann->neuron[l];
    weight = bpann->weight[lp];
    for (j=bpann->nn[l]; j--; neuronj++)
      {	/* note we don't touch the "bias" neuron (=-1) */
      neuroni = bpann->neuron[lp];
      u = *(weight++)**(neuroni++);
      for (i=bpann->nn[lp]; i--;)	/* The last one is the bias */
        u += *(weight++)**(neuroni++);
      if (l == ll)
        *(outvec++)= lflag?u:SIGMOID(u);
      else
        *neuronj = SIGMOID(u);
      }
    }

  return;
  }


/******************************* loadtab_bpann *******************************/
/*
Load the relevant ANN structure (using the LDACTools).
*/
bpannstruct	*loadtab_bpann(tabstruct *tab, char *filename)
  {
   bpannstruct	*bpann;
   keystruct	*key;
   char		*head, str[80];
   int		l;

/* OK, we now allocate memory for the ANN structure itself */
  QCALLOC(bpann, bpannstruct, 1);
/* Load important scalars (which are stored as FITS keywords) */
  head = tab->headbuf;
  if (fitsread(head, "BPNLAYER", &bpann->nlayers, H_INT, T_LONG) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: incorrect BP-ANN header in ", filename);
  if (fitsread(head, "BPLINEAR",&bpann->linearoutflag, H_INT,T_LONG)!=RETURN_OK)
    bpann->linearoutflag = 0;
/* Load all vectors!! */
  read_keys(tab, NULL, NULL, 0, NULL);
/* Now interpret the result */
  if (!(key = name_to_key(tab, "NNEUR_PER_LAYER")))
    error(EXIT_FAILURE, "*Error*: incorrect BP-ANN header in ", filename);
  bpann->nn = key->ptr; key->ptr = 0;
  QMALLOC(bpann->neuron, NFLOAT *, bpann->nlayers);
  QMALLOC(bpann->weight, NFLOAT *, bpann->nlayers-1);
  for (l=0; l<bpann->nlayers-1; l++)
    {
    QMALLOC(bpann->neuron[l], NFLOAT, bpann->nn[l]+1);
    bpann->neuron[l][bpann->nn[l]] = -1.0;
    sprintf(str, "WEIGHT_LAYER%d", l+1);
    if (!(key = name_to_key(tab, str)))
      error(EXIT_FAILURE, "*Error*: incorrect BP-ANN header in ", filename);
    bpann->weight[l] = key->ptr; key->ptr = 0;
    }

  QMALLOC(bpann->neuron[l], NFLOAT, bpann->nn[l]); /* no bias in this layer */

  return bpann;
  }


/******************************** free_bpann *********************************/
/*
Free all memory modules allocated for a Back-Propagation ANN structure.*/
void    free_bpann(bpannstruct *bpann)

  {
   int          i;

/* Loop over the "true" layers */
  for (i=0; i<bpann->nlayers-1; i++)
    {
    free(bpann->neuron[i]);
    free(bpann->weight[i]);
    }

  free(bpann->neuron[i]);       /* Because of the input layer */

/* Then free pointers of pointers */
  free(bpann->neuron);
  free(bpann->weight);
  free(bpann->nn);

/* And finally free the ANN structure itself */
  free(bpann);

  return;
  }


