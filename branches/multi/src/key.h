/*
*				key.h
*
* Keyword structure definition.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------------- constants -------------------------------*/

#define         FIND_STRICT     0
#define         FIND_NOSTRICT   1

/*--------------------------- structure definitions -------------------------*/
/* Preference keyword */
typedef struct
  {
  char		name[32];
  enum  {P_FLOAT, P_INT, P_STRING, P_BOOL, P_KEY, P_INTLIST, P_FLOATLIST,
	P_BOOLLIST, P_KEYLIST, P_STRINGLIST} type;
  void		*ptr;			/* Pointer to the keyword value */
  int		imin, imax;		/* Range for int's */
  double	dmin, dmax;		/* Range for floats */
  char		keylist[32][32];	/* List of keywords */
  int           nlistmin;		/* Minimum number of list members */
  int           nlistmax; 		/* Maximum number of list members */
  int		*nlistptr;		/* Ptr to store the nb of read params*/
  int		flag;
  }	pkeystruct;

/*---------------------------------- protos --------------------------------*/

int	findkeys(char *str, char key[][32], int mode);

