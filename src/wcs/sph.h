/*
*				sph.h
*
* Include file for sph.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 2000-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*				(C) 1995-1999 Mark Calabretta (original version)
*
*	Licenses:		GNU General Public License
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
/*=============================================================================
*
*   WCSLIB - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995-1999, Mark Calabretta
*
*   This library is free software; you can redistribute it and/or modify it
*   under the terms of the GNU Library General Public License as published
*   by the Free Software Foundation; either version 2 of the License, or (at
*   your option) any later version.
*
*   This library is distributed in the hope that it will be useful, but
*   WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library
*   General Public License for more details.
*
*   You should have received a copy of the GNU Library General Public License
*   along with this library; if not, write to the Free Software Foundation,
*   Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*   Correspondence concerning WCSLIB may be directed to:
*      Internet email: mcalabre@atnf.csiro.au
*      Postal address: Dr. Mark Calabretta,
*                      Australia Telescope National Facility,
*                      P.O. Box 76,
*                      Epping, NSW, 2121,
*                      AUSTRALIA
*
*   Author: Mark Calabretta, Australia Telescope National Facility
*   $Id: sph.h,v 1.1.1.1 2002/03/15 16:33:26 bertin Exp $
*===========================================================================*/

#ifndef WCSLIB_SPH
#define WCSLIB_SPH

#ifdef __cplusplus
extern "C" {
#endif

#if __STDC__  || defined(__cplusplus)
   int sphfwd(const double, const double,
              const double [],
              double *, double *);
   int sphrev(const double, const double,
              const double [],
              double *, double *);
#else
   int sphfwd(), sphrev();
#endif

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_SPH */
