/*
*				field.h
*
* Include file for field.c.
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

/*------------------------------ field flags -------------------------------*/
#define		DETECT_FIELD	0x01	/* Detection */
#define		MEASURE_FIELD	0x02	/* Measurement */
#define		FLAG_FIELD	0x04	/* Flagging */
#define		RMS_FIELD	0x08	/* Weighting with std deviations */
#define		VAR_FIELD	0x10	/* Weighting with variances */
#define		WEIGHT_FIELD	0x20	/* Weighting with weights */
#define		BACKRMS_FIELD	0x40	/* Weighting from a backrms matrix */
#define		INTERP_FIELD	0x80	/* Purely interpolated data */


