/*
*				wcscelsys.h
*
* Celestial system definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1998-2015 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		23/11/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------- constants --------------------------------*/

/* Equatorial coordinates of origin and pole and rotation sign of equatorial,*/
/* galactic, ecliptic and supergalactic reference frames, from Allen Astron. */
/* Quantities, 4th ed. */

char	celsysname[][2][8] = {  {"RA--", "DEC-"},
				{"GLON", "GLAT"},
				{"ELON", "ELAT"},
				{"SLON", "SLAT"},
				{""}};
double	celsysorig[][2] = {	{0.0, 0.0},
				{266.40499625, -28.93617242},
				{0.0, 0.0},
				{42.308333, 59.528333}},
	celsyspole[][2] = {	{0.0, 90.0},
				{192.85948123, 27.12825120},
				{270.00000000, 66.560709},
				{283.754167, 15.708889}},
/* Note: the code to handle the rotation sign is not yet implemented!!! */
	celsyssign[]	= {	 1.0,
				 1.0,
				 1.0,
				 1.0};

