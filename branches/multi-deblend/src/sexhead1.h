/*
*				sexhead1.h
*
* Keyword list for FITS_1.0 catalogue headers.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1996-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		21/03/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int	idummy;
double	ddummy;

keystruct	headkey1[] = {
  {"EPOCH   ", "",
	&thefield.epoch, H_FLOAT, T_DOUBLE, "%7.2f"},
  {"OBJECT  ", "",
	thefield.ident, H_STRING, T_STRING, "%18s"},
  {"ORIGIN  ", "",
	"SExtractor", H_STRING, T_STRING, "%18s"},
  {"CRVAL1", "WORLD X COORD. OF REFERENCE PIXEL",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"CRVAL2", "WORLD Y COORD. OF REFERENCE PIXEL",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"CRPIX1", "IMAGE X COORD. OF REFERENCE PIXEL",
	&idummy, H_INT, T_LONG, "%5d"},
  {"CRPIX2", "IMAGE Y COORD. OF REFERENCE PIXEL",
	&idummy, H_INT, T_LONG, "%5d"},
  {"CDELT1", "WORLD PIXEL STEP ALONG X",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"CDELT2", "WORLD PIXEL STEP ALONG Y",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"CROTA1", "CCW ANGLE FROM X-IMAGE TO X-WORLD",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"CROTA2", "CCW ANGLE FROM Y-IMAGE TO Y-WORLD",
	&ddummy, H_EXPO, T_DOUBLE, "%15G"},
  {"FITSFILE", "File name of the analysed image",
	thecat.image_name, H_STRING, T_STRING, "%-18s"},
  {"FITSEXT ", "FITS Extension number",
	&thecat.currext, H_INT, T_LONG, "%3d"},
  {"FITSNEXT", "Number of FITS image extensions in file",
	&thecat.next, H_INT, T_LONG, "3d"},
  {"SEXIMASX", "IMAGE WIDTH (PIXELS)",
	&thefield.width, H_INT, T_LONG, "%10d"},
  {"SEXIMASY", "IMAGE HEIGHT (PIXELS)",
	&thefield.height, H_INT, T_LONG, "%10d"},
  {"SEXSTRSY", "STRIP HEIGHT (LINES)",
	&thefield.stripheight, H_INT, T_LONG, "%10d"},
  {"SEXIMABP", "FITS IMAGE BITPIX",
	&thefield.bitpix, H_INT, T_LONG, "%3d"},
  {"SEXPIXS", "PIXEL SCALE (ARCSEC)",
	&thefield.pixscale, H_EXPO, T_DOUBLE, "%12G"},
  {"SEXSFWHM", "SEEING FWHM (ARCSEC)",
	&prefs.seeing_fwhm, H_EXPO, T_DOUBLE, "%12G"},
  {"SEXNNWF ", "CLASSIFICATION NNW FILENAME",
	thecat.nnw_name, H_STRING, T_STRING, "%18s"},
  {"SEXGAIN ", "GAIN (IN E- PER ADU)",
	&thefield.gain, H_EXPO, T_DOUBLE, "%7.3F"},
  {"SEXBKGND", "MEDIAN BACKGROUND (ADU)",
	&thefield.backmean, H_EXPO, T_FLOAT, "%12G"},
  {"SEXBKDEV", "MEDIAN RMS (ADU)",
	&thefield.backsig, H_EXPO, T_FLOAT, "%12G"},
  {"SEXBKTHD", "EXTRACTION THRESHOLD (ADU)",
	&thefield.thresh, H_EXPO, T_FLOAT, "%12G"},
  {"SEXCONFF", "CONFIGURATION FILENAME",
	thecat.prefs_name, H_STRING, T_STRING, "%18s"},
  {"SEXDETT ", "DETECTION TYPE",
	"CCD", H_STRING, T_STRING, "%s"},
  {"SEXTHLDT", "THRESHOLD TYPE",
	"SIGMA", H_STRING, T_STRING, "%s"},
  {"SEXTHLD ", "THRESHOLD",
	&prefs.dthresh[0], H_EXPO, T_DOUBLE, "%12G"},
  {"SEXMINAR", "EXTRACTION MINIMUM AREA (PIXELS)",
	&prefs.ext_minarea, H_INT, T_LONG, "%6d"},
  {"SEXCONV ", "CONVOLUTION FLAG",
	&prefs.filter_flag, H_BOOL, T_LONG, "%1s"},
  {"SEXCONVN", "CONVOLUTION NORM. FLAG",
	&prefs.filter_flag, H_BOOL, T_LONG, "%1s"},
  {"SEXCONVF", "CONVOLUTION FILENAME",
	thecat.filter_name, H_STRING, T_STRING, "%18s"},
  {"SEXDBLDN", "NUMBER OF SUB-THRESHOLDS",
	&prefs.deblend_nthresh, H_INT, T_LONG, "%3d"},
  {"SEXDBLDC", "CONTRAST PARAMETER",
	&prefs.deblend_mincont, H_FLOAT, T_DOUBLE, "%8f"},
  {"SEXCLN  ", "CLEANING FLAG",
	&prefs.clean_flag, H_BOOL, T_LONG, "%1s"},
  {"SEXCLNPA", "CLEANING PARAMETER",
	&prefs.clean_param, H_FLOAT, T_DOUBLE, "%8f"},
  {"SEXCLNST", "CLEANING OBJECT-STACK",
	&prefs.deblend_nthresh, H_INT, T_LONG, "%6d"},
  {"SEXAPERD", "1ST APERTURE DIAMETER (PIXELS)",
	&prefs.apert[0], H_FLOAT, T_DOUBLE, "%8.2f"},
  {"SEXAPEK1", "KRON PARAMETER",
	&prefs.autoparam[0], H_FLOAT, T_DOUBLE, "%4.1f"},
  {"SEXAPEK2", "KRON ANALYSIS RADIUS",
	&prefs.autoparam[0], H_FLOAT, T_DOUBLE, "%4.1f"},
  {"SEXAPEK3", "KRON MINIMUM RADIUS",
	&prefs.autoparam[1], H_FLOAT, T_DOUBLE, "%4.1f"},
  {"SEXSATLV", "SATURATION LEVEL (ADU)",
	&thefield.satur_level, H_EXPO, T_DOUBLE, "%12G"},
  {"SEXMGZPT", "MAGNITUDE ZERO-POINT",
	&prefs.mag_zeropoint, H_FLOAT, T_DOUBLE, "%8.4f"},
  {"SEXMGGAM", "MAGNITUDE GAMMA",
	&prefs.mag_gamma, H_FLOAT, T_DOUBLE, "%4.2f"},
  {"SEXBKGSX", "BACKGROUND MESH WIDTH (PIXELS)",
	&thefield.backw, H_INT, T_LONG, "%5d"},
  {"SEXBKGSY", "BACKGROUND MESH HEIGHT (PIXELS)",
	&thefield.backh, H_INT, T_LONG, "%5d"},
  {"SEXBKGFX", "BACKGROUND FILTER WIDTH",
	&thefield.nbackfx, H_INT, T_LONG, "%3d"},
  {"SEXBKGFY", "BACKGROUND FILTER HEIGHT",
	&thefield.nbackfy, H_INT, T_LONG, "%3d"},
  {"SEXPBKGT", "PHOTOM BACKGROUND TYPE",
	"GLOBAL", H_STRING, T_STRING, "%s"},
  {"SEXPBKGS", "LOCAL AREA THICKNESS (PIXELS)",
	&prefs.pback_size, H_INT, T_LONG, "%3d"},
  {"SEXPIXSK", "PIXEL STACKSIZE (PIXELS)",
	&prefs.mem_pixstack, H_INT, T_LONG, "%8d"},
  {"SEXFBUFS", "FRAME-BUFFER SIZE (LINES)",
	&prefs.mem_bufsize, H_INT, T_LONG, "%5d"},
  {"SEXISAPR", "ISO-APER RATIO",
	 &ddummy, H_FLOAT, T_DOUBLE, "%4.2f"},
  {"SEXNDET ", "NB OF DETECTIONS",
	&thecat.ndetect, H_INT, T_LONG, "%9d"},
  {"SEXNFIN ", "NB OF FINAL EXTRACTED OBJECTS",
	&thecat.ntotal, H_INT, T_LONG, "%9d"},
  {"SEXNPARA", "NB OF PARAMETERS PER OBJECT",
	&thecat.nparam, H_INT, T_LONG, "%3d"},
  {""}};

