/*
*				param.h
*
* List of regular measurement parameters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		03/08/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

objstruct	flagobj;
obj2struct	flagobj2;

/*--------------------------------- initialization --------------------------*/
keystruct	objkey[] = {
  {"NUMBER", "Running object number",
	&flagobj2.number, H_INT, T_LONG, "%10d", "",
	"meta.record", ""},
  {"EXT_NUMBER", "FITS extension number",
	&flagobj2.ext_number, H_INT, T_SHORT, "%3d", "",
	"meta.id;meta.dataset", ""},
  {"FLUX_ISO", "Isophotal flux",
	&flagobj2.flux_iso, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_ISO", "RMS error for isophotal flux",
	&flagobj2.fluxerr_iso, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_ISO", "Isophotal magnitude",
	&flagobj2.mag_iso, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_ISO", "RMS error for isophotal magnitude",
	&flagobj2.magerr_iso, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},

  {"FLUX_ISOCOR", "Corrected isophotal flux",
	&flagobj2.flux_isocor, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_ISOCOR", "RMS error for corrected isophotal flux",
	&flagobj2.fluxerr_isocor, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_ISOCOR", "Corrected isophotal magnitude",
	&flagobj2.mag_isocor, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_ISOCOR", "RMS error for corrected isophotal magnitude",
	&flagobj2.magerr_isocor, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},

  {"FLUX_APER", "Flux vector within fixed circular aperture(s)",
	&flagobj2.flux_aper, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct", 1, &prefs.flux_apersize},
  {"FLUXERR_APER", "RMS error vector for aperture flux(es)",
	&flagobj2.fluxerr_aper, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct", 1, &prefs.fluxerr_apersize},
  {"MAG_APER", "Fixed aperture magnitude vector",
	&flagobj2.mag_aper, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag", 1, &prefs.mag_apersize},
  {"MAGERR_APER", "RMS error vector for fixed aperture mag.",
	&flagobj2.magerr_aper, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag", 1, &prefs.magerr_apersize},

  {"FLUX_AUTO", "Flux within a Kron-like elliptical aperture",
	&flagobj2.flux_auto, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux;meta.main", "ct"},
  {"FLUXERR_AUTO", "RMS error for AUTO flux",
	&flagobj2.fluxerr_auto, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux;meta.main", "ct"},
  {"MAG_AUTO", "Kron-like elliptical aperture magnitude",
	&flagobj2.mag_auto, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;meta.main", "mag"},
  {"MAGERR_AUTO", "RMS error for AUTO magnitude",
	&flagobj2.magerr_auto, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag;meta.main", "mag"},

  {"FLUX_PETRO", "Flux within a Petrosian-like elliptical aperture",
	&flagobj2.flux_petro, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_PETRO", "RMS error for PETROsian flux",
	&flagobj2.fluxerr_petro, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_PETRO", "Petrosian-like elliptical aperture magnitude",
	&flagobj2.mag_petro, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_PETRO", "RMS error for PETROsian magnitude",
	&flagobj2.magerr_petro, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},

  {"FLUX_BEST", "Best of FLUX_AUTO and FLUX_ISOCOR",
	&flagobj2.flux_best, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_BEST", "RMS error for BEST flux",
	&flagobj2.fluxerr_best, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_BEST", "Best of MAG_AUTO and MAG_ISOCOR",
	&flagobj2.mag_best, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_BEST", "RMS error for MAG_BEST",
	&flagobj2.magerr_best, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},

  {"FLUX_WIN", "Gaussian-weighted flux",
	&flagobj2.flux_win, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_WIN", "RMS error for WIN flux",
	&flagobj2.fluxerr_win, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_WIN", "Gaussian-weighted magnitude",
	&flagobj2.mag_win, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_WIN", "RMS error for MAG_WIN",
	&flagobj2.magerr_win, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},
  {"SNR_WIN", "Gaussian-weighted SNR",
	&flagobj2.snr_win, H_FLOAT, T_FLOAT, "%10.4g", "",
	"stat.snr", ""},

  {"FLUX_SOMFIT", "Flux derived from SOM fit",
	&flagobj2.flux_somfit, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux", "ct"},
  {"FLUXERR_SOMFIT", "RMS error for SOMFIT flux",
	&flagobj2.fluxerr_somfit, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.flux", "ct"},
  {"MAG_SOMFIT", "Magnitude derived from SOM fit",
	&flagobj2.mag_somfit, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_SOMFIT", "Magnitude error derived from SOM fit",
	&flagobj2.magerr_somfit, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},
  {"ERROR_SOMFIT", "Reduced Chi-square error of the SOM fit",
	&flagobj2.stderr_somfit, H_FLOAT, T_FLOAT, "%10.4g", "",
	"stat.fit.chi2", ""},
  {"VECTOR_SOMFIT", "Position vector of the winning SOM node",
	&flagobj2.vector_somfit, H_FLOAT, T_FLOAT, "%5.2f", "",
	"src.morph.param", "", 1, &prefs.somfit_vectorsize},

  {"KRON_RADIUS", "Kron apertures in units of A or B",
	&flagobj2.kronfactor, H_FLOAT, T_FLOAT, "%5.2f", "",
	"arith.factor;arith.ratio", ""},
  {"PETRO_RADIUS", "Petrosian apertures in units of A or B",
	&flagobj2.petrofactor, H_FLOAT, T_FLOAT, "%5.2f", "",
	"arith.factor;arith.ratio", ""},
  {"BACKGROUND", "Background at centroid position",
	&flagobj.bkg, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"instr.skyLevel", "ct"},
  {"THRESHOLD", "Detection threshold above background",
	&flagobj.dthresh, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"instr.sensitivity;phot.flux.sb", "ct"},
  {"FLUX_MAX", "Peak flux above background",
	&flagobj.peak, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux.sb;stat.max", "ct"},
  {"ISOAREA_IMAGE", "Isophotal area above Analysis threshold",
	&flagobj.npix, H_INT, T_LONG, "%9d", "pixel**2",
	"phys.area", "pix2"},
  {"ISOAREAF_IMAGE", "Isophotal area (filtered) above Detection threshold",
	&flagobj.fdnpix, H_INT, T_LONG, "%9d", "pixel**2",
	"phys.area", "pix2"},

  {"XMIN_IMAGE", "Minimum x-coordinate among detected pixels",
	&flagobj.xmin, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.x;stat.min", "pix"},
  {"YMIN_IMAGE", "Minimum y-coordinate among detected pixels",
	&flagobj.ymin, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.y;stat.min", "pix"},
  {"XMAX_IMAGE", "Maximum x-coordinate among detected pixels",
	&flagobj.xmax, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.x;stat.max", "pix"},
  {"YMAX_IMAGE", "Maximum y-coordinate among detected pixels",
	&flagobj.ymax, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.y;stat.max", "pix"},

  {"XPEAK_IMAGE", "x-coordinate of the brightest pixel",
	&flagobj.peakx, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.x", "pix"},
  {"YPEAK_IMAGE", "y-coordinate of the brightest pixel",
	&flagobj.peaky, H_INT, T_LONG, "%10d", "pixel",
	"pos.cartesian.y", "pix"},
  {"XPEAK_FOCAL", "Focal-plane x coordinate of the brightest pixel",
	&flagobj2.peakxf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.x", ""},
  {"YPEAK_FOCAL", "Focal-plane y coordinate of the brightest pixel",
	&flagobj2.peakyf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.y", ""},
  {"XPEAK_WORLD", "World-x coordinate of the brightest pixel",
	&flagobj2.peakxw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.ra", "deg"},
  {"YPEAK_WORLD", "World-y coordinate of the brightest pixel",
	&flagobj2.peakyw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAPEAK_SKY", "Right ascension of brightest pix (native)",
	&flagobj2.peakalphas, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAPEAK_SKY", "Declination of brightest pix (native)",
	&flagobj2.peakdeltas, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAPEAK_J2000", "Right ascension of brightest pix (J2000)",
	&flagobj2.peakalpha2000, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAPEAK_J2000", "Declination of brightest pix (J2000)",
	&flagobj2.peakdelta2000, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAPEAK_B1950", "Right ascension of brightest pix (B1950)",
	&flagobj2.peakalpha1950, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAPEAK_B1950", "Declination of brightest pix (B1950)",
	&flagobj2.peakdelta1950, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"X_IMAGE", "Object position along x",
	&flagobj2.sposx, H_FLOAT, T_FLOAT, "%11.4f", "pixel",
	"pos.cartesian.x;pos.barycenter;instr.det;meta.main", "pix"},
  {"Y_IMAGE", "Object position along y",
	&flagobj2.sposy, H_FLOAT, T_FLOAT, "%11.4f", "pixel",
	"pos.cartesian.y;pos.barycenter;instr.det;meta.main", "pix"},
  {"X_IMAGE_DBL", "Object position along x (double precision)",
	&flagobj2.posx, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.x;pos.barycenter;instr.det", "pix"},
  {"Y_IMAGE_DBL", "Object position along y (double precision)",
	&flagobj2.posy, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.x;pos.barycenter;instr.det", "pix"},
  {"X_FOCAL", "Barycenter position along focal-plane x axis",
	&flagobj2.mxf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.x", ""},
  {"Y_FOCAL", "Barycenter position along focal-plane y axis",
	&flagobj2.myf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.y", ""},
  {"X_WORLD", "Barycenter position along world x axis",
	&flagobj2.mxw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.ra", "deg"},
  {"Y_WORLD", "Barycenter position along world y axis",
	&flagobj2.myw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.dec", "deg"},
  {"X_MAMA", "Barycenter position along MAMA x axis",
	&flagobj2.mamaposx, H_FLOAT, T_DOUBLE, "%8.1f", "m**(-6)",
	"pos.cartesian.x;instr.det;pos.barycenter", "um"},
  {"Y_MAMA", "Barycenter position along MAMA y axis",
	&flagobj2.mamaposy, H_FLOAT, T_DOUBLE, "%8.1f", "m**(-6)",
	"pos.cartesian.y;instr.det;pos.barycenter", "um"},

  {"ALPHA_SKY", "Right ascension of barycenter (native)",
	&flagobj2.alphas, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;pos.barycenter", "deg"},
  {"DELTA_SKY", "Declination of barycenter (native)",
	&flagobj2.deltas, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;pos.barycenter", "deg"},

  {"ALPHA_J2000", "Right ascension of barycenter (J2000)",
	&flagobj2.alpha2000, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;pos.barycenter;meta.main", "deg"},
  {"DELTA_J2000", "Declination of barycenter (J2000)",
	&flagobj2.delta2000, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;pos.barycenter;meta.main", "deg"},

  {"ALPHA_B1950", "Right ascension of barycenter (B1950)",
	&flagobj2.alpha1950, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;pos.barycenter", "deg"},
  {"DELTA_B1950", "Declination of barycenter (B1950)",
	&flagobj2.delta1950, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;pos.barycenter", "deg"},

  {"X2_IMAGE", "Variance along x",
	&flagobj.mx2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;instr.det", "pix2"},
  {"Y2_IMAGE", "Variance along y",
	&flagobj.my2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;instr.det", "pix2"},
  {"XY_IMAGE", "Covariance between x and y",
	&flagobj.mxy, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;instr.det", "pix2"},
  {"X2_WORLD", "Variance along X-WORLD (alpha)",
	&flagobj2.mx2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},
  {"Y2_WORLD", "Variance along Y-WORLD (delta)",
	&flagobj2.my2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},
  {"XY_WORLD", "Covariance between X-WORLD and Y-WORLD",
	&flagobj2.mxyw, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},

  {"CXX_IMAGE", "Cxx object ellipse parameter",
	&flagobj.cxx, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CYY_IMAGE", "Cyy object ellipse parameter",
	&flagobj.cyy, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CXY_IMAGE", "Cxy object ellipse parameter",
	&flagobj.cxy, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CXX_WORLD", "Cxx object ellipse parameter (WORLD units)",
	&flagobj2.cxxw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},
  {"CYY_WORLD", "Cyy object ellipse parameter (WORLD units)",
	&flagobj2.cyyw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},
  {"CXY_WORLD", "Cxy object ellipse parameter (WORLD units)",
	&flagobj2.cxyw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},

  {"A_IMAGE", "Profile RMS along major axis",
	&flagobj.a, H_FLOAT, T_FLOAT, "%9.3f", "pixel",
	"phys.size.smajAxis;instr.det;meta.main", "pix"},
  {"B_IMAGE", "Profile RMS along minor axis",
	&flagobj.b, H_FLOAT, T_FLOAT, "%9.3f", "pixel",
	"phys.size.sminAxis;instr.det;meta.main", "pix"},
  {"THETA_IMAGE", "Position angle (CCW/x)",
	&flagobj.theta, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;instr.det;meta.main", "deg"},
  {"A_WORLD", "Profile RMS along major axis (world units)",
	&flagobj2.aw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.size.smajAxis;meta.main", "deg"},
  {"B_WORLD", "Profile RMS along minor axis (world units)",
	&flagobj2.bw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.size.sminAxis;meta.main", "deg"},
  {"THETA_WORLD", "Position angle (CCW/world-x)",
	&flagobj2.thetaw, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng", "deg"},
  {"THETA_SKY", "Position angle (east of north) (native)",
	&flagobj2.thetas, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng", "deg"},
  {"THETA_J2000", "Position angle (east of north) (J2000)",
	&flagobj2.theta2000, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng;meta.main", "deg"},
  {"THETA_B1950", "Position angle (east of north) (B1950)",
	&flagobj2.theta1950, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng", "deg"},

  {"ERRX2_IMAGE", "Variance of position along x",
	&flagobj.poserr_mx2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRY2_IMAGE", "Variance of position along y",
	&flagobj.poserr_my2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRXY_IMAGE", "Covariance of position between x and y",
	&flagobj.poserr_mxy, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.covariance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRX2_WORLD", "Variance of position along X-WORLD (alpha)",
	&flagobj2.poserr_mx2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRY2_WORLD", "Variance of position along Y-WORLD (delta)",
	&flagobj2.poserr_my2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRXY_WORLD", "Covariance of position X-WORLD/Y-WORLD",
	&flagobj2.poserr_mxyw, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.covariance;pos.errorEllipse", "deg2"},

  {"ERRCXX_IMAGE", "Cxx error ellipse parameter",
	&flagobj2.poserr_cxx, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCYY_IMAGE", "Cyy error ellipse parameter",
	&flagobj2.poserr_cyy, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXY_IMAGE", "Cxy error ellipse parameter",
	&flagobj2.poserr_cxy, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXX_WORLD", "Cxx error ellipse parameter (WORLD units)",
	&flagobj2.poserr_cxxw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCYY_WORLD", "Cyy error ellipse parameter (WORLD units)",
	&flagobj2.poserr_cyyw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCXY_WORLD", "Cxy error ellipse parameter (WORLD units)",
	&flagobj2.poserr_cxyw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},

  {"ERRA_IMAGE", "RMS position error along major axis",
	&flagobj2.poserr_a, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.max;pos.errorEllipse;instr.det;meta.main", "pix"},
  {"ERRB_IMAGE", "RMS position error along minor axis",
	&flagobj2.poserr_b, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.min;pos.errorEllipse;instr.det;meta.main", "pix"},
  {"ERRTHETA_IMAGE", "Error ellipse position angle (CCW/x)",
	&flagobj2.poserr_theta, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;instr.det;meta.main", "deg"},
  {"ERRA_WORLD", "World RMS position error along major axis",
	&flagobj2.poserr_aw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"ERRB_WORLD", "World RMS position error along minor axis",
	&flagobj2.poserr_bw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"ERRTHETA_WORLD", "Error ellipse pos. angle (CCW/world-x)",
	&flagobj2.poserr_thetaw, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETA_SKY", "Native error ellipse pos. angle (east of north)",
	&flagobj2.poserr_thetas, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETA_J2000", "J2000 error ellipse pos. angle (east of north)",
	&flagobj2.poserr_theta2000, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;meta.main", "deg"},
  {"ERRTHETA_B1950", "B1950 error ellipse pos. angle (east of north)",
	&flagobj2.poserr_theta1950, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},

  {"XWIN_IMAGE", "Windowed position estimate along x",
	&flagobj2.winpos_x, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.x;instr.det", "pix"},
  {"YWIN_IMAGE", "Windowed position estimate along y",
	&flagobj2.winpos_y, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.y;instr.det", "pix"},

  {"XWIN_FOCAL", "Windowed position along focal-plane x axis",
	&flagobj2.winpos_xf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.x", ""},
  {"YWIN_FOCAL", "Windowed position along focal-plane y axis",
	&flagobj2.winpos_yf, H_FLOAT, T_DOUBLE, "%18.10e", "",
	"pos.cartesian.y", ""},

  {"XWIN_WORLD", "Windowed position along world x axis",
	&flagobj2.winpos_xw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.ra", "deg"},
  {"YWIN_WORLD", "Windowed position along world y axis",
	&flagobj2.winpos_yw, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAWIN_SKY", "Windowed right ascension  (native)",
	&flagobj2.winpos_alphas, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAWIN_SKY", "Windowed declination (native)",
	&flagobj2.winpos_deltas, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAWIN_J2000", "Windowed right ascension (J2000)",
	&flagobj2.winpos_alpha2000, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAWIN_J2000", "windowed declination (J2000)",
	&flagobj2.winpos_delta2000, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"ALPHAWIN_B1950", "Windowed right ascension (B1950)",
	&flagobj2.winpos_alpha1950, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra", "deg"},
  {"DELTAWIN_B1950", "Windowed declination (B1950)",
	&flagobj2.winpos_delta1950, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec", "deg"},

  {"X2WIN_IMAGE", "Windowed variance along x",
	&flagobj2.win_mx2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
 	"src.impactParam;instr.det", "pix2"},
  {"Y2WIN_IMAGE", "Windowed variance along y",
	&flagobj2.win_my2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;instr.det", "pix2"},
  {"XYWIN_IMAGE", "Windowed covariance between x and y",
	&flagobj2.win_mxy, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;instr.det", "pix2"},
  {"X2WIN_WORLD", "Windowed variance along X-WORLD (alpha)",
	&flagobj2.win_mx2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},
  {"Y2WIN_WORLD", "Windowed variance along Y-WORLD (delta)",
	&flagobj2.win_my2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},
  {"XYWIN_WORLD", "Windowed covariance between X-WORLD and Y-WORLD",
	&flagobj2.win_mxyw, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"src.impactParam", "deg2"},

  {"CXXWIN_IMAGE", "Windowed Cxx object ellipse parameter",
	&flagobj2.win_cxx, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CYYWIN_IMAGE", "Windowed Cyy object ellipse parameter",
	&flagobj2.win_cyy, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CXYWIN_IMAGE", "Windowed Cxy object ellipse parameter",
	&flagobj2.win_cxy, H_EXPO, T_FLOAT, "%15.7e", "pixel**(-2)",
	"src.impactParam;instr.det", "pix-2"},
  {"CXXWIN_WORLD", "Windowed Cxx object ellipse parameter (WORLD units)",
	&flagobj2.win_cxxw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},
  {"CYYWIN_WORLD", "Windowed Cyy object ellipse parameter (WORLD units)",
	&flagobj2.win_cyyw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},
  {"CXYWIN_WORLD", "Windowed Cxy object ellipse parameter (WORLD units)",
	&flagobj2.win_cxyw, H_EXPO, T_FLOAT, "%15.7e", "deg**(-2)",
	"src.impactParam", "deg-2"},

  {"AWIN_IMAGE", "Windowed profile RMS along major axis",
	&flagobj2.win_a, H_FLOAT, T_FLOAT, "%9.3f", "pixel",
	"phys.size.smajAxis;instr.det", "pix"},
  {"BWIN_IMAGE", "Windowed profile RMS along minor axis",
	&flagobj2.win_b, H_FLOAT, T_FLOAT, "%9.3f", "pixel",
	"phys.size.sminAxis;instr.det", "pix"},
  {"THETAWIN_IMAGE", "Windowed position angle (CCW/x)",
	&flagobj2.win_theta, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;instr.det", "deg"},
  {"AWIN_WORLD", "Windowed profile RMS along major axis (world units)",
	&flagobj2.win_aw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.size.smajAxis", "deg"},
  {"BWIN_WORLD", "Windowed profile RMS along minor axis (world units)",
	&flagobj2.win_bw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.size.sminAxis", "deg"},
  {"THETAWIN_WORLD", "Windowed position angle (CCW/world-x)",
	&flagobj2.win_thetaw, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng", "deg"},
  {"THETAWIN_SKY", "Windowed position angle (east of north) (native)",
	&flagobj2.win_thetas, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng", "deg"},
  {"THETAWIN_J2000", "Windowed position angle (east of north) (J2000)",
	&flagobj2.win_theta2000, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng", "deg"},
  {"THETAWIN_B1950", "Windowed position angle (east of north) (B1950)",
	&flagobj2.win_theta1950, H_FLOAT, T_FLOAT, "%+6.2f", "deg",
	"pos.posAng", "deg"},

  {"ERRX2WIN_IMAGE", "Variance of windowed pos along x",
	&flagobj2.winposerr_mx2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRY2WIN_IMAGE", "Variance of windowed pos along y",
	&flagobj2.winposerr_my2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRXYWIN_IMAGE", "Covariance of windowed pos between x and y",
	&flagobj2.winposerr_mxy, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.covariance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRX2WIN_WORLD", "Variance of windowed pos along X-WORLD (alpha)",
	&flagobj2.winposerr_mx2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRY2WIN_WORLD", "Variance of windowed pos along Y-WORLD (delta)",
	&flagobj2.winposerr_my2w, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRXYWIN_WORLD", "Covariance of windowed pos X-WORLD/Y-WORLD",
	&flagobj2.winposerr_mxyw, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.covariance;pos.errorEllipse", "deg2"},

  {"ERRCXXWIN_IMAGE", "Cxx windowed error ellipse parameter",
	&flagobj2.winposerr_cxx, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCYYWIN_IMAGE", "Cyy windowed error ellipse parameter",
	&flagobj2.winposerr_cyy, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXYWIN_IMAGE", "Cxy windowed error ellipse parameter",
	&flagobj2.winposerr_cxy, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXXWIN_WORLD", "Cxx windowed error ellipse parameter (WORLD units)",
	&flagobj2.winposerr_cxxw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCYYWIN_WORLD", "Cyy windowed error ellipse parameter (WORLD units)",
	&flagobj2.winposerr_cyyw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCXYWIN_WORLD", "Cxy windowed error ellipse parameter (WORLD units)",
	&flagobj2.winposerr_cxyw, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},

  {"ERRAWIN_IMAGE", "RMS windowed pos error along major axis",
	&flagobj2.winposerr_a, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.max;pos.errorEllipse;instr.det", "pix"},
  {"ERRBWIN_IMAGE", "RMS windowed pos error along minor axis",
	&flagobj2.winposerr_b, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.min;pos.errorEllipse;instr.det", "pix"},
  {"ERRTHETAWIN_IMAGE", "Windowed error ellipse pos angle (CCW/x)",
	&flagobj2.winposerr_theta, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;instr.det", "deg"},
  {"ERRAWIN_WORLD", "World RMS windowed pos error along major axis",
	&flagobj2.winposerr_aw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.max;pos.errorEllipse", "deg"},
  {"ERRBWIN_WORLD", "World RMS windowed pos error along minor axis",
	&flagobj2.winposerr_bw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.min;pos.errorEllipse", "deg"},
  {"ERRTHETAWIN_WORLD", "Windowed error ellipse pos. angle (CCW/world-x)",
	&flagobj2.winposerr_thetaw, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAWIN_SKY",
	"Native windowed error ellipse pos. angle (east of north)",
	&flagobj2.winposerr_thetas, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAWIN_J2000",
	"J2000 windowed error ellipse pos. angle (east of north)",
	&flagobj2.winposerr_theta2000, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAWIN_B1950",
	"B1950 windowed error ellipse pos. angle (east of north)",
	&flagobj2.winposerr_theta1950, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},

  {"NITER_WIN", "Number of iterations for WIN centering",
	&flagobj2.winpos_niter, H_INT, T_SHORT, "%3d", "",
	"meta.number", ""},

  {"MU_THRESHOLD", "Detection threshold above background",
	&flagobj2.threshmu, H_FLOAT, T_FLOAT, "%8.4f", "mag * arcsec**(-2)",
	"instr.sensitivity;phot.mag.sb", "mag.arcsec-2"},
  {"MU_MAX", "Peak surface brightness above background",
	&flagobj2.maxmu, H_FLOAT, T_FLOAT, "%8.4f", "mag * arcsec**(-2)",
	"phot.mag.sb;stat.max", "mag.arcsec-2"},
  {"ISOAREA_WORLD", "Isophotal area above Analysis threshold",
	&flagobj2.npixw, H_FLOAT, T_FLOAT, "%12.7g", "deg**2",
	"phys.angArea", "deg2"},
  {"ISOAREAF_WORLD", "Isophotal area (filtered) above Detection threshold",
	&flagobj2.fdnpixw, H_FLOAT, T_FLOAT, "%12.7g", "deg**2",
	"phys.angArea", "deg2"},
  {"ISO0", "Isophotal area at level 0",
	&flagobj.iso[0], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO1", "Isophotal area at level 1",
	&flagobj.iso[1], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO2", "Isophotal area at level 2",
	&flagobj.iso[2], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO3", "Isophotal area at level 3",
	&flagobj.iso[3], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO4", "Isophotal area at level 4",
	&flagobj.iso[4], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO5", "Isophotal area at level 5",
	&flagobj.iso[5], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO6", "Isophotal area at level 6",
	&flagobj.iso[6], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},
  {"ISO7", "Isophotal area at level 7",
	&flagobj.iso[7], H_INT, T_LONG, "%8d", "pixel**2",
	"phys.area;instr.det", "pix2"},

  {"FLAGS", "Extraction flags",
	&flagobj.flag, H_INT, T_SHORT, "%3d", "",
	"meta.code.qual", ""},
  {"FLAGS_WEIGHT", "Weighted extraction flags",
	&flagobj.wflag, H_INT, T_SHORT, "%1d", "",
	"meta.code.qual", ""},
   {"FLAGS_WIN", "Flags for WINdowed parameters",
	&flagobj2.win_flag, H_INT, T_SHORT, "%3d", "",
	"meta.code.qual", ""},
   {"IMAFLAGS_ISO", "FLAG-image flags OR'ed over the iso. profile",
	flagobj.imaflag, H_INT, T_LONG, "%9u", "",
	"meta.code.qual", "", 1, &prefs.imaflag_size},
  {"NIMAFLAGS_ISO", "Number of flagged pixels entering IMAFLAGS_ISO",
	flagobj.imanflag, H_INT, T_LONG, "%9d", "",
	"meta.number", "", 1, &prefs.imanflag_size},
  {"NLOWWEIGHT_ISO", "Nb of pixels with low weight over the iso. profile",
	&flagobj.nzwpix, H_INT, T_LONG, "%9d", "",
	"meta.number", ""},
  {"NLOWDWEIGHT_ISO", "Nb of pixels with low det. weight over the iso. profile",
	&flagobj.nzdwpix, H_INT, T_LONG, "%9d", "",
	"meta.number", ""},

  {"FWHM_IMAGE", "FWHM assuming a gaussian core",
	&flagobj.fwhm, H_FLOAT, T_FLOAT, "%8.2f", "pixel",
	"phys.size.diameter;instr.det.psf", "pix"},
  {"FWHM_WORLD", "FWHM assuming a gaussian core",
	&flagobj2.fwhmw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.angSize;instr.det.psf", "deg"},
  {"ELONGATION", "A_IMAGE/B_IMAGE",
	&flagobj2.elong, H_FLOAT, T_FLOAT, "%8.3f", "",
	"src.ellipticity;arith.ratio;instr.det", ""},
  {"ELLIPTICITY", "1 - B_IMAGE/A_IMAGE",
	&flagobj2.ellip, H_FLOAT, T_FLOAT, "%8.3f", "",
	"src.ellipticity;instr.det	", ""},
  {"POLAR_IMAGE", "(A_IMAGE^2 - B_IMAGE^2)/(A_IMAGE^2 + B_IMAGE^2)",
	&flagobj2.polar, H_FLOAT, T_FLOAT, "%7.5f", "",
	"src.ellipticity;instr.det", ""},
  {"POLAR_WORLD", "(A_WORLD^2 - B_WORLD^2)/(A_WORLD^2 + B_WORLD^2)",
	&flagobj2.polarw, H_FLOAT, T_FLOAT, "%7.5f", "",
	"src.ellipticity", ""},
  {"POLARWIN_IMAGE", "(AWIN^2 - BWIN^2)/(AWIN^2 + BWIN^2)",
	&flagobj2.win_polar, H_FLOAT, T_FLOAT, "%7.5f", "",
	"src.ellipticity;instr.det", ""},
  {"POLARWIN_WORLD", "(AWIN^2 - BWIN^2)/(AWIN^2 + BWIN^2)",
	&flagobj2.win_polarw, H_FLOAT, T_FLOAT, "%7.5f", "",
	"src.ellipticity", ""},
  {"CLASS_STAR", "S/G classifier output",
	&flagobj2.sprob, H_FLOAT, T_FLOAT, "%6.3f", "",
	"src.class.starGalaxy", ""},
  {"VIGNET", "Pixel data around detection",
	&flagobj2.vignet, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"obs.image", "ct", 2, prefs.vignetsize},
  {"VIGNET_SHIFT", "Pixel data around detection, corrected for shift",
	&flagobj2.vigshift, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"obs.image", "ct", 2, prefs.vigshiftsize},
  {"VECTOR_ASSOC", "ASSOCiated parameter vector",
	&flagobj2.assoc, H_FLOAT, T_DOUBLE, "%12.7g", "",
	"src", "", 1, &prefs.assoc_size},
  {"NUMBER_ASSOC", "Number of ASSOCiated IDs",
	&flagobj2.assoc_number, H_INT, T_LONG, "%10d", "",
	"meta.number;src", ""},

  {"THRESHOLDMAX", "Maximum threshold possible for detection",
	&flagobj.dthresh, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.flux.sb;stat.max", "ct"},

  {"FLUX_GROWTH", "Cumulated growth-curve",
	&flagobj2.flux_growth, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count", "ct", 1, &prefs.flux_growthsize},
  {"FLUX_GROWTHSTEP", "Step for growth-curves",
	&flagobj2.flux_growthstep, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.distance", "pix"},
  {"MAG_GROWTH", "Cumulated magnitude growth-curve",
	&flagobj2.mag_growth, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag", 1, &prefs.mag_growthsize},
  {"MAG_GROWTHSTEP", "Step for growth-curves",
	&flagobj2.mag_growthstep, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.distance", "pix"},
  {"FLUX_RADIUS", "Fraction-of-light radii",
	&flagobj2.flux_radius, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"phys.size.radius;instr.det", "pix",  1, &prefs.flux_radiussize},

  {"FWHMPSF_IMAGE", "FWHM of the local PSF model",
	&flagobj2.fwhm_psf, H_FLOAT, T_FLOAT, "%8.3f", "pixel",
	"phys.size.diameter;instr.det.psf", "pix"},
  {"FWHMPSF_WORLD", "FWHM of the local PSF model (world units)",
	&flagobj2.fwhmw_psf, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"phys.angSize;instr.det.psf", "deg"},

  {"XPSF_IMAGE", "X coordinate from PSF-fitting",
	&flagobj2.x_psf, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.x;stat.fit.param;instr.det", "pix"},
  {"YPSF_IMAGE", "Y coordinate from PSF-fitting",
	&flagobj2.y_psf, H_FLOAT, T_DOUBLE, "%11.4f", "pixel",
	"pos.cartesian.y;stat.fit.param;instr.det", "pix"},
  {"XPSF_WORLD", "PSF position along world x axis",
	&flagobj2.xw_psf, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"YPSF_WORLD", "PSF position along world y axis",
	&flagobj2.yw_psf, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAPSF_SKY", "Right ascension of the fitted PSF (native)",
	&flagobj2.alphas_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAPSF_SKY", "Declination of the fitted PSF (native)",
	&flagobj2.deltas_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAPSF_J2000", "Right ascension of the fitted PSF (J2000)",
	&flagobj2.alpha2000_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAPSF_J2000", "Declination of the fitted PSF (J2000)",
	&flagobj2.delta2000_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAPSF_B1950", "Right ascension of the fitted PSF (B1950)",
	&flagobj2.alpha1950_psf, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAPSF_B1950", "Declination of the fitted PSF (B1950)",
	&flagobj2.delta1950_psf, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"FLUX_PSF", "Flux from PSF-fitting",
	&flagobj2.flux_psf, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_PSF", "RMS flux error for PSF-fitting",
	&flagobj2.fluxerr_psf, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.stdev;phot.count", "ct"},
  {"MAG_PSF", "Magnitude from PSF-fitting",
	&flagobj2.mag_psf, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag", "mag"},
  {"MAGERR_PSF", "RMS magnitude error from PSF-fitting",
	&flagobj2.magerr_psf, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.stdev;phot.mag", "mag"},

  {"NITER_PSF", "Number of iterations for PSF-fitting",
	&flagobj2.niter_psf, H_INT, T_SHORT, "%3d", "",
	"meta.number", ""},
  {"CHI2_PSF", "Reduced chi2 from PSF-fitting",
	&flagobj2.chi2_psf, H_FLOAT, T_FLOAT, "%9.4g", "",
	"stat.fit.chi2", ""},

  {"ERRX2PSF_IMAGE", "Variance of PSF position along x",
	&flagobj2.poserrmx2_psf, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRY2PSF_IMAGE", "Variance of PSF position along y",
	&flagobj2.poserrmy2_psf, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRXYPSF_IMAGE", "Covariance of PSF position between x and y",
	&flagobj2.poserrmxy_psf, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.covariance;pos.errorEllipse;instr.det", "pix2"},
  {"ERRX2PSF_WORLD", "Variance of PSF position along X-WORLD (alpha)",
	&flagobj2.poserrmx2w_psf, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRY2PSF_WORLD", "Variance of PSF position along Y-WORLD (delta)",
	&flagobj2.poserrmy2w_psf, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse", "deg2"},
  {"ERRXYPSF_WORLD", "Covariance of PSF position X-WORLD/Y-WORLD",
	&flagobj2.poserrmxyw_psf, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.covariance;pos.errorEllipse", "deg2"},

  {"ERRCXXPSF_IMAGE", "Cxx PSF error ellipse parameter",
	&flagobj2.poserrcxx_psf, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCYYPSF_IMAGE", "Cyy PSF error ellipse parameter",
	&flagobj2.poserrcyy_psf, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXYPSF_IMAGE", "Cxy PSF error ellipse parameter",
	&flagobj2.poserrcxy_psf, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;instr.det", "pix-2"},
  {"ERRCXXPSF_WORLD", "Cxx PSF error ellipse parameter (WORLD units)",
	&flagobj2.poserrcxxw_psf, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCYYPSF_WORLD", "Cyy PSF error ellipse parameter (WORLD units)",
	&flagobj2.poserrcyyw_psf, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},
  {"ERRCXYPSF_WORLD", "Cxy PSF error ellipse parameter (WORLD units)",
	&flagobj2.poserrcxyw_psf, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse", "deg-2"},

  {"ERRAPSF_IMAGE", "PSF RMS position error along major axis",
	&flagobj2.poserra_psf, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.max;pos.errorEllipse;instr.det", "pix"},
  {"ERRBPSF_IMAGE", "PSF RMS position error along minor axis",
	&flagobj2.poserrb_psf, H_FLOAT, T_FLOAT, "%9.5f", "pixel",
	"stat.stdev;stat.min;pos.errorEllipse;instr.det", "pix"},
  {"ERRTHETAPSF_IMAGE", "PSF error ellipse position angle (CCW/x)",
	&flagobj2.poserrtheta_psf, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;instr.det", "deg"},
  {"ERRAPSF_WORLD", "World PSF RMS position error along major axis",
	&flagobj2.poserraw_psf, H_FLOAT, T_FLOAT, "%12.7g", "pixel",
	"stat.stdev;stat.max;pos.errorEllipse", "deg"},
  {"ERRBPSF_WORLD", "World PSF RMS position error along minor axis",
	&flagobj2.poserrbw_psf, H_FLOAT, T_FLOAT, "%12.7g", "pixel",
	"stat.stdev;stat.min;pos.errorEllipse", "deg"},
  {"ERRTHETAPSF_WORLD", "PSF error ellipse pos. angle (CCW/world-x)",
	&flagobj2.poserrthetaw_psf, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAPSF_SKY", "Native PSF error ellipse pos. angle (east of north)",
	&flagobj2.poserrthetas_psf, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAPSF_J2000", "J2000 PSF error ellipse pos. angle (east of north)",
	&flagobj2.poserrtheta2000_psf, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},
  {"ERRTHETAPSF_B1950", "B1950 PSF error ellipse pos. angle (east of north)",
	&flagobj2.poserrtheta1950_psf, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse", "deg"},

  {"DURATION_ANALYSIS", "Duration of analysis for this source",
	&flagobj2.analtime, H_FLOAT, T_FLOAT, "%9.4g", "s",
	"time.duration;time.processing", "s"},

#ifdef USE_MODEL
#include "paramprofit.h"
#endif

  {""}
  };

