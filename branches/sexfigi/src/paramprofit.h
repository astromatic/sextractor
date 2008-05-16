/*
 				paramprofit.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Profile-fitting parameter list for catalog data.
*
*	Last modify:	10/05/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/



  {"VECTOR_PROF", "Profile-fitting coefficients",
	&outobj2.prof_vector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.fit.param;src.morph.param", "", 1, &prefs.prof_vectorsize},
  {"VECTOR_PROFERR", "Profile-fitting coefficient uncertainties",
	&outobj2.prof_errvector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.stdev;stat.fit;src.morph.param", "", 1,
	&prefs.prof_errvectorsize},
  {"CHI2_PROF", "Reduced Chi2 of the fit",
	&outobj2.prof_chi2, H_FLOAT, T_FLOAT, "%12.7g", "",
	"stat.fit.chi2;src.morph", ""},
  {"FLAGS_PROF", "Profile-fitting flags",
	&outobj2.prof_flag, H_INT, T_BYTE, "%3d", "",
	"meta.code;stat.fit;src.morph", ""},
  {"NITER_PROF", "Number of iterations for profile-fitting",
	&outobj2.prof_niter, H_INT, T_SHORT, "%3d", "",
	"meta.number;stat.fit;src.morph", ""},
  {"FLUX_PROF", "Flux from profile-fitting",
	&outobj2.flux_prof, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_PROF", "RMS error on profile-fitting flux",
	&outobj2.fluxerr_prof, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_PROF", "Magnitude from profile-fitting",
	&outobj2.mag_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_PROF", "RMS error on profile-fitting magnitude",
	&outobj2.mag_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},
  {"XPROF_IMAGE", "X coordinate from profile-fitting",
	&outobj2.x_prof, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.cartesian.x;stat.fit.param;instr.det;meta.main", "pix"},
  {"YPROF_IMAGE", "Y coordinate from profile-fitting",
	&outobj2.y_prof, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.cartesian.y;stat.fit.param;instr.det;meta.main", "pix"},
  {"X2PROF_IMAGE", "Variance along x from profile-fitting",
	&outobj2.prof_mx2, H_EXPO, T_DOUBLE, "%15.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"Y2PROF_IMAGE", "Variance along y from profile-fitting",
	&outobj2.prof_my2, H_EXPO, T_DOUBLE, "%15.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"XYPROF_IMAGE", "Covariance between x and y from profile-fitting",
	&outobj2.prof_mxy, H_EXPO, T_DOUBLE, "%15.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"E1PROF_IMAGE", "Ellipticity component from profile-fitting",
	&outobj2.prof_e1, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"E2PROF_IMAGE", "Ellipticity component from profile-fitting",
	&outobj2.prof_e2, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"EPS1PROF_IMAGE", "Ellipticity component (quadratic) from profile-fitting",
	&outobj2.prof_eps1, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"EPS2PROF_IMAGE", "Ellipticity component (quadratic) from profile-fitting",
	&outobj2.prof_eps2, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},

  {"FLUX_BACKOFFSET", "Background offset from fitting",
	&outobj2.prof_offset_flux, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"instr.skyLevel;arith.diff;stat.fit.param", "ct"},
  {"FLUXERR_BACKOFFSET", "RMS error on fitted background offset",
	&outobj2.prof_offset_fluxerr, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.error;instr.skyLevel;arith.diff;stat.fit.param", "ct"},

  {"FLUX_SPHEROID", "Spheroid total flux from fitting",
	&outobj2.prof_spheroid_flux, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_SPHEROID", "RMS error on fitted spheroid total flux",
	&outobj2.prof_spheroid_fluxerr, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_SPHEROID", "Spheroid total magnitude from fitting",
	&outobj2.prof_spheroid_mag, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_SPHEROID", "RMS error on fitted spheroid total magnitude",
	&outobj2.prof_spheroid_magerr, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},
  {"SPHEROID_REFF_IMAGE", "Spheroid effective radius from fitting",
	&outobj2.prof_spheroid_reff, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"SPHEROID_REFFERR_IMAGE", "RMS error on fitted spheroid effective radius",
	&outobj2.prof_spheroid_refferr, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"stat.error;src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"SPHEROID_REFF_WORLD", "Spheroid effective radius from fitting",
	&outobj2.prof_spheroid_reffw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"SPHEROID_REFFERR_WORLD", "RMS error on fitted spheroid effective radius",
	&outobj2.prof_spheroid_refferrw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.error;src.morph.scLength;stat.fit.param", "deg"},
  {"SPHEROID_ASPECT_IMAGE", "Spheroid aspect ratio from fitting",
	&outobj2.prof_spheroid_aspect, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"SPHEROID_ASPECTERR_IMAGE", "RMS error on fitted spheroid aspect ratio",
	&outobj2.prof_spheroid_aspecterr, H_FLOAT, T_FLOAT, "%6.4f", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"SPHEROID_ASPECT_WORLD", "Spheroid aspect ratio from fitting",
	&outobj2.prof_spheroid_aspectw, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"SPHEROID_ASPECTERR_WORLD", "RMS error on fitted spheroid aspect ratio",
	&outobj2.prof_spheroid_aspecterrw, H_FLOAT, T_FLOAT, "%6.4f", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"SPHEROID_THETA_IMAGE", "Spheroid position angle (CCW/x) from fitting",
	&outobj2.prof_spheroid_theta, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"SPHEROID_THETAERR_IMAGE", "RMS error on spheroid position angle",
	&outobj2.prof_spheroid_thetaerr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"SPHEROID_THETA_WORLD", "Spheroid position angle (CCW/world-x)",
	&outobj2.prof_spheroid_thetaw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"SPHEROID_THETAERR_WORLD", "RMS error on spheroid position angle",
	&outobj2.prof_spheroid_thetaerrw, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param", "deg"},
  {"SPHEROID_THETA_SKY", "Spheroid position angle (east of north, native)",
	&outobj2.prof_spheroid_thetas, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"SPHEROID_THETA_J2000", "Spheroid position angle (east of north, J2000)",
	&outobj2.prof_spheroid_theta2000, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"SPHEROID_THETA_B1950", "Spheroid position angle (east of north, B1950)",
	&outobj2.prof_spheroid_theta1950, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"SPHEROID_SERSICN", "Spheroid Sersic index from fitting",
	&outobj2.prof_spheroid_sersicn, H_FLOAT, T_FLOAT, "%6.3f", "",
	"src.morph;stat.fit.param", ""},
  {"SPHEROID_SERSICNERR", "RMS error on fitted spheroid Sersic index",
	&outobj2.prof_spheroid_sersicnerr, H_FLOAT, T_FLOAT, "%6.3f", "",
	"stat.error;src.morph;stat.fit.param", ""},

  {"FLUX_DISK", "Disk total flux from fitting",
	&outobj2.prof_disk_flux, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_DISK", "RMS error on fitted disk total flux",
	&outobj2.prof_disk_fluxerr, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_DISK", "Disk total magnitude from fitting",
	&outobj2.prof_disk_mag, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_DISK", "RMS error on fitted disk total magnitude",
	&outobj2.prof_disk_magerr, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},
  {"DISK_SCALE_IMAGE", "Disk scalelength from fitting",
	&outobj2.prof_disk_scale, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"DISK_SCALEERR_IMAGE", "RMS error on fitted disk scalelength",
	&outobj2.prof_disk_scaleerr, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"stat.error;src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"DISK_SCALE_WORLD", "Disk scalelength from fitting",
	&outobj2.prof_disk_scalew, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"DISK_SCALEERR_WORLD", "RMS error on fitted disk scalelength",
	&outobj2.prof_disk_scaleerrw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.error;src.morph.scLength;stat.fit.param", "deg"},
  {"DISK_ASPECT_IMAGE", "Disk aspect ratio from fitting",
	&outobj2.prof_disk_aspect, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"DISK_ASPECTERR_IMAGE", "RMS error on fitted disk aspect ratio",
	&outobj2.prof_disk_aspecterr, H_FLOAT, T_FLOAT, "%6.4f", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"DISK_ASPECT_WORLD", "Disk aspect ratio from fitting",
	&outobj2.prof_disk_aspectw, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"DISK_ASPECTERR_WORLD", "RMS error on disk aspect ratio",
	&outobj2.prof_disk_aspecterrw, H_FLOAT, T_FLOAT, "%6.4f", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"DISK_INCLINATION_IMAGE", "Disk inclination from fitting",
	&outobj2.prof_disk_inclination, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"src.morph;stat.fit.param;instr.det", "deg"},
  {"DISK_INCLINATION_WORLD", "Disk inclination from fitting",
	&outobj2.prof_disk_inclinationw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"src.morph;stat.fit.param", "deg"},
  {"DISK_THETA_IMAGE", "Disk position angle (CCW/x) from fitting",
	&outobj2.prof_disk_theta, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"DISK_THETAERR_IMAGE", "RMS error on fitted disk position angle",
	&outobj2.prof_disk_thetaerr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"DISK_THETA_WORLD", "Disk position angle (CCW/world-x)",
	&outobj2.prof_disk_thetaw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"DISK_THETAERR_WORLD", "RMS error on disk position angle",
	&outobj2.prof_disk_thetaerrw, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param", "deg"},
  {"DISK_THETA_SKY", "Disk position angle (east of north, native)",
	&outobj2.prof_disk_thetas, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"DISK_THETA_J2000", "Disk position angle (east of north, J2000)",
	&outobj2.prof_disk_theta2000, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"DISK_THETA_B1950", "Disk position angle (east of north, B1950)",
	&outobj2.prof_disk_theta1950, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"FLUX_BAR", "Bar total flux from fitting",
	&outobj2.prof_bar_flux, H_FLOAT, T_FLOAT, "%12.g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_BAR", "RMS error on fitted total bar flux",
	&outobj2.prof_bar_fluxerr, H_FLOAT, T_FLOAT, "%12.g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_BAR", "Bar total magnitude from fitting",
	&outobj2.prof_bar_mag, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_BAR", "RMS error on fitted total bar magnitude",
	&outobj2.prof_bar_magerr, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},
  {"BAR_LENGTH_IMAGE", "Bar length from fitting",
	&outobj2.prof_bar_length, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"BAR_LENGTH_WORLD", "Bar length from fitting",
	&outobj2.prof_bar_lengthw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"BAR_ASPECT_IMAGE", "Bar aspect ratio from fitting",
	&outobj2.prof_bar_aspect, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"BAR_ASPECT_WORLD", "Bar aspect ratio from fitting",
	&outobj2.prof_bar_aspectw, H_FLOAT, T_FLOAT, "%12.7g", "",
	"phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"BAR_POSANG", "Bar true position angle (CCW/disk maj.axis) from fitting",
	&outobj2.prof_bar_posang, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.bodyrc.long;src.morph;stat.fit.param", "deg"},
  {"BAR_THETA_IMAGE", "Bar projected angle (CCW/x) from fitting",
	&outobj2.prof_bar_theta, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"BAR_THETA_WORLD", "Bar projected angle (CCW/world-x) from fitting",
	&outobj2.prof_bar_thetaw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"BAR_THETA_SKY", "Bar projected angle (east of north, native) from fitting",
	&outobj2.prof_bar_thetas, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"BAR_THETA_J2000", "Bar projected angle (east of north, J2000) from fitting",
	&outobj2.prof_bar_theta2000, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"BAR_THETA_B1950", "Bar projected angle (east of north, B1950) from fitting",
	&outobj2.prof_bar_theta1950, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"FLUX_ARMS", "Total flux in spiral arms from fitting",
	&outobj2.prof_arms_flux, H_FLOAT, T_FLOAT, "%12.g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"MAG_ARMS", "Total magnitude in spiral arms from fitting",
	&outobj2.prof_arms_mag, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"ARMS_SCALE_IMAGE", "Spiral arms scale length from fitting",
	&outobj2.prof_arms_scale, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"ARMS_SCALE_WORLD", "Spiral arms scale length from fitting",
	&outobj2.prof_arms_scalew, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"ARMS_POSANG", "Pos. angle (CCW/disk maj.axis) of spiral arms from fitting",
	&outobj2.prof_arms_posang, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.bodyrc.long;src.morph;stat.fit.param", "deg"},
/*
  {"ARMS_THETA_WORLD", "Pos. angle (CCW/world-x) of spiral arms",
	&outobj2.prof_arms_thetaw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"ARMS_THETA_SKY", "Pos. angle (east of north, native) of spiral arms",
	&outobj2.prof_arms_thetas, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"ARMS_THETA_J2000", "Pos. angle (east of north, J2000) of spiral arms",
	&outobj2.prof_arms_theta2000, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"ARMS_THETA_B1950", "Pos. angle (east of north, B1950) of spiral arms",
	&outobj2.prof_arms_theta1950, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
*/
  {"ARMS_PITCH", "Pitch angle of spiral arms from fitting",
	&outobj2.prof_arms_pitch, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"ARMS_START_IMAGE", "Starting radius of spiral arms from fitting",
	&outobj2.prof_arms_start, H_FLOAT, T_FLOAT, "%6.4f", "pixel",
	"pos.distance;src.morph;stat.fit.param;instr.det", "pix"},
  {"ARMS_START_WORLD", "Starting radius of spiral arms from fitting",
	&outobj2.prof_arms_startw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"pos.distance;src.morph;stat.fit.param", "deg"},
  {"ARMS_QUADFRAC", "Fraction of spiral arms in quadrature from fitting",
	&outobj2.prof_arms_quadfrac, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"phot.count;arith.ratio;src.morph;stat.fit.param", "deg"},
