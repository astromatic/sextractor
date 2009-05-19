/*
 				paramprofit.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Model-fitting parameter list for catalog data.
*
*	Last modify:	19/05/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  {"VECTOR_MODEL", "Model-fitting coefficients",
	&outobj2.prof_vector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.fit.param;src.morph.param", "", 1, &prefs.prof_vectorsize},
  {"VECTOR_MODELERR", "Model-fitting coefficient uncertainties",
	&outobj2.prof_errvector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.stdev;stat.fit;src.morph.param", "", 1,
	&prefs.prof_errvectorsize},
  {"CHI2_MODEL", "Reduced Chi2 of the fit",
	&outobj2.prof_chi2, H_FLOAT, T_FLOAT, "%12.7g", "",
	"stat.fit.chi2;src.morph", ""},
  {"FLAGS_MODEL", "Model-fitting flags",
	&outobj2.prof_flag, H_INT, T_BYTE, "%3d", "",
	"meta.code;stat.fit;src.morph", ""},
  {"NITER_MODEL", "Number of iterations for model-fitting",
	&outobj2.prof_niter, H_INT, T_SHORT, "%3d", "",
	"meta.number;stat.fit;src.morph", ""},
  {"FLUX_MODEL", "Flux from model-fitting",
	&outobj2.flux_prof, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"phot.count;stat.fit.param", "ct"},
  {"FLUXERR_MODEL", "RMS error on model-fitting flux",
	&outobj2.fluxerr_prof, H_FLOAT, T_FLOAT, "%12.7g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_MODEL", "Magnitude from model-fitting",
	&outobj2.mag_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_MODEL", "RMS error on model-fitting magnitude",
	&outobj2.mag_prof, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},

  {"XMODEL_IMAGE", "X coordinate from model-fitting",
	&outobj2.x_prof, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.cartesian.x;stat.fit.param;instr.det;meta.main", "pix"},
  {"YMODEL_IMAGE", "Y coordinate from model-fitting",
	&outobj2.y_prof, H_FLOAT, T_FLOAT, "%10.3f", "pixel",
	"pos.cartesian.y;stat.fit.param;instr.det;meta.main", "pix"},

  {"XMODEL_WORLD", "Fitted position along world x axis",
	&outobj2.xw_prof, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"YMODEL_WORLD", "Fitted position along world y axis",
	&outobj2.yw_prof, H_FLOAT, T_DOUBLE, "%18.10e", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAMODEL_SKY", "Fitted position along right ascension  (native)",
	&outobj2.alphas_prof, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAMODEL_SKY", "Fitted position along declination (native)",
	&outobj2.deltas_prof, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAMODEL_J2000", "Fitted position along right ascension (J2000)",
	&outobj2.alpha2000_prof, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAMODEL_J2000", "Fitted position along declination (J2000)",
	&outobj2.delta2000_prof, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ALPHAMODEL_B1950", "Fitted position along right ascension (B1950)",
	&outobj2.alpha1950_prof, H_FLOAT, T_DOUBLE, "%11.7f", "deg",
	"pos.eq.ra;stat.fit.param", "deg"},
  {"DELTAMODEL_B1950", "Fitted position along declination (B1950)",
	&outobj2.delta1950_prof, H_FLOAT, T_DOUBLE, "%+11.7f", "deg",
	"pos.eq.dec;stat.fit.param", "deg"},

  {"ERRX2MODEL_IMAGE", "Variance of fitted position along x",
	&outobj2.poserrmx2_prof, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;stat.fit.param;instr.det", "pix2"},
  {"ERRY2MODEL_IMAGE", "Variance of fitted position along y",
	&outobj2.poserrmy2_prof, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.variance;pos.errorEllipse;stat.fit.param;instr.det", "pix2"},
  {"ERRXYMODEL_IMAGE", "Covariance of fitted position between x and y",
	&outobj2.poserrmxy_prof, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"stat.covariance;pos.errorEllipse;stat.fit.param;instr.det", "pix2"},
  {"ERRX2MODEL_WORLD", "Variance of fitted position along X-WORLD (alpha)",
	&outobj2.poserrmx2w_prof, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse;stat.fit.param", "deg2"},
  {"ERRY2MODEL_WORLD", "Variance of fitted position along Y-WORLD (delta)",
	&outobj2.poserrmy2w_prof, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.variance;pos.errorEllipse;stat.fit.param", "deg2"},
  {"ERRXYMODEL_WORLD", "Covariance of fitted position X-WORLD/Y-WORLD",
	&outobj2.poserrmxyw_prof, H_EXPO, T_DOUBLE, "%18.10e", "deg**2",
	"stat.covariance;pos.errorEllipse;stat.fit.param", "deg2"},

  {"ERRCXXMODEL_IMAGE", "Cxx error ellipse parameter of fitted position",
	&outobj2.poserrcxx_prof, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;stat.fit.param;instr.det", "pix-2"},
  {"ERRCYYMODEL_IMAGE", "Cyy error ellipse parameter of fitted position",
	&outobj2.poserrcyy_prof, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;stat.fit.param;instr.det", "pix-2"},
  {"ERRCXYMODEL_IMAGE", "Cxy error ellipse parameter of fitted position",
	&outobj2.poserrcxy_prof, H_EXPO, T_FLOAT, "%12.7g", "pixel**(-2)",
	"src.impactParam;pos.errorEllipse;stat.fit.param;instr.det", "pix-2"},
  {"ERRCXXMODEL_WORLD", "Cxx fitted error ellipse parameter (WORLD units)",
	&outobj2.poserrcxxw_prof, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse;stat.fit.param", "deg-2"},
  {"ERRCYYMODEL_WORLD", "Cyy fitted error ellipse parameter (WORLD units)",
	&outobj2.poserrcyyw_prof, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipse;stat.fit.param", "deg-2"},
  {"ERRCXYMODEL_WORLD", "Cxy fitted error ellipse parameter (WORLD units)",
	&outobj2.poserrcxyw_prof, H_EXPO, T_FLOAT, "%12.7g", "deg**(-2)",
	"src.impactParam;pos.errorEllipsestat.fit.param", "deg-2"},

  {"ERRAMODEL_IMAGE", "RMS error of fitted position along major axis",
	&outobj2.poserra_prof, H_FLOAT, T_FLOAT, "%8.4f", "pixel",
	"stat.stdev;stat.max;pos.errorEllipse;stat.fit.param;instr.det", "pix"},
  {"ERRBMODEL_IMAGE", "RMS error of fitted position along minor axis",
	&outobj2.poserrb_prof, H_FLOAT, T_FLOAT, "%8.4f", "pixel",
	"stat.stdev;stat.min;pos.errorEllipse;stat.fit.param;instr.det", "pix"},
  {"ERRTHETAMODEL_IMAGE", "Error ellipse pos.angle of fitted position (CCW/x)",
	&outobj2.poserrtheta_prof, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;stat.fit.param;instr.det", "deg"},
  {"ERRAMODEL_WORLD", "World RMS error of fitted position along major axis",
	&outobj2.poserraw_prof, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.max;pos.errorEllipse;stat.fit.param", "deg"},
  {"ERRBMODEL_WORLD", "World RMS error of fitted position along minor axis",
	&outobj2.poserrbw_prof, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.stdev;stat.min;pos.errorEllipse;stat.fit.param", "deg"},
  {"ERRTHETAMODEL_WORLD", "Error ellipse pos.angle of fitted position (CCW/world-x)",
	&outobj2.poserrthetaw_prof, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;stat.fit.param", "deg"},
  {"ERRTHETAMODEL_SKY", "Native fitted error ellipse pos. angle (east of north)",
	&outobj2.poserrthetas_prof, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;stat.fit.param", "deg"},
  {"ERRTHETAMODEL_J2000", "J2000 fitted error ellipse pos. angle (east of north)",
	&outobj2.poserrtheta2000_prof, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;stat.fit.param", "deg"},
  {"ERRTHETAMODEL_B1950", "B1950 fitted error ellipse pos. angle (east of north)",
	&outobj2.poserrtheta1950_prof, H_FLOAT, T_FLOAT, "%6.2f", "deg",
	"pos.posAng;pos.errorEllipse;stat.fit.param", "deg"},


  {"X2MODEL_IMAGE", "Variance along x from model-fitting",
	&outobj2.prof_mx2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"Y2MODEL_IMAGE", "Variance along y from model-fitting",
	&outobj2.prof_my2, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"XYMODEL_IMAGE", "Covariance between x and y from model-fitting",
	&outobj2.prof_mxy, H_EXPO, T_DOUBLE, "%18.10e", "pixel**2",
	"src.impactParam;stat.fit;instr.det", "pix2"},
  {"E1MODEL_IMAGE", "Ellipticity component from model-fitting",
	&outobj2.prof_e1, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"E2MODEL_IMAGE", "Ellipticity component from model-fitting",
	&outobj2.prof_e2, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"EPS1MODEL_IMAGE", "Ellipticity component (quadratic) from model-fitting",
	&outobj2.prof_eps1, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},
  {"EPS2MODEL_IMAGE", "Ellipticity component (quadratic) from model-fitting",
	&outobj2.prof_eps2, H_FLOAT, T_FLOAT, "%10.6f", "",
	"src.ellipticity;stat.fit;instr.det", ""},

  {"CONCENTRATION_MODEL", "Concentration parameter from model-fitting",
	&outobj2.prof_concentration, H_FLOAT, T_FLOAT, "%7.4f", "",
	"src.morph.param", ""},
  {"CLASS_STAR_MODEL", "S/G classifier from model-fitting",
	&outobj2.prof_class_star, H_FLOAT, T_FLOAT, "%7.4f", "",
	"src.class.starGalaxy", ""},

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
  {"DISK_SCALE_WORLD", "Disk scalelength from fitting (world coords)",
	&outobj2.prof_disk_scalew, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"DISK_SCALEERR_WORLD", "RMS error on fitted disk scalelength (world coords)",
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
  {"DISK_INCLINATION", "Disk inclination from fitting",
	&outobj2.prof_disk_inclination, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"src.morph;stat.fit.param;instr.det", "deg"},
  {"DISK_INCLINATIONERR", "RMS error on disk inclination from fitting",
	&outobj2.prof_disk_inclinationerr, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"stat.error;src.morph;stat.fit.param;instr.det", "deg"},
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
  {"DISK_PATTERN_VECTOR", "Disk pattern fitted coefficients",
	&outobj2.prof_disk_patternvector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.fit.param;src.morph.param", "", 1,
	&prefs.prof_disk_patternvectorsize},
  {"DISK_PATTERNMOD_VECTOR", "Disk pattern fitted moduli",
	&outobj2.prof_disk_patternmodvector, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.fit.param;src.morph.param", "", 1,
	&prefs.prof_disk_patternmodvectorsize},
  {"DISK_PATTERNARG_VECTOR", "Disk pattern fitted arguments",
	&outobj2.prof_disk_patternargvector, H_FLOAT, T_FLOAT, "%12.4g", "deg",
	"stat.fit.param;src.morph.param", "deg", 1,
	&prefs.prof_disk_patternargvectorsize},
  {"DISK_PATTERN_SPIRAL", "Disk pattern spiral index",
	&outobj2.prof_disk_patternspiral, H_FLOAT, T_FLOAT, "%12.4g", "",
	"stat.fit.param;src.morph.param", ""},
/*
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
  {"BAR_LENGTHERR_IMAGE", "RMS error on fitted bar length",
	&outobj2.prof_bar_lengtherr, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"stat.error;src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"BAR_LENGTH_WORLD", "Bar length from fitting",
	&outobj2.prof_bar_lengthw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"BAR_LENGTHERR_WORLD", "RMS error on fitted bar length (world coords)",
	&outobj2.prof_bar_lengtherrw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.error;src.morph.scLength;stat.fit.param", "deg"},
  {"BAR_ASPECT_IMAGE", "Bar aspect ratio from fitting",
	&outobj2.prof_bar_aspect, H_FLOAT, T_FLOAT, "%6.4f", "",
	"phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"BAR_ASPECTERR_IMAGE", "RMS error on fitted bar aspect ratio",
	&outobj2.prof_bar_aspecterr, H_FLOAT, T_FLOAT, "%6.4f", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param;instr.det", ""},
  {"BAR_ASPECT_WORLD", "Bar aspect ratio from fitting",
	&outobj2.prof_bar_aspectw, H_FLOAT, T_FLOAT, "%12.7g", "",
	"phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"BAR_ASPECTERR_WORLD", "RMS error on fitted bar aspect ratio",
	&outobj2.prof_bar_aspecterrw, H_FLOAT, T_FLOAT, "%12.7g", "",
	"stat.error;phys.size.axisRatio;src.morph;stat.fit.param", ""},
  {"BAR_POSANG", "Bar true position angle (CCW/disk maj.axis) from fitting",
	&outobj2.prof_bar_posang, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.bodyrc.long;src.morph;stat.fit.param", "deg"},
  {"BAR_POSANGERR", "RMS error on fitted true bar position angle",
	&outobj2.prof_bar_posangerr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.bodyrc.long;src.morph;stat.fit.param", "deg"},
  {"BAR_THETA_IMAGE", "Bar projected angle (CCW/x) from fitting",
	&outobj2.prof_bar_theta, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"BAR_THETAERR_IMAGE", "RMS error on fitted bar projected angle",
	&outobj2.prof_bar_thetaerr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"BAR_THETA_WORLD", "Bar projected angle (CCW/world-x) from fitting",
	&outobj2.prof_bar_thetaw, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param", "deg"},
  {"BAR_THETAERR_WORLD", "RMS error on fitted bar projected angle",
	&outobj2.prof_bar_thetaerrw, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param", "deg"},
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
  {"FLUXERR_ARMS", "RMS error on fitted total flux in spiral arms",
	&outobj2.prof_arms_fluxerr, H_FLOAT, T_FLOAT, "%12.g", "count",
	"stat.error;phot.count;stat.fit.param", "ct"},
  {"MAG_ARMS", "Total magnitude in spiral arms from fitting",
	&outobj2.prof_arms_mag, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"phot.mag;stat.fit.param", "mag"},
  {"MAGERR_ARMS", "RMS error on fitted total magnitude in spiral arms",
	&outobj2.prof_arms_magerr, H_FLOAT, T_FLOAT, "%8.4f", "mag",
	"stat.error;phot.mag;stat.fit.param", "mag"},
  {"ARMS_SCALE_IMAGE", "Spiral arms scale length from fitting",
	&outobj2.prof_arms_scale, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"ARMS_SCALEERR_IMAGE", "RMS error on fitted spiral arms scale length",
	&outobj2.prof_arms_scaleerr, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"stat.error;src.morph.scLength;stat.fit.param;instr.det", "pix"},
  {"ARMS_SCALE_WORLD", "Spiral arms scale length from fitting",
	&outobj2.prof_arms_scalew, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"src.morph.scLength;stat.fit.param", "deg"},
  {"ARMS_SCALEERR_WORLD", "RMS error on fitted spiral arm scale length",
	&outobj2.prof_arms_scaleerrw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.error;src.morph.scLength;stat.fit.param", "deg"},
  {"ARMS_POSANG", "Pos. angle (CCW/disk maj.axis) of spiral arms from fitting",
	&outobj2.prof_arms_posang, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.bodyrc.long;src.morph;stat.fit.param", "deg"},
  {"ARMS_POSANGERR", "RMS error on fitted spiral arm position angle",
	&outobj2.prof_arms_posangerr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.bodyrc.long;src.morph;stat.fit.param", "deg"},

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
  {"ARMS_PITCH", "Pitch angle of spiral arms from fitting",
	&outobj2.prof_arms_pitch, H_FLOAT, T_FLOAT, "%+7.3f", "deg",
	"pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"ARMS_PITCHERR", "RMS error on fitted spiral arm pitch angle",
	&outobj2.prof_arms_pitcherr, H_FLOAT, T_FLOAT, "%7.3f", "deg",
	"stat.error;pos.posAng;src.morph;stat.fit.param;instr.det", "deg"},
  {"ARMS_START_IMAGE", "Starting radius of spiral arms from fitting",
	&outobj2.prof_arms_start, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"pos.distance;src.morph;stat.fit.param;instr.det", "pix"},
  {"ARMS_STARTERR_IMAGE", "RMS error on fitted spiral arm starting radius",
	&outobj2.prof_arms_starterr, H_FLOAT, T_FLOAT, "%10.4f", "pixel",
	"stat.error;pos.distance;src.morph;stat.fit.param;instr.det", "pix"},
  {"ARMS_START_WORLD", "Starting radius of spiral arms from fitting",
	&outobj2.prof_arms_startw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"pos.distance;src.morph;stat.fit.param", "deg"},
  {"ARMS_STARTERR_WORLD", "RMS error on spiral arm starting radius",
	&outobj2.prof_arms_starterrw, H_FLOAT, T_FLOAT, "%12.7g", "deg",
	"stat.error;pos.distance;src.morph;stat.fit.param", "deg"},
  {"ARMS_QUADFRAC", "Fraction of spiral arms in quadrature from fitting",
	&outobj2.prof_arms_quadfrac, H_FLOAT, T_FLOAT, "%6.4f", "deg",
	"phot.count;arith.ratio;src.morph;stat.fit.param", "deg"},
  {"ARMS_QUADFRACERR", "RMS error on fitted spiral arms quadrature fraction",
	&outobj2.prof_arms_quadfracerr, H_FLOAT, T_FLOAT, "%6.4f", "deg",
	"stat.error;phot.count;arith.ratio;src.morph;stat.fit.param", "deg"},
*/
