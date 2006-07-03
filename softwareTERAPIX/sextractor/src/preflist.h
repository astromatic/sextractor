 /*
 				preflist.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Keywords for the configuration file.
*
*	Last modify:	03/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/


#include "key.h"

#ifndef _XML_H_
#include "xml.h"
#endif

#ifdef  USE_THREADS
#define THREADS_PREFMAX THREADS_NMAX
#else
#define THREADS_PREFMAX 65535
#endif

/*-------------------------------- initialization ---------------------------*/
int	idummy;

pkeystruct key[] =
 {
  {"ANALYSIS_THRESH", P_FLOATLIST, prefs.thresh, 0,0, -BIG, BIG,
    {""}, 1, 2, &prefs.nthresh},
  {"ASSOC_DATA", P_INTLIST, prefs.assoc_data, 0, 1000000,0.0,0.0,
    {""}, 1,MAXLIST, &prefs.nassoc_data},
  {"ASSOC_NAME", P_STRING, prefs.assoc_name},
  {"ASSOC_PARAMS", P_INTLIST, prefs.assoc_param, 1, 1000000,0.0,0.0,
    {""}, 2,3, &prefs.nassoc_param},
  {"ASSOC_RADIUS", P_FLOAT, &prefs.assoc_radius, 0,0, 1e-10,1e+10},
  {"ASSOC_TYPE", P_KEY, &prefs.assoc_type, 0,0, 0.0,0.0,
   {"FIRST", "NEAREST", "MEAN", "MAG_MEAN", "SUM", "MAG_SUM",
   "MIN", "MAX", ""}},
  {"ASSOCSELEC_TYPE", P_KEY, &prefs.assocselec_type, 0,0, 0.0,0.0,
   {"ALL","MATCHED","-MATCHED",""}},
  {"BACK_FILTERSIZE", P_INTLIST, prefs.backfsize, 1,11, 0.0,0.0,
    {""}, 1,2, &prefs.nbackfsize},
  {"BACK_FILTTHRESH", P_FLOAT, &prefs.backfthresh, 0,0, 0.0,BIG},
  {"BACKPHOTO_THICK", P_INT, &prefs.pback_size, 1, 256},
  {"BACKPHOTO_TYPE", P_KEY, &prefs.pback_type, 0,0, 0.0,0.0,
   {"GLOBAL","LOCAL",""}},
  {"BACK_SIZE", P_INTLIST, prefs.backsize, 1,2000000000, 0.0,0.0,
    {""}, 1,2, &prefs.nbacksize},
  {"BACK_TYPE", P_KEYLIST, prefs.back_type, 0,0, 0.0,0.0,
   {"AUTO","MANUAL",""},
    1, 2, &prefs.nback_type},
  {"BACK_VALUE", P_FLOATLIST, prefs.back_val, 0,0, -BIG,BIG,
   {""}, 1, 2, &prefs.nback_val},
  {"CATALOG_NAME", P_STRING, prefs.cat_name},
  {"CATALOG_TYPE", P_KEY, &prefs.cat_type, 0,0, 0.0,0.0,
   {"NONE", "ASCII","ASCII_HEAD", "ASCII_SKYCAT", "FITS_LDAC", "FITS_TPX",
	"FITS_1.0",""}},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "IDENTICAL",
   "BACKGROUND", "BACKGROUND_RMS", "MINIBACKGROUND",
   "MINIBACK_RMS", "-BACKGROUND",
   "FILTERED", "OBJECTS", "APERTURES", "SEGMENTATION", "ASSOC",
   "-OBJECTS", "-PSF_PROTOS", "PSF_PROTOS",
   "-PC_CONVPROTOS", "PC_CONVPROTOS", "PC_PROTOS", ""},
   0, 17, &prefs.ncheck_type},
  {"CLEAN", P_BOOL, &prefs.clean_flag},
  {"CLEAN_PARAM", P_FLOAT, &prefs.clean_param, 0,0, 0.1,10.0},
  {"DEBLEND_MINCONT", P_FLOAT, &prefs.deblend_mincont, 0,0, 0.0,1.0},
  {"DEBLEND_NTHRESH", P_INT, &prefs.deblend_nthresh, 1,64},
  {"DETECT_MINAREA", P_INT, &prefs.ext_minarea, 1,1000000},
  {"DETECT_THRESH", P_FLOATLIST, prefs.dthresh, 0,0, -BIG, BIG,
   {""}, 1, 2, &prefs.ndthresh},
  {"DETECT_TYPE", P_KEY, &prefs.detect_type, 0,0, 0.0,0.0,
   {"CCD","PHOTO",""}},
  {"FILTER", P_BOOL, &prefs.filter_flag},
  {"FILTER_NAME", P_STRING, prefs.filter_name},
  {"FILTER_THRESH", P_FLOATLIST, prefs.filter_thresh, 0,0,-BIG,BIG,
   {""}, 0, 2, &prefs.nfilter_thresh},
  {"FITS_UNSIGNED", P_BOOL, &prefs.fitsunsigned_flag},
  {"FLAG_IMAGE", P_STRINGLIST, prefs.fimage_name, 0,0,0.0,0.0,
    {""}, 0, MAXFLAG, &prefs.nfimage_name},
  {"FLAG_TYPE",  P_KEYLIST, prefs.flag_type, 0,0, 0.0,0.0,
   {"OR","AND","MIN", "MAX", "MOST",""}, 0, MAXFLAG, &idummy},
  {"GAIN", P_FLOAT, &prefs.gain, 0,0, 0.0, 1e+30},
  {"INTERP_MAXXLAG", P_INTLIST, prefs.interp_xtimeout, 1,1000000, 0.0,0.0,
   {""}, 1, 2, &prefs.ninterp_xtimeout},
  {"INTERP_MAXYLAG", P_INTLIST, prefs.interp_ytimeout, 1,1000000, 0.0,0.0,
   {""}, 1, 2, &prefs.ninterp_ytimeout},
  {"INTERP_TYPE", P_KEYLIST, prefs.interp_type, 0,0, 0.0,0.0,
   {"NONE","VAR_ONLY","ALL",""}, 1, 2, &prefs.ninterp_type},
  {"MAG_GAMMA", P_FLOAT, &prefs.mag_gamma, 0,0, 1e-10,1e+30},
  {"MAG_ZEROPOINT", P_FLOAT, &prefs.mag_zeropoint, 0,0, -100.0, 100.0},
  {"MAMA_CORFLEX", P_FLOAT, &prefs.mama_corflex, 0,0, -1.0,1.0},
  {"MASK_TYPE", P_KEY, &prefs.mask_type, 0,0, 0.0,0.0,
   {"NONE","BLANK","CORRECT",""}},
  {"MEMORY_BUFSIZE", P_INT, &prefs.mem_bufsize, 8, 65534},
  {"MEMORY_OBJSTACK", P_INT, &prefs.clean_stacksize, 16,65536},
  {"MEMORY_PIXSTACK", P_INT, &prefs.mem_pixstack, 1000, 10000000},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PARAMETERS_NAME", P_STRING, prefs.param_name},
  {"PHOT_APERTURES", P_FLOATLIST, prefs.apert, 0,0, 0.0,2*MAXPICSIZE,
   {""}, 1, MAXNAPER, &prefs.naper},
  {"PHOT_AUTOPARAMS", P_FLOATLIST, prefs.autoparam, 0,0, 0.0,10.0,
   {""}, 2,2, &prefs.nautoparam},
  {"PHOT_AUTOAPERS", P_FLOATLIST, prefs.autoaper, 0,0, 0.0,1e6,
   {""}, 2,2, &prefs.nautoaper},
  {"PHOT_FLUXFRAC", P_FLOATLIST, prefs.flux_frac, 0,0, 1e-6, 1.0,
   {""}, 1, MAXNAPER, &prefs.nflux_frac},
  {"PHOT_PETROPARAMS", P_FLOATLIST, prefs.petroparam, 0,0, 0.0,10.0,
   {""}, 2,2, &prefs.npetroparam},
  {"PIXEL_SCALE", P_FLOAT, &prefs.pixel_scale, 0,0, 0.0, 1e+10},
  {"PSF_NAME", P_STRINGLIST, prefs.psf_name, 0,0, 0.0,0.0,
   {""}, 1, 2, &prefs.npsf_name},	/*?*/
  {"PSF_NMAX", P_INT, &prefs.psf_npsfmax, 1, PSF_NPSFMAX},
  {"PSFDISPLAY_TYPE", P_KEY, &prefs.psfdisplay_type, 0,0, 0.0,0.0,
   {"SPLIT","VECTOR",""}},
  {"SATUR_LEVEL", P_FLOAT, &prefs.satur_level, 0,0, -1e+30, 1e+30},
  {"SEEING_FWHM", P_FLOAT, &prefs.seeing_fwhm, 0,0, 1e-10, 1e+10},
  {"SOM_NAME", P_STRING, prefs.som_name},
  {"STARNNW_NAME", P_STRING, prefs.nnw_name},
  {"THRESH_TYPE", P_KEYLIST, prefs.thresh_type, 0,0, 0.0,0.0,
   {"RELATIVE","ABSOLUTE"},
    1, 2, &prefs.nthresh_type},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL", "EXTRA_WARNINGS", "FULL",""}},
  {"WEIGHT_GAIN", P_BOOL, &prefs.weightgain_flag},
  {"WEIGHT_IMAGE", P_STRINGLIST, prefs.wimage_name, 0,0,0.0,0.0,
    {""}, 0, MAXIMAGE, &prefs.nwimage_name},
  {"WEIGHT_THRESH", P_FLOATLIST, prefs.weight_thresh, 0,0, 0.0, BIG,
   {""}, 0, 2, &prefs.nweight_thresh},
  {"WEIGHT_TYPE", P_KEYLIST, prefs.weight_type, 0,0, 0.0,0.0,
   {"NONE","BACKGROUND", "MAP_RMS", "MAP_VAR","MAP_WEIGHT", ""},
   0, MAXIMAGE, &prefs.nweight_type},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {""}
 };

char		keylist[sizeof(key)/sizeof(pkeystruct)][32];
const char	notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#-------------------------------- Catalog ------------------------------------",
" ",
"CATALOG_NAME     test.cat       # name of the output catalog",
"CATALOG_TYPE     ASCII_HEAD     # \"NONE\",\"ASCII_HEAD\",\"ASCII\",\"FITS_1.0\"",
"                                # \"FITS_LDAC\" or \"FITS_TPX\"",
" ",
"PARAMETERS_NAME  default.param  # name of the file containing catalog contents",
" ",
"#------------------------------- Extraction ----------------------------------",
" ",
"DETECT_TYPE      CCD            # \"CCD\" or \"PHOTO\"",
"DETECT_MINAREA   5              # minimum number of pixels above threshold",
"*THRESH_TYPE      RELATIVE       # threshold type: \"RELATIVE\" (in sigmas)",
"*                                # or \"ABSOLUTE\" (in ADUs)",
"DETECT_THRESH    1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2",
"ANALYSIS_THRESH  1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2",
" ",
"FILTER           Y              # apply filter for detection (\"Y\" or \"N\")?",
"FILTER_NAME      default.conv   # name of the file containing the filter",
"*FILTER_THRESH                   # Threshold[s] for retina filtering",

" ",
"DEBLEND_NTHRESH  32             # Number of deblending sub-thresholds",
"DEBLEND_MINCONT  0.005          # Minimum contrast parameter for deblending",
" ",
"CLEAN            Y              # Clean spurious detections? (Y or N)?",
"CLEAN_PARAM      1.0            # Cleaning efficiency",
" ",
"MASK_TYPE        CORRECT        # type of detection MASKing: can be one of",
"                                # \"NONE\", \"BLANK\" or \"CORRECT\"",
" ",
"*#-------------------------------- WEIGHTing ----------------------------------",
"*",
"*WEIGHT_TYPE      NONE           # type of WEIGHTing: \"NONE\", \"BACKGROUND\",",
"*                                # \"MAP_RMS\", \"MAP_VAR\" or \"MAP_WEIGHT\"",
"*WEIGHT_IMAGE     weight.fits    # weight-map filename",
"*WEIGHT_GAIN      Y              # modulate gain (E/ADU) with weights? (Y/N)",
"*WEIGHT_THRESH                   # weight threshold[s] for bad pixels",
"*",
"*#-------------------------------- FLAGging -----------------------------------",
"*",
"*FLAG_IMAGE       flag.fits      # filename for an input FLAG-image",
"*FLAG_TYPE        OR             # flag pixel combination: \"OR\", \"AND\",",
"*                                # \"MIN\", \"MAX\", \"MOST\"",
"*",
"#------------------------------ Photometry -----------------------------------",
" ",
"PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels",
"PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>",
"PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,",
"                                # <min_radius>",
"*PHOT_AUTOAPERS   0.0,0.0        # <estimation>,<measurement> minimum apertures",
"*                                # for MAG_AUTO and MAG_PETRO",
"*PHOT_FLUXFRAC    0.5            # flux fraction[s] used for FLUX_RADIUS",
" ",
"SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation",
" ",
"MAG_ZEROPOINT    0.0            # magnitude zero-point",
"MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)",
"GAIN             0.0            # detector gain in e-/ADU",
"PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)",
" ",
"#------------------------- Star/Galaxy Separation ----------------------------",
" ",
"SEEING_FWHM      1.2            # stellar FWHM in arcsec",
"STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename",
" ",
"#------------------------------ Background -----------------------------------",
" ",
"*BACK_TYPE        AUTO           # \"AUTO\" or \"MANUAL\"",
"*BACK_VALUE       0.0            # Default background value in MANUAL mode",
"BACK_SIZE        64             # Background mesh: <size> or <width>,<height>",
"BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>",
" ",
"BACKPHOTO_TYPE   GLOBAL         # can be \"GLOBAL\" or \"LOCAL\"",
"*BACKPHOTO_THICK  24             # thickness of the background LOCAL annulus",
"*BACK_FILTTHRESH  0.0            # Threshold above which the background-",
"*                                # map filter operates",
" ",
"#------------------------------ Check Image ----------------------------------",
" ",
"CHECKIMAGE_TYPE  NONE           # can be one of \"NONE\", \"BACKGROUND\",",
"                                # \"MINIBACKGROUND\", \"-BACKGROUND\", \"OBJECTS\",",
"                                # \"-OBJECTS\", \"SEGMENTATION\", \"APERTURES\",",
"                                # or \"FILTERED\"",
"CHECKIMAGE_NAME  check.fits     # Filename for the check-image",
" ",
"#--------------------- Memory (change with caution!) -------------------------",
" ",
"MEMORY_OBJSTACK  3000           # number of objects in stack",
"MEMORY_PIXSTACK  300000         # number of pixels in stack",
"MEMORY_BUFSIZE   1024           # number of lines in buffer",
" ",
"*#------------------------------- ASSOCiation ---------------------------------",
"*",
"*ASSOC_NAME       sky.list       # name of the ASCII file to ASSOCiate",
"*ASSOC_DATA       2,3,4          # columns of the data to replicate (0=all)",
"*ASSOC_PARAMS     2,3,4          # columns of xpos,ypos[,mag]",
"*ASSOC_RADIUS     2.0            # cross-matching radius (pixels)",
"*ASSOC_TYPE       MAG_SUM        # ASSOCiation method: FIRST, NEAREST, MEAN,",
"*                                # MAG_MEAN, SUM, MAG_SUM, MIN or MAX",
"*ASSOCSELEC_TYPE  MATCHED        # ASSOC selection type: ALL, MATCHED or -MATCHED",
"*",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
"VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL",
"WRITE_XML        N              # Write XML file (Y/N)?",
"XML_NAME         sex.xml        # Filename for XML output",
"*XSL_URL          " XSL_URL,
"*                                # Filename for XSL style-sheet",
#ifdef USE_THREADS
"*NTHREADS         0              # Number of simultaneous threads for",
"*                                # the SMP version of " BANNER,
"*                                # 0 = automatic",
#else
"*NTHREADS         1              # 1 single thread",
#endif
"*",
"*FITS_UNSIGNED    N              # Treat FITS integer values as unsigned (Y/N)?",
"*INTERP_MAXXLAG   16             # Max. lag along X for 0-weight interpolation",
"*INTERP_MAXYLAG   16             # Max. lag along Y for 0-weight interpolation",
"*INTERP_TYPE      ALL            # Interpolation type: NONE, VAR_ONLY or ALL",
"*",
"*#--------------------------- Experimental Stuff -----------------------------",
"*",
"*MAMA_CORFLEX     3.3e-5         # MAMA correction factor",
"*PSF_NAME         default.psf    # File containing the PSF model",
"*PSF_NMAX         9              # Max.number of PSFs fitted simultaneously",
"*PSFDISPLAY_TYPE  SPLIT          # Catalog type for PSF-fitting: SPLIT or VECTOR",
"*SOM_NAME         default.som    # File containing Self-Organizing Map weights",
""
 };

