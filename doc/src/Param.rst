.. File Param.rst

.. include:: global.rst

The measurement (or catalog) parameter file
===========================================

In addition to the configuration file detailed above, |SExtractor| requires a file containing the list of measurements ("catalog parameters") that will be listed in the output catalog for every detection. This allows the software to compute only the measurements that are needed. The name of this catalog parameter file is traditionally suffixed with ``.param``, and must be specified using the :param:`PARAMETERS_NAME` config parameter. The full set of parameters can be queried with the command

.. code-block:: console

 $ sex -dp

Format
------

The format of the catalog parameter list is ASCII, and there must be
*a single keyword per line*. Presently two kinds of keywords are
recognized by |SExtractor|: scalars and vectors. Scalars, like :param:`X_IMAGE`,
produce single numbers in the output catalog. Vectors, like :param:`MAG_APER(4)`
or :param:`VIGNET(15,15)`, produce arrays of numbers. The ordering of measurements
in the output catalog is identical to that of the keywords in the parameter
list. Comments are allowed, they must begin with a :param:`#`.

Variants
--------

For many catalog parameters, especially those related to flux,
position, or shape, several variants of the same measurement are
available:

Fluxes
~~~~~~

Fluxes may be expressed in linear (ADU) units or as Pogson :cite:`1856MNRAS_17_12P` magnitudes. Flux measurements in ADUs are prefixed with :param:`FLUX_`, for example: :param:`FLUX_AUTO`, :param:`FLUX_ISO`, etc. Magnitudes are prefixed with :param:`MAG_` e.g., :param:`MAG_AUTO`, :param:`MAG_ISO`, ... In
|SExtractor| the magnitude :math:`m` of a source is derived from the flux
:math:`f`:

.. math::

   m = \left\{\begin{array}{ll}
   m_{ZP} -2.5 \log_{10} f\ &\mbox{if } f > 0\\
   99.0 &\mbox{otherwise},
   \end{array}\right.

where :math:`m_{ZP}` is the magnitude zero-point set with the
``MAG_ZEROPOINT`` configuration parameter.

Flux uncertainties
~~~~~~~~~~~~~~~~~~

Flux uncertainties follow a scheme similar to that of fluxes. Flux uncertainties are prefixed with :param:`FLUXERR_`, as in :param:`FLUXERR_AUTO` or :param:`FLUXERR_ISO`. Magnitude uncertainties start with :param:`MAGERR_`, for instance: :param:`MAGERR_AUTO`, :param:`MAGERR_ISO`,... Magnitude uncertainties :math:`\sigma_m` are derived from the estimated 1-\ :math:`\sigma` flux error :math:`\sigma_f`:

.. math::

   \sigma_m = \left\{\begin{array}{ll}
   (2.5/\ln 10) (\sigma_f/f)\ &\mbox{if } f > 0\\
   99.0 &\mbox{otherwise}.
   \end{array}\right.

.. _coord_suffix:

Positions and shapes
~~~~~~~~~~~~~~~~~~~~

Positions, distances and position angles are computed in pixel coordinates. They may be expressed in image pixels, world coordinates, or in celestial coordinates, depending on the suffix:

:param:`_IMAGE`
  Measurements are given in pixel coordinates, in units of pixels. For example: :param:`Y_IMAGE`, :param:`ERRAWIN_IMAGE`, :param:`THETA_IMAGE` etc. Following the FITS convention, in |SExtractor| the center of the first image pixel has coordinates (1.0,1.0). Position angles are counted from the *x* axis (axis 1), positive towards the *y* axis (axis 2)

:param:`_WORLD`
  Measurements are given in so-called “world coordinates”, converted from pixel coordinates using the local Jacobian of the transformation between both systems. This requires World Coordinate System (|WCS|_) metadata :cite:`2002AA_395_1061G` to be present in the FITS image header(s). Position angles are counted from the first world axis, positive towards the second world axis.

:param:`_SKY`, :param:`_J2000`, :param:`_B1950`
  Measurements are given in celestial (equatorial) coordinates, converted from pixel coordinates using the local Jacobian of the transformation between both systems. Positions and distances are in units of degrees. This requires celestial |WCS| metadata :cite:`2002AA_395_1077C` to be present in the FITS image header(s). :param:`_SKY` measurements are given in the native world coordinate system. :param:`_J2000` and :param:`_B1950` measurements are automatically converted from the native |WCS|, taking into account the change of reference frame. In all cases, positions angles are counted East-of-North.

Measurement parameter list
--------------------------

Below is an exhaustive list of all the measurement parameters known to
|SExtractor|. Please refer to the next sections for a detailed description
of their meaning.

.. csv-table:: |SExtractor| measurement parameters
  :class: targetparam
  :header: "Name", "Unit", "Description"
  :widths: 15 10 30

  _`NUMBER`,, Running object number
  _`ID_PARENT`,..., Parent ID (before deblending)
  _`EXT_NUMBER`,..., FITS extension number
  _`FLAGS`,..., Extraction flags
  _`FLUX_ISO`, count, :ref:`Isophotal flux <flux_iso_def>`
  _`FLUXERR_ISO`, count, :ref:`RMS error estimate for the isophotal flux <flux_iso_def>`
  _`MAG_ISO`, magnitude, :ref:`Isophotal magnitude <flux_iso_def>`
  _`MAGERR_ISO`, magnitude, :ref:`RMS error estimate for the isophotal magnitude <flux_iso_def>`
  _`FLUX_ISOCOR`, count, :ref:`Corrected isophotal flux <mag_isocor_def>`
  _`FLUXERR_ISOCOR`, count, :ref:`RMS error estimate for the corrected isophotal flux <mag_isocor_def>`
  _`MAG_ISOCOR`, magnitude, :ref:`Corrected isophotal magnitude <mag_isocor_def>`
  _`MAGERR_ISOCOR`, magnitude, :ref:`RMS error estimate for the corrected isophotal magnitude <mag_isocor_def>`
  _`FLUX_APER`, count, :ref:`Flux(es) within fixed circular aperture(s) <flux_aper_def>`
  _`FLUXERR_APER`, count, :ref:`RMS error estimate(s) for the flux(es) within fixed circular aperture(s) <flux_aper_def>`
  _`MAG_APER`, magnitude, :ref:`Circular aperture magnitude(s) <flux_aper_def>`
  _`MAGERR_APER`, magnitude, :ref:`RMS error estimate(s) for circular aperture magnitude(s) <flux_aper_def>`
  _`FLUX_AUTO`, count, :ref:`Kron-like automated aperture flux <flux_auto_def>`
  _`FLUXERR_AUTO`, count, :ref:`RMS error estimate for Kron-like automated aperture flux <flux_auto_def>`
  _`MAG_AUTO`, magnitude, :ref:`Kron-like automated aperture magnitude <flux_auto_def>`
  _`MAGERR_AUTO`, magnitude, :ref:`RMS error estimate for Kron-like automated aperture magnitude <flux_auto_def>`
  _`X_IMAGE`, pixel, :ref:`Isophotal image centroid along axis 1 (x) <pos_iso_def>`
  _`Y_IMAGE`, pixel, :ref:`Isophotal image centroid along axis 2 (y) <pos_iso_def>`
  _`ERRX2_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated variance of isophotal image centroid x coordinate <poserr_iso_def>`
  _`ERRY2_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated variance of isophotal image centroid y coordinate <poserr_iso_def>`
  _`ERRXY_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated covariance of isophotal image centroid x and y coordinates <poserr_iso_def>`
  _`ERRA_IMAGE`, pixel, :ref:`Major axis of the isophotal image centroid error ellipse <poserr_iso_def>`
  _`ERRB_IMAGE`, pixel, :ref:`Minor axis of the isophotal image centroid error ellipse <poserr_iso_def>`
  _`ERRTHETA_IMAGE`, degree, :ref:`Position angle of the isophotal image centroid ellipse <poserr_iso_def>`
  _`ERRCXX_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image centroid Cxx error ellipse parameter <poserr_iso_def>`
  _`ERRCYY_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image centroid Cyy error ellipse parameter <poserr_iso_def>`
  _`ERRCXY_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image centroid Cxy error ellipse parameter <poserr_iso_def>`
  _`XPEAK_IMAGE`, pixel, :ref:`x coordinate of the brightest pixel <pospeak_def>`
  _`YPEAK_IMAGE`, pixel, :ref:`y coordinate of the brightest pixel <pospeak_def>`
  _`XMIN_IMAGE`, pixel, :ref:`Minimum x coordinate among detected pixels <xyminmax_def>`
  _`YMIN_IMAGE`, pixel, :ref:`Minimum y coordinate among detected pixels <xyminmax_def>`
  _`XMAX_IMAGE`, pixel, :ref:`Maximum x coordinate among detected pixels <xyminmax_def>`
  _`YMAX_IMAGE`, pixel, :ref:`Maximum y coordinate among detected pixels <xyminmax_def>`
  _`XWIN_IMAGE`, pixel, :ref:`x coordinate of windowed image centroid <pos_win_def>`
  _`YWIN_IMAGE`, pixel, :ref:`y coordinate of windowed image centroid <pos_win_def>`
  _`ERRX2WIN_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated variance of windowed image centroid x coordinate <poserr_win_def>`
  _`ERRY2WIN_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated variance of windowed image centroid y coordinate <poserr_win_def>`
  _`ERRXYWIN_IMAGE`, pixel\ :sup:`2`, :ref:`Estimated covariance of windowed image centroid x and y coordinates <poserr_win_def>`
  _`ERRAWIN_IMAGE`, pixel, :ref:`Major axis of the windowed image centroid error ellipse <poserr_win_def>`
  _`ERRBWIN_IMAGE`, pixel, :ref:`Minor axis of the windowed image centroid error ellipse <poserr_win_def>`
  _`ERRTHETAWIN_IMAGE`, degree, :ref:`Position angle of the windowed image centroid ellipse <poserr_win_def>`
  _`ERRCXXWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image centroid Cxx error ellipse parameter <poserr_win_def>`
  _`ERRCYYWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image centroid Cyy error ellipse parameter <poserr_win_def>`
  _`ERRCXYWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image centroid Cxy error ellipse parameter <poserr_win_def>`
  _`X2_IMAGE`, pixel\ :sup:`2`, :ref:`Isophotal image 2nd order central moment in x <moments_iso_def>`
  _`Y2_IMAGE`, pixel\ :sup:`2`, :ref:`Isophotal image 2nd order central moment in y <moments_iso_def>`
  _`XY_IMAGE`, pixel\ :sup:`2`, :ref:`Isophotal image 2nd order central cross-moment in xy <moments_iso_def>`
  _`A_IMAGE`, pixel, :ref:`Isophotal image major axis <shape_iso_def>`
  _`B_IMAGE`, pixel, :ref:`Isophotal image minor axis <shape_iso_def>`
  _`THETA_IMAGE`, degree, :ref:`Isophotal image position angle<shape_iso_def>`
  _`ELONGATION`, ..., :ref:`A_IMAGE / B_IMAGE <elong_iso_def>`
  _`ELLIPTICITY`, ..., :ref:`1 - B_IMAGE / A_IMAGE <elong_iso_def>`
  _`CXX_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image Cxx ellipse parameter <ellipse_iso_def>`
  _`CYY_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image Cyy ellipse parameter <ellipse_iso_def>`
  _`CXY_IMAGE`, pixel\ :sup:`-2`, :ref:`Isophotal image Cxy ellipse parameter <ellipse_iso_def>`
  _`X2WIN_IMAGE`, pixel\ :sup:`2`, :ref:`Windowed image 2nd order central moment in x <moments_win_def>`
  _`Y2WIN_IMAGE`, pixel\ :sup:`2`, :ref:`Windowed image 2nd order central moment in y <moments_win_def>`
  _`XYWIN_IMAGE`, pixel\ :sup:`2`, :ref:`Windowed image 2nd order central cross-moment in xy <moments_win_def>`
  _`CXXWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image Cxx ellipse parameter <ellipse_win_def>`
  _`CYYWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image Cyy ellipse parameter <ellipse_win_def>`
  _`CXYWIN_IMAGE`, pixel\ :sup:`-2`, :ref:`Windowed image Cxy ellipse parameter <ellipse_win_def>`
  _`AWIN_IMAGE`, pixel, :ref:`Windowed image major axis <shape_win_def>`
  _`BWIN_IMAGE`, pixel, :ref:`Windowed image minor axis <shape_win_def>`
  _`THETAWIN_IMAGE`, degree, :ref:`Windowed image position angle <shape_win_def>`
  _`VECTOR_MODEL`, ..., :ref:`Model-fitting coefficients <models_def>`
  _`VECTOR_MODELERR`, ..., :ref:`Model-fitting coefficient uncertainties <models_def>`
  _`MATRIX_MODELERR`, ..., :ref:`Model-fitting covariance matrix <model_minimization_def>`
  _`CHI2_MODEL`, ..., :ref:`Reduced modified Chi2 of the fit <model_minimization_def>`
  _`FLAGS_MODEL`, ..., :ref:`Model-fitting flags <model_minimization_def>`
  _`NITER_MODEL`, ..., :ref:`Number of model-fitting iterations <model_minimization_def>`  
  _`FLUX_MODEL`, count, :ref:`Flux from model-fitting <models_def>`
  _`FLUXERR_MODEL`, count, :ref:`RMS error estimate for the model-fitting flux <models_def>`
  _`MAG_MODEL`, magnitude, :ref:`Magnitude from model-fitting <models_def>`
  _`MAGERR_MODEL`, count, :ref:`RMS error estimate for the model-fitting magnitude <models_def>`
  _`FLUX_MAX_MODEL`, count, :ref:`Peak model flux above the background <models_def>`
  _`FLUX_EFF_MODEL`, count, :ref:`Effective model flux above the background <models_def>`
  _`FLUX_EFF_MODEL`, count, :ref:`Mean effective model flux above the background <models_def>`
  _`MU_MAX_MODEL`, mag.arcsec\ :sup:`-2`, :ref:`Peak model surface brightness above the background <models_def>`
  _`MU_EFF_MODEL`, mag.arcsec\ :sup:`-2`, :ref:`Effective model surface brightness above the background <models_def>`
  _`MU_MEAN_MODEL`, mag.arcsec\ :sup:`-2`, :ref:`Mean effective model surface brightness above the background <models_def>`
  _`XMODEL_IMAGE`, pixel, :ref:`x coordinate from model-fitting <models_def>`
  _`YMODEL_IMAGE`, pixel, :ref:`y coordinate from model-fitting <models_def>`

..
  #XMODEL_WORLD             Fitted position along world x axis                        [deg]
  #YMODEL_WORLD             Fitted position along world y axis                        [deg]
  #ALPHAMODEL_SKY           Fitted position along right ascension  (native)           [deg]
  #DELTAMODEL_SKY           Fitted position along declination (native)                [deg]
  #ALPHAMODEL_J2000         Fitted position along right ascension (J2000)             [deg]
  #DELTAMODEL_J2000         Fitted position along declination (J2000)                 [deg]
  #ALPHAMODEL_B1950         Fitted position along right ascension (B1950)             [deg]
  #DELTAMODEL_B1950         Fitted position along declination (B1950)                 [deg]
  #ERRX2MODEL_IMAGE         Variance of fitted position along x                       [pixel**2]
  #ERRY2MODEL_IMAGE         Variance of fitted position along y                       [pixel**2]
  #ERRXYMODEL_IMAGE         Covariance of fitted position between x and y             [pixel**2]
  #ERRX2MODEL_WORLD         Variance of fitted position along X-WORLD (alpha)         [deg**2]
  #ERRY2MODEL_WORLD         Variance of fitted position along Y-WORLD (delta)         [deg**2]
  #ERRXYMODEL_WORLD         Covariance of fitted position X-WORLD/Y-WORLD             [deg**2]
  #ERRCXXMODEL_IMAGE        Cxx error ellipse parameter of fitted position            [pixel**(-2)]
  #ERRCYYMODEL_IMAGE        Cyy error ellipse parameter of fitted position            [pixel**(-2)]
  #ERRCXYMODEL_IMAGE        Cxy error ellipse parameter of fitted position            [pixel**(-2)]
  #ERRCXXMODEL_WORLD        Cxx fitted error ellipse parameter (WORLD units)          [deg**(-2)]
  #ERRCYYMODEL_WORLD        Cyy fitted error ellipse parameter (WORLD units)          [deg**(-2)]
  #ERRCXYMODEL_WORLD        Cxy fitted error ellipse parameter (WORLD units)          [deg**(-2)]
  #ERRAMODEL_IMAGE          RMS error of fitted position along major axis             [pixel]
  #ERRBMODEL_IMAGE          RMS error of fitted position along minor axis             [pixel]
  #ERRTHETAMODEL_IMAGE      Error ellipse pos.angle of fitted position (CCW/x)        [deg]
  #ERRAMODEL_WORLD          World RMS error of fitted position along major axis       [deg]
  #ERRBMODEL_WORLD          World RMS error of fitted position along minor axis       [deg]
  #ERRTHETAMODEL_WORLD      Error ellipse pos.angle of fitted position (CCW/world-x)  [deg]
  #ERRTHETAMODEL_SKY        Native fitted error ellipse pos. angle (east of north)    [deg]
  #ERRTHETAMODEL_J2000      J2000 fitted error ellipse pos. angle (east of north)     [deg]
  #ERRTHETAMODEL_B1950      B1950 fitted error ellipse pos. angle (east of north)     [deg]
  #X2MODEL_IMAGE            Variance along x from model-fitting                       [pixel**2]
  #Y2MODEL_IMAGE            Variance along y from model-fitting                       [pixel**2]
  #XYMODEL_IMAGE            Covariance between x and y from model-fitting             [pixel**2]
  #ELLIP1MODEL_IMAGE        Ellipticity component from model-fitting                 
  #ELLIP2MODEL_IMAGE        Ellipticity component from model-fitting                 
  #POLAR1MODEL_IMAGE        Ellipticity component (quadratic) from model-fitting     
  #POLAR2MODEL_IMAGE        Ellipticity component (quadratic) from model-fitting     
  #ELLIP1ERRMODEL_IMAGE     Ellipticity component std.error from model-fitting       
  #ELLIP2ERRMODEL_IMAGE     Ellipticity component std.error from model-fitting       
  #ELLIPCORRMODEL_IMAGE     Corr.coeff between ellip.components from model-fitting   
  #POLAR1ERRMODEL_IMAGE     Polarisation component std.error from model-fitting      
  #POLAR2ERRMODEL_IMAGE     Polarisation component std.error from model-fitting      
  #POLARCORRMODEL_IMAGE     Corr.coeff between polar. components from fitting        
  #X2MODEL_WORLD            Variance along X-WORLD (alpha) from model-fitting         [deg**2]
  #Y2MODEL_WORLD            Variance along Y_WORLD (delta) from model-fitting         [deg**2]
  #XYMODEL_WORLD            Covariance between X-WORLD and Y-WORLD from model-fitting [deg**2]
  #ELLIP1MODEL_WORLD        Ellipticity component from model-fitting                 
  #ELLIP2MODEL_WORLD        Ellipticity component from model-fitting                 
  #POLAR1MODEL_WORLD        Polarisation component from model-fitting                
  #POLAR2MODEL_WORLD        Polarisation component from model-fitting                
  #ELLIP1ERRMODEL_WORLD     Ellipticity component std.error from model-fitting       
  #ELLIP2ERRMODEL_WORLD     Ellipticity component std.error from model-fitting       
  #ELLIPCORRMODEL_WORLD     Corr.coeff between ellip.components from model-fitting   
  #POLAR1ERRMODEL_WORLD     Polarisation component std.error from model-fitting      
  #POLAR2ERRMODEL_WORLD     Polarisation component std.error from model-fitting      
  #POLARCORRMODEL_WORLD     Corr.coeff between polar. components from fitting        
  #CXXMODEL_IMAGE           Cxx ellipse parameter from model-fitting                  [pixel**(-2)]
  #CYYMODEL_IMAGE           Cyy ellipse parameter from model-fittinh                  [pixel**(-2)]
  #CXYMODEL_IMAGE           Cxy ellipse parameter from model-fitting                  [pixel**(-2)]
  #CXXMODEL_WORLD           Cxx ellipse parameter (WORLD) from model-fitting          [deg**(-2)]
  #CYYMODEL_WORLD           Cyy ellipse parameter (WORLD) from model-fitting          [deg**(-2)]
  #CXYMODEL_WORLD           Cxy ellipse parameter (WORLD) from model-fitting          [deg**(-2)]
  #AMODEL_IMAGE             Model RMS along major axis                                [pixel]
  #BMODEL_IMAGE             Model RMS along minor axis                                [pixel]
  #THETAMODEL_IMAGE         Model position angle (CCW/x)                              [deg]
  #AMODEL_WORLD             Model RMS along major axis (WORLD units)                  [deg]
  #BMODEL_WORLD             Model RMS along minor axis (WORLD units)                  [deg]
  #THETAMODEL_WORLD         Model position angle (CCW/WORLD-x)                        [deg]
  #THETAMODEL_SKY           Model position angle (east of north) (native)             [deg]
  #THETAMODEL_J2000         Model position angle (east of north) (J2000)              [deg]
  #THETAMODEL_B1950         Model position angle (east of north) (B1950)              [deg]
  #SPREAD_MODEL             Spread parameter from model-fitting                      

