.. File Weighting.rst

.. include:: global.rst

Weighting
=========

The noise level in astronomical images is often fairly constant, that is, constant values for the gain, the background noise and the detection thresholds can be used over the whole frame.
Unfortunately in some cases, like strongly vignetted or composited images, this approximation is no longer good enough.
This leads to detecting clusters of detected noise peaks in the noisiest parts of the image, or missing obvious objects in the most sensitive ones.
|SExtractor| is able to handle images with variable noise.
It does it through *weight maps*, which are frames having the same size as the images where objects are detected or measured, and which describe the noise intensity at each pixel.
These maps are internally stored in units of *absolute variance* (in :math:`ADU^2`)
We employ the generic term *weight map* because these maps can also be interpreted as quality index maps: infinite variance (:math:`\ge 10^{30}` by definition in |SExtractor|) means that the related pixel in the science frame is totally unreliable and should be ignored.
The variance format was adopted as it linearizes most of the operations done over weight maps (see below).

This means that the noise covariances between pixels are ignored. Although raw CCD images have essentially white noise, this is not the case for warped images, for which resampling may induce a strong correlation between neighboring pixels.
In theory, all non-zero covariances within the geometrical limits of the analysed patterns should be taken into account to derive thresholds or error estimates.
Fortunately, the correlation length of the noise is often smaller than the patterns to be detected or measured, and constant over the image.
In that case one can apply a simple *fudge factor* to the estimated variance to account for correlations on small scales.
This proves to be a good approximation in general, although it certainly leads to underestimations for the smallest patterns.

.. _weight-map_format:

Weight-map formats
------------------

|SExtractor| accepts in input, and converts to its internal variance format, several types of weight-maps.
This is controlled through the ``WEIGHT_TYPE`` configuration keyword.
Those weight-maps can either be read from a FITS file, whose name is specified by the ``WEIGHT_IMAGE`` keyword, or computed internally.
Valid ``WEIGHT_TYPE`` values are:

* ``NONE`` :
    No weighting is applied.
    The related ``WEIGHT_IMAGE`` and ``WEIGHT_THRESH`` (see below) parameters are ignored.


* ``BACKGROUND`` :
    The science image itself is used to compute internally a variance map (the related ``WEIGHT_IMAGE`` parameter is ignored).
    Robust (:math:`3\sigma`-clipped) variance estimates are first computed within the same background meshes as the :ref:`background model <background_model>` [#f1]_.
    The resulting low-resolution variance map is then bicubic-spline-interpolated on the fly to produce the actual full-size variance map.
    A check-image with ``CHECKIMAGE_TYPE`` ``MINIBACK_RMS`` can be requested to examine the low-resolution variance map.


* ``MAP_RMS`` :
    The FITS image specified by the ``WEIGHT_IMAGE`` file name must contain a weight-map in units of absolute standard deviations (in ADUs per pixel).


* ``MAP_VAR`` :
    The FITS image specified by the ``WEIGHT_IMAGE`` file name must contain a weight-map in units of relative variance.
    A robust scaling to the appropriate absolute level is then performed by comparing this variance map to an internal, low-resolution, absolute variance map built from the science image itself.


* ``MAP_WEIGHT`` :
    The FITS image specified by the ``WEIGHT_IMAGE`` file name must contain a weight-map in units of relative weights.
    The data are converted to variance units (by definition variance :math:`\propto\frac{1}{weight}`), and scaled as for ``MAP_VAR``.
    ``MAP_WEIGHT`` is the most commonly used type of weight-map: a flat-field, for example, is generally a good approximation to a perfect weight-map.

.. _effect_of_weighting:

Effect of weighting
-------------------

Weight-maps modify the working of |SExtractor| in the following respects:

#. Bad pixels are discarded from the background statistics.
   If more than :math:`50\%` of the pixels in a background mesh are bad, the local background value and the background standard deviation are replaced by interpolation of the nearest valid meshes.

#. The detection threshold *t* above the local sky background is adjusted for each pixel *i* with variance :math:`\sigma^2_i`: :math:`t_i=`\ ``DETECT_THRESH``\ :math:`\times\sqrt{\sigma^2_i}`, where ``DETECT_THRESH`` is expressed in units of standard deviations of the background noise.
   Pixels with variance above the threshold set with the ``WEIGHT_THRESH`` parameter are therefore simply not detected.
   This may result in splitting objects crossed by a group of bad pixels.
   Interpolation should be used to avoid this problem.
   If convolution filtering is applied for detection, the variance map is convolved too.
   This yields optimum scaling of the detection threshold in the case where noise is uncorrelated from pixel to pixel.
   Non-linear filtering operations (like those offered by artificial retinae) are not affected.

#. The cleaning process takes into account the exact individual thresholds assigned to each pixel for deciding about the fate of faint detections.

#. Error estimates like ``FLUXISO_ERR`` , ``ERRA_IMAGE`` , etc. also make use of individual variances.
   Local background-noise standard deviation is simply set to :math:`\sqrt{\sigma^2_i}`.
   In addition, if the ``WEIGHT_GAIN`` parameter is set to ``Y`` (which is the default), it is assumed that the local pixel gain (i.e., the conversion factor from photo-electrons to ADUs) is inversely proportional to :math:`\sigma^2_i`, the median value over the image being set by the ``GAIN`` configuration parameter.
   In other words, it is then supposed that the changes in noise intensities seen over the images are due to gain changes.
   This is the most common case: correction for vignetting, or coverage depth.
   When this is not the case, for instance when changes are purely dominated by those of the read-out noise, ``WEIGHT_GAIN`` shall be set to ``N``.

#. Finally, pixels with weights beyond ``WEIGHT_THRESH`` are treated just like pixels discarded by the masking process.

Combining weight maps
---------------------

All the weighting options listed in :ref:`weight-map_format` can be applied separately to detection and measurement images (:ref:`using-sextractor`), even if some combinations may not always make sense.
For instance, the following set of configuration lines:

``WEIGHT_IMAGE`` *rms.fits*, *weight.fits*

``WEIGHT_TYPE`` ``MAP_RMS``, ``MAP_WEIGHT``

will load the FITS file *rms.fits* and use it as an RMS map for adjusting the detection threshold and cleaning, while the *weight.fits* weight map will only be used for scaling the error estimates on measurements.
This can be done in single- as well as in dual-image mode (see :ref:`using-sextractor`).
``WEIGHT_IMAGE`` can be ignored for ``BACKGROUND`` ``WEIGHT_TYPE``.
It is of course possible to use weight-maps for detection or for measurement only.
The following configuration:

``WEIGHT_IMAGE`` *weight.fits*

``WEIGHT_TYPE`` ``NONE``, ``MAP_WEIGHT``

will apply weighting only for measurements; detection and cleaning operations will remain unaffected.

.. _flags_weight_def:

Weight-map flags: :param:`FLAGS_WEIGHT`
---------------------------------------

The :param:`FLAGS_WEIGHT` catalog parameter flags various issues which may happen during the weighting process (see the :ref:`flagging` section for details on how flags are managed in |SExtractor|):

.. _flags_weight_table:

.. csv-table:: :param:`FLAGS_WEIGHT` description
  :header: "Value", "Meaning"
  :widths: 3 60

  1, "at least one measurement-image weight with a value below ``WEIGHT_THRESH`` *overlaps* the isophotal footprint of the detected object"
  2, "at least one measurement-image weight with a value below ``WEIGHT_THRESH`` *touches* the isophotal footprint of the detected object"

.. rubric:: Footnotes

.. [#f1] The mesh-filtering procedures act on the variance map, too.


