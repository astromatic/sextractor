.. File Background.rst

.. include:: global.rst

.. _background_model:

Modeling the background
=======================

On linear detectors, the value measured at each pixel is the sum of a "background" signal and light coming from the sources of interest.
To be able to detect the faintest objects and make accurate measurements, |SExtractor| needs first computing a precise estimate of the background level at any position of the image: a *background map*.
Strictly speaking, there should be one background map per source, that is, what would the image look like if that very source was missing.
But, at least for detection, one can start by assuming that most discrete sources do not overlap too severely — which is generally the case for high galactic latitude fields —, and that the background varies smoothly across the field.

Background estimation
---------------------

To compute the background map, |SExtractor| makes a first pass through the pixel data, estimating the local background in each mesh of a rectangular grid that covers the whole frame.
The background estimator is a combination of :math:`\kappa\,\sigma` clipping and mode estimation, similar to Stetson’s |DAOPHOT|_ program :cite:`1987PASP_99_191S,1992ASPC_23_90D`.

Briefly, the local background histogram is clipped iteratively until convergence at :math:`\pm 3\sigma` around its median. The mode of the histogram is estimated using:

.. math::
  :label: sexbackmode

  \mbox{Mode} = 2.5 \times \mbox{Median} - 1.5 \times \mbox{Mean}.

Using simulated images, the expression above was found more accurate with clipped distributions :cite:`1996AAS_117_393B` than the usual approximation (e.g., :cite:`stuart2009kendall`):

.. math::
  :label: ksbackmode

  \mbox{Mode} = 3 \times \mbox{Median} - 2 \times \mbox{Mean}.

:numref:`fig_modevsmean` shows that the mode estimation in :eq:`sexbackmode` is considerably less affected by source crowding than a simple clipped mean :cite:`1981AJ_86_476J,1987AA_183_177I` but it is :math:`\approx 30\%` noisier. 
Obviously :eq:`sexbackmode` is not valid for any distribution; |SExtractor| falls back to a simple median for estimating the local background value if the mode and the median disagree by more than 30%.

.. _fig_modevsmean:

.. figure:: figures/modevsmean.*
   :figwidth: 100%
   :align: center

   Simulations of :math:`32\times32` pixels background meshes contamined by random Gaussian profiles.
   The true background lies at 0 ADUs.
   While being a bit noisier, the clipped "mode" gives a more robust estimate than the clipped mean in crowded regions. 

Background map
--------------

Once the grid is set up, a median filter can be applied to suppress possible local overestimations due to bright stars.
The final background map is simply a (natural) bicubic-spline interpolation between the meshes of the grid.
Median filtering helps reducing possible ringing effects of the bicubic-spline around bright features.

In parallel with the making of the background map, an *RMS background map*, that is, a map of the background noise standard deviation in the image is produced.
It may be used as an internal weight map if the ``WEIGHT_TYPE`` configuration parameter is to ``BACKGROUND`` (see :ref:`weight-map_format`).

Configuration and tuning
------------------------

.. note::
  All background configuration parameters also affect background-RMS maps.

The choice of the mesh size ``BACK_SIZE`` is very important.
If it is too small, the background estimation is affected by the presence of
objects and random noise.
Most importantly, part of the flux of the most extended objects can be absorbed into the background map.
If the mesh size is too large, it cannot reproduce the small scale variations of the background.
Therefore a good compromise must be found by the user.
Typically, for reasonably sampled images, a width [#recmesh]_ of 32 to 512 pixels works well.

The user has some control over background map filtering by specifying the size of the median filter ``BACK_FILTERSIZE``.
A width and height of 1 means that no filtering is applied to the background grid.
A :math:`3\times3` filtering is sufficient in most cases.
Larger dimensions may occasionally be used to compensate for small background mesh sizes, or in the presence of large image artifacts. 
In some specific cases it is desirable to median-filter only background meshes that have values exceeding some threshold above the filtered value.
This differential threshold is set by the ``BACK_FILTERTHRESH`` parameter, in ADUs.

By default, the computed background maps are automatically subtracted from input images.
However there are some situations where it is more appropriate to subtract a *constant* from the image (e.g., images with strongly non-Gaussian background noise |pdf|\ s).
The ``BACK_TYPE`` configuration parameter (set to ``AUTO`` by default) may be switched to ``MANUAL`` to force the value specified by the ``BACK_VALUE`` parameter to be subtracted from the input image, instead of the background map. ``BACK_VALUE`` is ``0`` by default.

Computing cost
--------------

The background estimation operation is generally |I/O|\ -bound, unless the image file already resides in the disk cache.

.. [#recmesh]
   It is possible to specify rectangular background meshes, although it is advised to use square ones, except in special cases (background varying rapidly along the x or y axis).
