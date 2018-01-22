.. File Photom.rst

.. include:: global.rst

.. _photometry:

Aperture photometry
===================

Besides :ref:`isophotal <flux_iso_def>`, |PSF| and :ref:`model-fitting <models_def>` flux estimates, |SExtractor| can currently perform two types of flux measurements: :ref:`fixed-aperture <flux_aper_def>` and :ref:`adaptive-aperture <flux_auto_def>`.
For every :param:`FLUX_` measurement, an error estimate :param:`FLUXERR_`, a magnitude :param:`MAG_` and a magnitude error estimate :param:`MAGERR_` are also available: see :ref:`fluxes_and_magnitudes`.

An estimate of the error is available for each type of flux.
For aperture fluxes, the flux uncertainty is computed using

.. math::
  :label: fluxerr

  {\tt FLUXERR} = \sqrt{\sum_{i\in{\cal A}}\, (\sigma_i^2 + \frac{p_i}{g_i})}

where :math:`{\cal A}` is the set of pixels defining the photometric aperture, and :math:`\sigma_i`, :math:`p_i`, :math:`g_i` respectively the standard deviation of noise (in ADU) estimated from the local background, :math:`p_i` the measurement image pixel value subtracted from the background, and :math:`g_i` the effective detector gain in :math:`e^- / \mbox{ADU}` at pixel :math:`i`.
Note that this error estimate provides a lower limit of the true uncertainty, as it only takes into account photon and detector noise.

.. _flux_aper_def:

Fixed-aperture flux: :param:`FLUX_APER`
---------------------------------------

:param:`FLUX_APER` estimates the flux from the measurement image above the background inside a circular aperture.
The diameter of the aperture in pixels is defined by the ``PHOTOM_APERTURES`` configuration parameter.
It does not have to be an integer: each "regular" pixel is subdivided in :math:`5\times 5` sub-pixels before measuring the flux within the aperture.
If :param:`FLUX_APER` is provided as a vector :param:`FLUX_APER[n]`, at least :math:`n` apertures must be specified with the ``PHOTOM_APERTURES`` configuration parameter.

.. _flux_auto_def:

Automatic aperture flux: :param:`FLUX_AUTO`
-------------------------------------------

:param:`FLUX_AUTO` provides an estimate of the “total flux” by integrating pixel values within an adaptively scaled aperture.
|SExtractor|’s automatic aperture photometry routine derives from Kron’s “first moment” algorithm :cite:`1980ApJS_43_305K`:

#. An elliptical aperture is :ref:`defined by the second order moments of the object’s light distribution <ellipse_iso_def>`, with semi-major axis :math:`a={\tt A\_IMAGE}`, semi-minor axis :math:`b={\tt B\_IMAGE}`, and position angle :param:`THETA_IMAGE`.
#. The ellipse's major and minor axes are multiplied by 6 (which corresponds roughly to twice the size of the isophotal footprint on each axis).
#. Inside this elliptical aperture :math:`{\cal E}`, an analog of Kron's "first moment" is computed:

.. math::
  :label: kron_radius

  r_{\rm Kron} = \frac{\sum_{i\in\cal E} r_i\,p^{(d)}_i}{\sum_{i\in\cal E} p^{(d)}_i},

where :math:`p^{(d)}_i` is the pixel value *in the detection image*. :math:`r_i` is what we shall call the "reduced pseudo-radius" at pixel :math:`i` 

.. math::
  :label: reduced_radius

  r_i \equiv \sqrt{{\tt CXX\_IMAGE} \times \Delta x_i^2 + {\tt CYY\_IMAGE} \times \Delta y_i^2 + {\tt CXY\_IMAGE} \times \Delta x_i \Delta y_i},

where :math:`\Delta x_i` and  :math:`\Delta y_i` are the pixel coordinates relative to the detection centroid:

.. math::

  \begin{aligned}
  \Delta x_i & = x_i - {\tt X\_IMAGE}\\
  \Delta y_i & = y_i - {\tt Y\_IMAGE}.
  \end{aligned}


:cite:`1980ApJS_43_305K` and :cite:`1987AA_183_177I` have shown that for stars and galaxy profiles convolved with Gaussian seeing, :math:`\ge 90\%` of the flux is expected to lie inside a circular aperture of radius :math:`k r_{\rm Kron}` with :math:`k = 2`, almost independently of the magnitude.
Experiments have shown :cite:`1996AAS_117_393B` that this conclusion remains unchanged if one replaces the circular aperture with the "Kron elliptical aperture" :math:`{\cal K}` with reduced pseudo-radius :math:`k r_{\rm Kron}`.

:param:`FLUX_AUTO` is the sum of pixel values from the measurement image, subtracted from the local background, inside the Kron ellipse:

.. math::
  :label: flux_auto

  {\tt FLUX\_AUTO} = \sum_{i\in\cal K} p_i.

The quantity :math:`k r_{\rm Kron}`, known as the *Kron radius* (which in |SExtractor| is actually a "reduced pseudo-radius") is provided by the :param:`KRON_RADIUS`.
:math:`k = 2` defines a sort of balance between systematic and random errors.
By choosing a larger :math:`k = 2.5`, the mean fraction of flux lost drops from about 10% to 6%, at the expense of |SNR| in the measurement.
Very noisy objects may sometimes end up with a Kron ellipse being too small, even smaller that the isophotal footprint of the object itself. For this reason, |SExtractor| imposes a minimum size for the Kron radius, which cannot be less than :math:`r_{\rm Kron,min}`.
The user has full control over the parameters :math:`k` and :math:`r_{\rm Kron,min}` through the ``PHOT_AUTOPARAMS`` configuration parameters. ``PHOT_AUTOPARAMS`` is set by default to ``2.5,3.5``.

..
   .. figure:: ps/simlostflux.ps
   :alt:  Flux lost (expressed as a mean magnitude difference) with different faint-object photometry techniques as a function of total magnitude (see text). Only isolated galaxies (no blends) of the simulations have been considered.
   :width: 15.00000cm

   Flux lost (expressed as a mean magnitude difference) with different
   faint-object photometry techniques as a function of total magnitude
   (see text). Only isolated galaxies (no blends) of the simulations
   have been considered. 

.. hint::
  Aperture magnitudes are sensitive to crowding.
  In |SExtractor| v1, :param:`MAG_AUTO` measurements were not very robust in that respect.
  It was therefore suggested to replace the aperture magnitude by the corrected-isophotal one when an object is too close to its neighbors (2 isophotal radii for instance).
  This was done automatically when using the :param:`MAG_BEST` magnitude: :param:`MAG_BEST`\=\ :param:`MAG_AUTO` when it is sure that no neighbor can bias :param:`MAG_AUTO` by more than 10%, and :math:`{\tt MAG\_BEST} = {\tt MAG\_ISOCOR}` otherwise.
  Experience showed that the :param:`MAG_ISOCOR` and :param:`MAG_AUTO` magnitude would loose about the same fraction of flux on stars or compact galaxy profiles: around 0.06 % for default extraction parameters.
  The use of :param:`MAG_BEST` is now deprecated as :param:`MAG_AUTO` measurements are much more robust in versions 2.x of |SExtractor|.
  The first improvement is a crude subtraction of all the neighbors that have been detected around the measured source (``MASK_TYPE BLANK`` option).
  The second improvement is an automatic correction of parts of the aperture that are suspected to be contaminated by a neighbor.
  This is done by mirroring the opposite, cleaner side of the measurement ellipse if available (``MASK_TYPE CORRECT`` option, which is also the default).

..
  Figure [figphot] shows the mean loss of flux measured with isophotal (threshold 24.4 magnitude.arsec\ :sup:`-2`), corrected isophotal and automatic aperture photometry for simulated galaxies on a typical Schmidt-survey B\ :sub:`J` plate image.
  The automatic adaptive aperture photometry leads to the lowest loss of flux.

.. _flux_petro_def:

Petrosian aperture flux: :param:`FLUX_PETRO`
--------------------------------------------

Similar to :param:`FLUX_AUTO`, :param:`FLUX_PETRO` provides an estimate of the “total flux” by integrating pixel values within an adaptively scaled elliptical aperture. :param:`FLUX_PETRO`\ 's algorithm derives from Petrosian’s photometric estimator :cite:`1976ApJ_209L_1P,2001AJ_121_2358B,2001AJ_122_1104Y`:

#. An elliptical aperture is :ref:`defined by the second order moments of the object’s light distribution <ellipse_iso_def>`, with semi-major axis :math:`a={\tt A\_IMAGE}`, semi-minor axis :math:`b={\tt B\_IMAGE}`, and position angle :param:`THETA_IMAGE`.
#. The ellipse's major and minor axes are multiplied by 6 (which corresponds roughly to twice the size of the isophotal footprint on each axis).
#. Within this elliptical aperture :math:`{\cal E}`, the *Petrosian ratio* :math:`R_{\rm P}(r)` is computed:

.. math::
  :label: petrosian_ratio

  R_{\rm P}(r) = \frac{\sum_{0.9r < r_i < 1.1r} p^{(d)}_i}{\sum_{r_i < r} p^{(d)}_i} \frac{N_{r_i < r}}{N_{0.9r < r_i < 1.1r}},

where :math:`p^{(d)}_i` is the pixel value *in the detection image*. :math:`r_i` is the "reduced pseudo-radius" at pixel :math:`i` as defined in :eq:`reduced_radius`.
The *Petrosian ellipse* :math:`{\cal P}` is the ellipse with reduced pseudo-radius :math:`N_{\rm P}r_{\rm P}`, where :math:`r_{\rm P}` is defined by

.. math::
  :label: petrosian_radius

  R_{\rm P}(r_{\rm p}) \equiv 0.2

The quantity :math:`N_{\rm P}r_{\rm P}` is called *Petrosian radius* in |SExtractor|\ [#petro_radius]_ and is provided by the :param:`PETRO_RADIUS` catalog parameter.
The Petrosian factor :math:`N_{\rm P}` is set to 2.0 by default.
Very noisy objects may sometimes end up with a Petrosian ellipse being too small.
For this reason, |SExtractor| imposes a minimum size for the Petrosian radius, which cannot be less than :math:`r_{\rm P,min}`.
The user has full control over the parameters :math:`N_{\rm P}` and :math:`r_{\rm P,min}` through the ``PHOT_PETROPARAMS`` configuration parameters. ``PHOT_PETROPARAMS`` is set by default to ``2.0,3.5``.

The Petrosian flux is the sum of pixel values from the measurement image, subtracted from the local background, inside the Petrosian ellipse:

.. math::
  :label: flux_petro

  {\tt FLUX\_PETRO} = \sum_{i\in\cal P} p_i.

.. [#petro_radius]
   Some authors prefer to define the Petrosian radius as :math:`r_{\rm P}` instead of :math:`N_{\rm P}r_{\rm P}`.

Photographic photometry
-----------------------

In ``DETECT_TYPE PHOTO`` mode, SExtractor assumes that the response of the detector, over the dynamic range of the image, is logarithmic.
This is generally a good approximation for photographic density on deep exposures. 
Photometric procedures described above remain unchanged, except that for each pixel we apply first the transformation

.. math::
  :label: dtoi

  I = I_0\,10^{D/\gamma},

where :math:`\gamma` (``MAG_GAMMA``) is the contrast index of the emulsion, :math:`D` the original pixel value from the background-subtracted image, and :math:`I_0` is computed from the magnitude zero-point :math:`m_0`:

.. math::
  :label: m0toi0

  I_0 = \frac{\gamma}{\ln 10} \,10^{-0.4\, m_0}.

One advantage of using a density-to-intensity transformation relative to the local sky background is that it corrects (to some extent) large-scale inhomogeneities in sensitivity (see :cite:`1996PhDT_68B` for details).

..
  Magnitude uncertainties
  -----------------------

  In ``DETECT_TYPE PHOTO`` mode, things are slightly more complex.
  Making the assumption that plate-noise is the major contributor to photometric errors, and that it is roughly constant in density, one can write:

  .. math::
  :label: magerr_photo

   \Delta m = 1.0857 \,\ln 10\, {\sigma\over \gamma}\,
     \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}
   =2.5\,{\sigma\over \gamma}\,
   \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}

  where :math:`I(x,y)` is the contribution of pixel :math:`(x,y)` to the total flux :eq:`dtoi`. ``GAIN`` is ignored in ``PHOTO`` mode.

.. _background_def:

Local background
----------------

Almost all |SExtractor| measurements are done using background-subtracted pixel values.
In crowded fields, or in images where the background is irregular, the :ref:`background model <background_model>` may be significantly inaccurate, locally creating biases in photometric estimates.

The user has the possibility to force |SExtractor| to correct, for every detection, the background used to compute fluxes by setting the ``BACKPHOTO_TYPE`` configuration parameter to ``LOCAL`` (``GLOBAL`` is the default).
In ``LOCAL`` mode, a mean background residual level is estimated from background-subtracted pixel values within a "rectangular annulus" around the isophotal limits of the object.
The user can specify the thickness of the annulus, in pixels, with the ``BACKPHOTO_SIZE`` configuration parameter. The default thickness is ``24`` pixels.
The residual background level computed in ``LOCAL`` mode is used by |SExtractor| to correct all aperture photometry measurements, as well as basic surface brightness estimations such as :param:`FLUX_MAX`.
However in practice the ``BACKPHOTO_TYPE LOCAL`` option has not proven as being really beneficial to photometric accuracy, and it is generally advised to leave ``BACKPHOTO_TYPE`` set to ``GLOBAL``.

In both ``LOCAL`` and ``GLOBAL`` modes, the :param:`BACKGROUND` catalog parameter gives the value of the background estimated at the centroid of the object.


