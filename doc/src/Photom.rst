.. File Photom.rst

.. include:: global.rst

Photometry
==========

Besides |PSF| and :ref:`model-fitting <models_def>` flux estimates, |SExtractor| can currently perform four types of flux measurements: :ref:`isophotal <flux_iso_def>`, :ref:`corrected-isophotal <mag_isocor_def>`, :ref:`fixed-aperture <flux_aper_def>` and :ref:`adaptive-aperture <flux_auto_def>`.
For every :param:`FLUX_` measurement, an error estimate :param:`FLUXERR_`, a magnitude :param:`MAG_` and a magnitude error estimate :param:`MAGERR_` are also available.
The ``MAG_ZEROPOINT`` configuration parameter sets the magnitude zero-point
of magnitudes:

.. math::
 :label: mag

   {\tt MAG} = \mathrm{MAG\_ZEROPOINT} - 2.5 \log_{10} {\tt FLUX}

Magnitude uncertainties (error estimates) are computed using

.. math::
 :label: magerr

   {\tt MAGERR} = \frac{2.5}{\ln 10}\frac{\tt FLUXERR}{\tt FLUX}


.. _flux_iso_def:

Isophotal flux: :param:`FLUX_ISO`
---------------------------------

:param:`FLUX_ISO` is computed simply by integrating pixels values :math:`p_i` within the detection footprint, with the additional constraint that the background-subtracted, filtered value of detection image pixels must exceed the threshold set with the ``ANALYSIS_THRESH`` configuration parameter:
 
.. math::
 :label: fluxiso

   {\tt FLUX\_ISO} = F_{\rm iso} = \sum_{i \in {\cal D}} p_i.


.. _mag_isocor_def:

Corrected isophotal magnitude: :param:`MAG_ISOCOR`
--------------------------------------------------

.. note::
  Corrected isophotal magnitudes are now deprecated; they remain in |SExtractor| v2.x for compatibility with |SExtractor| v1.

:param:`MAG_ISOCOR` magnitudes are a quick-and-dirty way of retrieving the fraction of flux lost by isophotal magnitudes.

If one makes the assumption that the intensity profiles of faint objects recorded in the frame are roughly Gaussian because of atmospheric blurring, then the fraction :math:`\eta = \frac{F_{\rm iso}}{F_{\rm tot}}` of the total flux enclosed within a particular isophote reads :cite:`1990MNRAS_246_433M`:

.. math::
  :label: isocor

   \left(1-\frac{1}{\eta}\right ) \ln (1-\eta) = \frac{A\,t}{F_{\rm iso}},

where :math:`A` is the area and :math:`t` the threshold related to this isophote.
:eq:isocor is not analytically invertible, but a good approximation to :math:`\eta` (error :math:`< 10^{-2}` for :math:`\eta > 0.4`) can be done with the second-order polynomial fit:

.. math::
 :label: isocor2

   \eta \approx 1 - 0.1961 \frac{A\,t}{F_{\rm iso}} - 0.7512
   \left( \frac{A\,t}{F_{\rm iso}}\right)^2.

A “total” magnitude :param:`MAG_ISOCOR` estimate is then

.. math::
 :label: magisocor

   {\tt MAG\_ISOCOR} = {\tt MAG\_ISO} + 2.5 \log_{10} \eta.

Clearly this cheap correction works best with stars; and although it gives reasonably accurate results with most disk galaxies, it breaks down for ellipticals because of the broader wings in the profiles.

.. _flux_aper_def:

Fixed-aperture flux: :param:`FLUX_APER`
---------------------------------------

:param:`FLUX_APER` estimates the flux above the background within a circular aperture.
The diameter of the aperture in pixels is defined by the ``PHOTOM_APERTURES`` configuration parameter.
It does not have to be an integer: each "regular" pixel is subdivided in :math:`5\times 5` sub-pixels before measuring the flux within the aperture.
If :param:`FLUX_APER` is provided as a vector :param:`FLUX_APER[n]`, at least :math:`n` apertures must be specified with the ``PHOTOM_APERTURES`` configuration parameter.

.. _flux_auto_def:

Automatic aperture flux: :param:`FLUX_AUTO`
-------------------------------------------

:param:`FLUX_AUTO` provides an estimate of the “total flux” by integrating pixel values within an adaptively scaled aperture.
|SExtractor|’s automatic aperture photometry routine is inspired by Kron’s “first moment” algorithm :cite:`1980ApJS_43_305K`.

#. An elliptical aperture is :ref:`defined by the second order moments of the object’s light distribution <shape_iso_def>`, with semi-major axis :math:`a={\tt A\_IMAGE}`, semi-minor axis :math:`b={\tt B\_IMAGE}`, and position angle :param:`THETA_IMAGE`.
   The ellipse's major and minor axes are multiplied by 6 (which corresponds roughly to twice the size of the isophotal footprint on each axis).
#. Within this elliptical aperture :math:`{\cal E}`, Kron's "first moment" is computed:

.. math::
  :label: kron_r1

  r_1 = \frac{\sum_{i\in\cal E} d_i\,p_i}{\sum_{i\in\cal E} p_i},

where :math:`d_i` and :math:`p_i` are respectively the Euclidean distance to the (isophotal) centroid and the pixel value above the local background at image pixel :math:`i`.
:cite:`1980ApJS_43_305K` and :cite:`1987AA_183_177I` have shown that for stars and galaxy profiles convolved with Gaussian seeing, :math:`\ge 90\%` of the flux is expected to lie within a circular aperture of radius :math:`k r_1` if :math:`k = 2`, almost independently of the magnitude.
Experiments have shown :cite:`1996AAS_117_393B` that this conclusion remains unchanged if one replaces the circular aperture with a "Kron ellipse" :math:`{\cal K}`, which is simply the ellipse :math:`{\cal E}` above, scaled such that the principle axes become :math:`\sqrt{a/b}\,k r_1` and :math:`\sqrt{b/a}\,k r_1`.

:param:`FLUX_AUTO` is the flux above the background summed over the Kron ellipse:

.. math::
  :label: flux_auto

  {\tt FLUX\_AUTO} = \sum_{i\in\cal K} p_i.

The quantity :math:`kr_1` is known as the *Kron radius*, and is provided in |SExtractor| by the :param:`KRON_RADIUS` catalog parameter in units of :param:`A_IMAGE` (or :param:`B_IMAGE`):

.. math::
  :label: kron_radius

  {\tt KRON\_RADIUS} = \frac{1}{\sqrt{{\tt A\_IMAGE}\times{\tt B\_IMAGE}}}\,k\,r_1.


:math:`k = 2` defines a sort of balance between systematic and random errors.
By choosing a larger :math:`k = 2.5`, the mean fraction of flux lost drops from about 10% to 6%, at the expense of |SNR| in the measurement.
Very noisy objects may sometimes end up with a Kron radius being too small, even smaller that the isophotal footprint of the object itself. For this reason, |SExtractor| imposes a minimum size for the Kron radius, which must be larger than :math:`\sqrt{ab}\, r_{\rm min}`
The user has full control over the parameters :math:`k` and :math:`r_{\rm min}` through the ``PHOT_AUTOPARAMS`` configuration parameters. ``PHOT_AUTOPARAMS`` is set by default to ``2.5,3.5``.

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


Photographic photometry
-----------------------

In ``DETECT_TYPE PHOTO`` mode, SExtractor assumes that the response of the detector, over the dynamic range of the image, is logarithmic.
This is generally a good approximation for photographic density on deep exposures. 
Photometric procedures described above remain unchanged, except that for each pixel we apply first the transformation

.. math::

   I = I_0\,10^{D/\gamma} \ ,
   \label{eq:dtoi}

where :math:`\gamma` (``MAG_GAMMA``) is the contrast index of the emulsion, :math:`D` the original pixel value from the background-subtracted image, and :math:`I_0` is computed from the magnitude zero-point :math:`m_0`:

.. math:: I_0 = \frac{\gamma}{\ln 10} \,10^{-0.4\, m_0} \ .

One advantage of using a density-to-intensity transformation relative to the local sky background is that it corrects (to some extent) large-scale inhomogeneities in sensitivity (see :cite:`1996PhDT_68B` for details).


Magnitude uncertainties
-----------------------

An estimate of the error [#error]_ is available for each type of magnitude.
It is computed through

.. math:: \Delta m = 1.0857\, \frac{\sqrt{A\,\sigma^2 + F/g}}{F}

where :math:`A` is the area (in pixels) over which the total flux
:math:`F` (in ADU) is summed, :math:`\sigma` the standard deviation of
noise (in ADU) estimated from the background, and g the detector gain
(GAIN parameter [#gain]_ , in :math:`e^- / \mbox{ADU}`).
For corrected-isophotal magnitudes, a term, derived from Eq. [eq:isocor] is quadratically added to take into account the error on the correction itself.

In ``DETECT_TYPE PHOTO`` mode, things are slightly more complex.
Making the assumption that plate-noise is the major contributor to photometric errors, and that it is roughly constant in density, one can write:

.. math::

   \Delta m = 1.0857 \,\ln 10\, {\sigma\over \gamma}\,
     \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}
   =2.5\,{\sigma\over \gamma}\,
   \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}

where :math:`I(x,y)` is the contribution of pixel :math:`(x,y)` to the
total flux (Eq. [eq:dtoi]). ``GAIN`` is ignored in ``PHOTO`` mode.


..
  Background
  ----------

  Background is the last point relative to photometry. The assumption made in §[chap:backest] — that the "local" background associated to an object can be interpolated from the global background map — is no longer valid in crowded regions.
  An example is a globular cluster superimposed to a bulge of galaxy.
  SExtractor offers the possibility to estimate locally the background used to compute magnitudes.
  When this option is switched on (``BACKPHOTO_TYPE LOCAL`` instead of ``GLOBAL``), the "photometric"   background is estimated within a "rectangular annulus" around the isophotal limits of the object.
  The thickness of the annulus (in pixels) can be specified by the user with ``BACKPHOTO_SIZE``. A typical value is ``BACKPHOTO_SIZE``=``24``.

.. [#error]
   It is important to note that this error provides a lower limit, since
   it does not take into account the (complex) uncertainty on the local
   background estimate.

.. [#gain]
   Setting GAIN to 0 in the configuration file is equivalent to
   :math:`g = +\infty`

