.. File Photom.rst

Photometry
==========

Besides PSF and model-fitting flux estimates, |SExtractor| can currently
perform four types of flux measurements: isophotal, *corrected-isophotal*,
fixed-aperture and *adaptive-aperture*. For every ``FLUX`` measurement,
an error estimate ``FLUXERR``, a magnitude ``MAG`` and a magnitude
error estimate ``MAGERR`` are also available.
The ``MAG_ZEROPOINT`` configuration parameter sets the magnitude zero-point
of magnitudes:

.. math::
 :label: mag

   {\tt MAG} = \mathrm{MAG\_ZEROPOINT} - 2.5 \log_{10} {\tt FLUX}

Magnitude uncertainties (error estimates) are computed using

.. math::
 :label: magerr

   {\tt MAGERR} = \frac{2.5}{\ln 10}\frac{\tt FLUXERR}{\tt FLUX}

Isophotal flux
--------------

``FLUX_ISO`` is computed simply by integrating pixels values within the
detection footprint, with the additional constraint that the
background-subtracted, filtered value of detection image pixels must exceed
the threshold set with the ``ANALYSIS_THRESH`` configuration parameter.
 
.. math::
 :label: fluxiso

   {\tt FLUX\_ISO} = \sum_{i \in {\cal S}} I_i


Corrected isophotal magnitudes
------------------------------

``MAG_ISOCOR`` can be considered as a quick-and-dirty way for retrieving
the fraction of flux lost by isophotal magnitudes. Although their use is
now deprecated, they have been kept in |SExtractor| v2.x and above for
compatibility with |SExtractor| v1. If we make the assumption that the
intensity profiles of the faint objects recorded in the frame are
roughly Gaussian because of atmospheric blurring, then the fraction
:math:`\eta = \frac{I_{\rm iso}}{I_{\rm tot}}` of the total flux enclosed
within a particular isophote reads :cite:`1990MNRAS.246..433M`:

.. math::
  :label: isocor

   \left(1-\frac{1}{\eta}\right ) \ln (1-\eta) = \frac{A\,t}{I_{\rm iso}}

where :math:`A` is the area and :math:`t` the threshold related to this
isophote. :eq:isocor is not analytically invertible, but a good
approximation to :math:`\eta` (error :math:`< 10^{-2}` for :math:`\eta >
0.4`) can be done with the second-order polynomial fit:

.. math::
 :label: isocor2

   \eta \approx 1 - 0.1961 \frac{A\,t}{I_{\rm iso}} - 0.7512
   \left( \frac{A\,t}{I_{\rm iso}}\right)^2 \label{eq:isocor}

A “total” magnitude ``MAG_ISOCOR`` estimate is then

.. math::
 :label: magisocor

   {\tt MAG\_ISOCOR} = {\tt MAG\_ISO} + 2.5 \log_{10} \eta

Clearly this cheap correction works best with stars; and although it is
shown to give tolerably accurate results with most disk galaxies, it
fails with ellipticals because of the broader wings of their profiles.

Fixed-aperture flux
-------------------

``FLUX_APER`` estimates the flux above the background within a circular
aperture. The diameter of the aperture in pixels is defined by the
``PHOTOM_APERTURES`` configuration parameter. It does not have to be an
integer: each “normal” pixel is subdivided in :math:`5\times 5` sub-pixels
before measuring the flux within the aperture. If ``FLUX_APER`` is provided as a
vector ``FLUX_APER[n]``, at least :math:`n` apertures must be
specified with ``PHOTOM_APERTURES``.

Automatic aperture magnitudes
-----------------------------

(MAG\_AUTO) provides an estimate of the “total magnitude” by integrating
the source flux within an adaptively scaled aperture. SExtractor’s
automatic aperture photometry routine is inspired by Kron’s “first
moment” algorithm (1980). (1) We define an elliptical aperture whose
elongation :math:`\epsilon` and position angle :math:`\theta` are
defined by second order moments of the object’s light distribution. The
ellipse is scaled to :math:`R_{\rm max}.\sigma_{\rm iso}`
(:math:`6 \sigma_{\rm iso}`, which corresponds roughly to 2 isophotal
“radii”). (2) Within this aperture we compute the “first moment”:

.. math:: r_1 = \frac{\sum r\,I(r)}{\sum I(r)}

Kron (1980) and Infante (1987) have shown that for stars and galaxy
profiles convolved with Gaussian seeing, :math:`\ge 90\%` of the flux is
expected to lie within a circular aperture of radius :math:`k r_1` if
:math:`k =
2`, almost independently of their magnitude. This picture remains
unchanged if we consider an ellipse with :math:`\epsilon\, k r_1` and
:math:`k r_1 /
\epsilon` as principal axes. :math:`k = 2` defines a sort of balance
between systematic and random errors. By choosing a larger
:math:`k = 2.5`, the mean fraction of flux lost drops from about 10% to
6%. When Signal to Noise is low, it may appear that an erroneously small
aperture is taken by the algorithm. That’s why we have to bound the
smallest accessible aperture to :math:`R_{\rm min}` (typically
:math:`R_{\rm min} = 3 - 4\,
\sigma_{\rm iso}`). The user has full control over the parameters
:math:`k` and :math:`R_{\rm min}` through the configuration parameters
PHOT\_AUTOPARAMS; by default, PHOT\_AUTOPARAMS is set to 2.5,3.5.

.. figure:: ps/simlostflux.ps
   :alt:  Flux lost (expressed as a mean magnitude difference) with
   different faint-object photometry techniques as a function of total
   magnitude (see text). Only isolated galaxies (no blends) of the
   simulations have been considered.
   :width: 15.00000cm

    Flux lost (expressed as a mean magnitude difference) with different
   faint-object photometry techniques as a function of total magnitude
   (see text). Only isolated galaxies (no blends) of the simulations
   have been considered. 

Aperture magnitudes are sensitive to crowding. In SExtractor 1,
MAG\_AUTO measurements were not very robust in that respect. It was
therefore suggested to replace the aperture magnitude by the
corrected-isophotal one when an object is too close to its neighbours (2
isophotal radii for instance). This was done automatically when using
the MAG\_BEST magnitude: :math:`{\tt MAG\_BEST} = {\tt MAG\_AUTO}` when
it is sure that no neighbour can bias MAG\_AUTO by more than 10%, or
:math:`{\tt MAG\_BEST} = {\tt MAG\_ISOCOR}` otherwise. Experience showed
that the MAG\_ISOCOR and MAG\_AUTO magnitude would loose about the same
fraction of flux on stars or compact galaxy profiles: around 0.06 % for
default extraction parameters. The use of MAG\_BEST is now deprecated as
MAG\_AUTO measurements are much more robust in versions 2.x of
SExtractor. The first improvement is a crude subtraction of all the
neighbours which have been detected around the measured source (the
MASK\_TYPE BLANK option). The second improvement is an automatic
correction of parts of the aperture that are suspected to be
contaminated by a neighbour. This is done by mirroring the opposite,
cleaner side of the measurement ellipse if available (the MASK\_TYPE
CORRECT option, which is also the default). Figure [figphot] shows the
mean loss of flux measured with isophotal (threshold
:math:`= 24.4\ \mbox{\rm magnitude\,arsec}^{-2}`), corrected isophotal
and automatic aperture photometries for simulated galaxy :math:`B_J` on
a typical Schmidt-survey plate image. The automatic adaptive aperture
photometry leads to the lowest loss of flux.

Photographic photometry
-----------------------

In DETECT\_TYPE PHOTO mode, SExtractor assumes that the response of the
detector, over the dynamic range of the image, is logarithmic. This is
generally a good approximation for photographic density on deep
exposures. Photometric procedures described above remain unchanged,
except that for each pixel we apply first the transformation

.. math::

   I = I_0\,10^{D/\gamma} \ ,
   \label{eq:dtoi}

where :math:`\gamma` (MAG\_GAMMA) is the contrast index of the
emulsion, :math:`D` the original pixel value from the
background-subtracted image, and :math:`I_0` is computed from the
magnitude zero-point :math:`m_0`:

.. math:: I_0 = \frac{\gamma}{\ln 10} \,10^{-0.4\, m_0} \ .

 One advantage of using a density-to-intensity transformation relative
to the local sky background is that it corrects (to some extent)
large-scale inhomogeneities in sensitivity (see Bertin 1996 for
details).

Errors on magnitude
-------------------

An estimate of the error [1]_ is available for each type of magnitude.
It is computed through

.. math:: \Delta m = 1.0857\, \frac{\sqrt{A\,\sigma^2 + F/g}}{F}

where :math:`A` is the area (in pixels) over which the total flux
:math:`F` (in ADU) is summed, :math:`\sigma` the standard deviation of
noise (in ADU) estimated from the background, and g the detector gain
(GAIN parameter [2]_ , in :math:`e^- / \mbox{ADU}`). For
corrected-isophotal magnitudes, a term, derived from Eq. [eq:isocor] is
quadratically added to take into account the error on the correction
itself.

In DETECT\_TYPE PHOTO mode, things are slightly more complex. Making the
assumption that plate-noise is the major contributor to photometric
errors, and that it is roughly constant in density, we can write:

.. math::

   \Delta m = 1.0857 \,\ln 10\, {\sigma\over \gamma}\,
     \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}
   =2.5\,{\sigma\over \gamma}\,
   \frac{\sqrt{\sum_{x,y}{I^2(x,y)}}}{\sum_{x,y}I(x,y)}

where :math:`I(x,y)` is the contribution of pixel :math:`(x,y)` to the
total flux (Eq. [eq:dtoi]). The GAIN is ignored in PHOTO mode.

Background
----------

is the last point relative to photometry. The assumption made in
§[chap:backest] — that the “local” background associated to an object
can be interpolated from the global background map — is no longer valid
in crowded regions. An example is a globular cluster superimposed on a
bulge of galaxy. SExtractor offers the possibility to estimate locally
the background used to compute magnitudes. When this option is switched
on (BACKPHOTO\_TYPE LOCAL instead of GLOBAL), the “photometric”
background is estimated within a “rectangular annulus” around the
isophotal limits of the object. The thickness of the annulus (in pixels)
can be specified by the user with BACKPHOTO\_SIZE. A typical value is
BACKPHOTO\_SIZE=24.

.. [1]
   It is important to note that this error provides a lower limit, since
   it does not take into account the (complex) uncertainty on the local
   background estimate.

.. [2]
   Setting GAIN to 0 in the configuration file is equivalent to
   :math:`g = +\infty`

.. include:: keys.rst

