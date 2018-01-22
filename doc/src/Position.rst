.. File Position.rst

.. include:: global.rst

Isophotal measurements
======================

Position and shape
------------------

The following quantities are derived from the spatial distribution :math:`\cal S` of pixels detected above the detection threshold (see :ref:`description <isophotal_measurements>`).

.. important::
  Unless otherwise noted, the pixel values used for computing "isophotal" positions and shapes are taken from the filtered, background-subtracted detection image.

.. note::
  Unless otherwise noted, all parameter names given below are only prefixes.
  They must be followed by _IMAGE if the results shall be expressed in pixel coordinates or :param:`_WORLD`, :param:`_SKY`, :param:`_J2000` or :param:`_B1950` for |WCS|_ coordinates (see :ref:`coord_suffix`).


.. _xyminmax_def:

Limits: :param:`XMIN`, :param:`YMIN`, :param:`XMAX`, :param:`YMAX`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These coordinates define two corners of a rectangle which encloses the detected object:

.. math::
  :label: xminymax

   \begin{eqnarray}
   {\tt XMIN} & = & \min_{i \in {\cal S}} x_i,\\
   {\tt YMIN} & = & \min_{i \in {\cal S}} y_i,\\
   {\tt XMAX} & = & \max_{i \in {\cal S}} x_i,\\
   {\tt YMAX} & = & \max_{i \in {\cal S}} y_i,
   \end{eqnarray}

where :math:`x_i` and :math:`y_i` are respectively the x-coordinate and y-coordinate of pixel :math:`i`.


.. _pos_iso_def:

Barycenter: :param:`X`, :param:`Y`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Barycenter coordinates generally define the position of the “center” of a source, although this definition can be inadequate or inaccurate if its spatial profile shows a strong skewness or very large wings. X and Y are simply computed as the first order moments of the profile:

.. math::
  :label: xy

   \begin{eqnarray}
   {\tt X} & = & \overline{x} = \frac{\displaystyle \sum_{i \in {\cal S}}
   p^{(f)}_i x_i}{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i},\\
   {\tt Y} & = & \overline{y} = \frac{\displaystyle \sum_{i \in {\cal S}}
   p^{(f)}_i y_i}{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i},
   \end{eqnarray}

where :math:`p^{(f)}_i` is the value of the pixel :math:`i` in the (filtered) detection image.
In practice, :math:`x_i` and :math:`y_i` are summed relative to :param:`XMIN` and :param:`YMIN` in order to reduce roundoff errors in the summing.


.. _pospeak_def:

Position of the peak: :param:`XPEAK`, :param:`YPEAK`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is sometimes useful to have the position :param:`XPEAK`, :param:`YPEAK` of the pixel with maximum intensity in a detected object, for instance when working with likelihood maps, or when searching for artifacts. For better robustness, PEAK coordinates are computed on *filtered* profiles if available. On symmetrical profiles, PEAK positions and barycenters coincide within a fraction of pixel (:param:`XPEAK` and :param:`YPEAK` coordinates are quantized by steps of 1 pixel, hence :param:`XPEAK_IMAGE` and :param:`YPEAK_IMAGE` are integers). This is no longer true for skewed profiles, therefore a simple comparison between PEAK and barycenter coordinates can be used to identify asymmetrical objects on well-sampled images.


.. _moments_iso_def:

Second-order moments: :param:`X2`, :param:`Y2`, :param:`XY`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Centered) second-order moments are convenient for measuring the spatial spread of a source profile. In |SExtractor| they are computed with:

.. math::
  :label: x2y2

   \begin{eqnarray}
   {\tt X2} & = \overline{x^2} = & \frac{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i x_i^2}{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i} - \overline{x}^2,\\
   {\tt Y2} & = \overline{y^2} = & \frac{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i y_i^2}{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i} - \overline{y}^2,\\
   {\tt XY} & = \overline{xy} = & \frac{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i x_i y_i}{\displaystyle \sum_{i \in {\cal S}} p^{(f)}_i} - \overline{x}\,\overline{y},
   \end{eqnarray}

These expressions are more subject to roundoff errors than if the 1st-order moments were subtracted before summing, but allow both 1st and
2nd order moments to be computed in one pass.
Roundoff errors are however kept to a negligible value by measuring all positions relative here again to :param:`XMIN` and :param:`YMIN`.


.. _shape_iso_def:

Basic shape parameters: :param:`A`, :param:`B`, :param:`THETA`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters are intended to describe the detected object as an elliptical shape. :param:`A` and :param:`B` are the lengths of the semi-major and semi-minor axes, respectively.
More precisely, they represent the maximum and minimum spatial dispersion of the object profile along any direction.
:param:`THETA` is the position-angle of the :param:`A` axis relative to the first image axis.
It is counted positive in the direction of the second axis.

Here is how shape parameters are computed. 2nd-order moments can easily be expressed in a referential rotated from the :math:`x,y` image coordinate system by an angle +\ :math:`\theta`:

.. math::
  :label: varproj

   \begin{array}{lcrrr}
   \overline{x_{\theta}^2} & = & \cos^2\theta\:\overline{x^2} & +\,\sin^2\theta\:\overline{y^2}
               & -\,2 \cos\theta \sin\theta\:\overline{xy},\\
   \overline{y_{\theta}^2} & = & \sin^2\theta\:\overline{x^2} & +\,\cos^2\theta\:\overline{y^2}
               & +\,2 \cos\theta \sin\theta\:\overline{xy},\\
   \overline{xy_{\theta}} & = & \cos\theta \sin\theta\:\overline{x^2} &
   -\,\cos\theta \sin\theta\:\overline{y^2} & +\,(\cos^2\theta -
   \sin^2\theta)\:\overline{xy}.
   \end{array}

One can find the angle(s) :math:`\theta_0` for which the variance is minimized (or maximized) along :math:`x_{\theta}`:

.. math::
  :label: theta0

  {\left.\frac{\partial \overline{x_{\theta}^2}}{\partial \theta} \right|}_{\theta_0} = 0,

which leads to

.. math::
  :label: theta0_2

   2 \cos\theta \sin\theta_0\:(\overline{y^2} - \overline{x^2})
       + 2 (\cos^2\theta_0 - \sin^2\theta_0)\:\overline{xy} = 0.

If :math:`\overline{y^2} \neq \overline{x^2}`, this implies:

.. math::
   :label: theta0_3

   \tan 2\theta_0 = 2 \frac{\overline{xy}}{\overline{x^2} - \overline{y^2}},

a result which can also be obtained by requiring the covariance :math:`\overline{xy_{\theta_0}}` to be null.
Over the domain :math:`[-\pi/2, +\pi/2[`, two different angles — with opposite signs — satisfy :eq:`theta0_3`.
By definition, :param:`THETA` is the position angle for which :math:`\overline{x_{\theta}^2}` is maximized.
:param:`THETA` is therefore the solution to :eq:`theta0_3` that has the same sign as the
covariance :math:`\overline{xy}`.
:param:`A` and :param:`B` can now simply be expressed as:

.. math::
  :label: aimage

   \begin{eqnarray}
   {\tt A}^2 & = & \overline{x^2}_{\tt THETA},\ \ \ {\rm and}\\
   {\tt B}^2 & = & \overline{y^2}_{\tt THETA}.
   \end{eqnarray}

:param:`A` and :param:`B` can be computed directly from the 2nd-order moments, using the following equations derived from :eq:`varproj` after some algebra:

.. math::
  :label: aimage_2

   \begin{eqnarray}
   {\tt A}^2 & = & \frac{\overline{x^2}+\overline{y^2}}{2}
       + \sqrt{\left(\frac{\overline{x^2}-\overline{y^2}}{2}\right)^2 + \overline{xy}^2},\\
   {\tt B}^2 & = & \frac{\overline{x^2}+\overline{y^2}}{2},
       - \sqrt{\left(\frac{\overline{x^2}-\overline{y^2}}{2}\right)^2 + \overline{xy}^2}.
   \end{eqnarray}

Note that :param:`A` and :param:`B` are exactly halves the :math:`a` and :math:`b` parameters computed by the COSMOS image analyser :cite:`1980SPIE_264_208S`.
Actually, :math:`a` and :math:`b` are defined in :cite:`1980SPIE_264_208S` as the semi-major and semi-minor axes of an elliptical shape with constant surface brightness, which would have the same 2nd-order moments as the analyzed object.


.. _elong_iso_def:

By-products of shape parameters: :param:`ELONGATION` and :param:`ELLIPTICITY`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters [#elongation]_ are directly derived from :param:`A` and :param:`B`:

.. math::
  :label: elongation

   \begin{eqnarray}
   {\tt ELONGATION} & = & \frac{\tt A}{\tt B}\ \ \ \ \ \mbox{and}\\
   {\tt ELLIPTICITY} & = & 1 - \frac{\tt B}{\tt A}.
   \end{eqnarray}

.. [#elongation] These parameters are dimensionless, and for historical reasons do not accept suffixes such as :param:`_IMAGE` or :param:`_WORLD`.

.. _ellipse_iso_def:

Ellipse parameters: :param:`CXX`, :param:`CYY`, :param:`CXY`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:param:`A`, :param:`B` and :param:`THETA` are not very convenient to use when, for instance, one wants to know if a particular |SExtractor| detection extends over some
position.
For this kind of application, three other ellipse parameters are provided; :param:`CXX`, :param:`CYY`, and :param:`CXY`.
They do nothing more than describing the same ellipse, but in a different way: the elliptical shape associated to a detection is now parameterized as

.. math::
  :label: ellipse

   {\tt CXX} (x-\overline{x})^2 + {\tt CYY} (y-\overline{y})^2
       + {\tt CXY} (x-\overline{x})(y-\overline{y}) = R^2,

where :math:`R` is a parameter which scales the ellipse, in units of :param:`A`
(or :param:`B`). Generally, the isophotal limit of a detected object is well
represented by :math:`R\approx 3` (:numref:`fig_ellipse`). Ellipse
parameters can be derived from the 2nd order moments:

.. math::
  :label: ellipse_2

   \begin{eqnarray}
   {\tt CXX} & = & \frac{\cos^2 {\tt THETA}}{{\tt A}^2} + \frac{\sin^2
   {\tt THETA}}{{\tt B}^2} =
   \frac{\overline{y^2}}{\overline{x^2} \overline{y^2} - \overline{xy}^2}\\
   {\tt CYY} & = & \frac{\sin^2 {\tt THETA}}{{\tt
   A}^2} + \frac{\cos^2 {\tt THETA}}{{\tt B}^2} =
   \frac{\overline{x^2}}{\overline{x^2} \overline{y^2} - \overline{xy}^2}\\
   {\tt CXY} & = & 2 \,\cos {\tt THETA}\,\sin {\tt
   THETA} \left( \frac{1}{{\tt A}^2} - \frac{1}{{\tt B}^2}\right) = -2\,
   \frac{\overline{xy}}{\overline{x^2} \overline{y^2} - \overline{xy}^2}.
   \end{eqnarray}

.. _fig_ellipse:

.. figure:: figures/ellipse.*
   :figwidth: 100%
   :align: center

   Meaning of basic shape parameters.


.. _poserr_iso_def:

Position uncertainties: :param:`ERRX2`, :param:`ERRY2`, :param:`ERRXY`, :param:`ERRA`, :param:`ERRB`, :param:`ERRTHETA`, :param:`ERRCXX`, :param:`ERRCYY`, :param:`ERRCXY`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uncertainties on the position of the barycenter can be estimated using
photon statistics.
In practice, such estimates are a lower-value of the full uncertainties because they do not include, for instance, the contribution of detection biases or contamination by neighbors.
Furthermore, |SExtractor| does not currently take into account possible correlations of the noise between adjacent pixels. Hence variances simply write:

.. math::
  :label: errxy

   \begin{eqnarray}
   {\tt ERRX2} & = {\rm var}(\overline{x}) =
   & \frac{\displaystyle \sum_{i \in {\cal S}} \sigma^2_i (x_i-\overline{x})^2}
   {\displaystyle \left(\sum_{i \in {\cal S}} p^{(f)}_i\right)^2},\\
   {\tt ERRY2} & = {\rm var}(\overline{y}) =
   & \frac{\displaystyle \sum_{i \in {\cal S}} \sigma^2_i (y_i-\overline{y})^2}
   {\displaystyle \left(\sum_{i \in {\cal S}} p^{(f)}_i\right)^2},\\
   {\tt ERRXY} & = {\rm cov}(\overline{x},\overline{y}) =
   & \frac{\displaystyle \sum_{i \in {\cal S}} \sigma^2_i (x_i-\overline{x})(y_i-\overline{y})}
   {\displaystyle \left(\sum_{i \in {\cal S}} p^{(f)}_i\right)^2}.
   \end{eqnarray}

:math:`\sigma_i` is the flux uncertainty estimated for pixel :math:`i`:

.. math:: \sigma^2_i = {\sigma_B}^2_i + \frac{p^{(f)}_i}{g_i},

where :math:`{\sigma_B}_i` is the local background noise and :math:`g_i` the local gain — conversion factor — for pixel :math:`i` (see :ref:`effect_of_weighting` for more details). Semi-major axis :param:`ERRA`, semi-minor axis :param:`ERRB`, and position angle :param:`ERRTHETA` of the :math:`1\sigma` position error ellipse are computed from the covariance matrix exactly like for :ref:`basic shape parameters<shape_iso_def>`:

.. math::
  :label: errabtheta

   \begin{eqnarray}
   {\tt ERRA}^2 & = & \frac{{\rm var}(\overline{x})+{\rm var}(\overline{y})}{2}
       + \sqrt{\left(\frac{{\rm var}(\overline{x})-{\rm var}(\overline{y})}{2}\right)^2
       + {\rm cov}^2(\overline{x},\overline{y})},\\
   {\tt ERRB}^2 & = & \frac{{\rm var}(\overline{x})+{\rm var}(\overline{y})}{2}
       - \sqrt{\left(\frac{{\rm var}(\overline{x})-{\rm var}(\overline{y})}{2}\right)^2
       + {\rm cov}^2(\overline{x},\overline{y})},\\
   \tan (2{\tt ERRTHETA}) & = & 2 \,\frac{{\rm cov}(\overline{x},\overline{y})}
                       {{\rm var}(\overline{x}) - {\rm var}(\overline{y})}.
   \end{eqnarray}

And the error ellipse parameters are:

.. math::
  :label: errellipse

   \begin{eqnarray}
   {\tt ERRCXX} & = & \frac{\cos^2 {\tt ERRTHETA}}{{\tt ERRA}^2} +
   \frac{\sin^2 {\tt ERRTHETA}}{{\tt ERRB}^2} = \frac{{\rm
   var}(\overline{y})}{{\rm var}(\overline{x}) {\rm var}(\overline{y}) -
   {\rm cov}^2(\overline{x},\overline{y})},\\
   {\tt ERRCYY} & = & \frac{\sin^2 {\tt ERRTHETA}}{{\tt ERRA}^2} +
   \frac{\cos^2 {\tt ERRTHETA}}{{\tt ERRB}^2} =
   \frac{{\rm var}(\overline{x})}{{\rm var}(\overline{x}) {\rm var}(\overline{y}) -
   {\rm cov}^2(\overline{x},\overline{y})},\\
   {\tt ERRCXY} & = & 2 \cos {\tt
   ERRTHETA}\sin {\tt ERRTHETA} \left( \frac{1}{{\tt ERRA}^2} -
   \frac{1}{{\tt ERRB}^2}\right)\\ & = & -2 \frac{{\rm
   cov}(\overline{x},\overline{y})}{{\rm var}(\overline{x}) {\rm var}(\overline{y}) -
   {\rm cov}^2(\overline{x},\overline{y})}.
   \end{eqnarray}


Handling of “infinitely thin” detections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from the mathematical singularities that can be found in some of the above equations describing shape parameters (and which |SExtractor| is able to handle, of course), some detections with very specific shapes may yield unphysical parameters, namely null values for :param:`B`, :param:`ERRB`, or even :param:`A` and :param:`ERRA`. 
Such detections include single-pixel objects and horizontal, vertical or diagonal lines which are 1-pixel wide. They will generally originate from glitches; but strongly undersampled and/or low S/N genuine sources may also produce such shapes.

For basic shape parameters, the following convention was adopted: if the
light distribution of the object falls on one single pixel, or lies on a
sufficiently thin line of pixels, which we translate mathematically by

.. math::
  :label: singutest

   \overline{x^2}\,\overline{y^2} - \overline{xy}^2 < \rho^2,

then :math:`\overline{x^2}` and :math:`\overline{y^2}` are incremented by :math:`\rho`. |SExtractor| sets :math:`\rho=1/12`, which is the variance of a 1-dimensional top-hat distribution with unit width.
Therefore :math:`1/\sqrt{12}` represents the typical minor-axis values assigned (in pixels units) to undersampled sources in |SExtractor|.

Positional errors are more difficult to handle, as objects with very high signal-to-noise can yield extremely small position uncertainties, just like singular profiles do. Therefore |SExtractor| first checks that :eq:`singutest` is true. If this is the case, a new test is conducted:

.. math::
  :label: singutest2

   {\rm var}(\overline{x})\,{\rm var}(\overline{y}) - {\rm
   covar}^2(\overline{x}, \overline{y}) < \rho^2_e,

where :math:`\rho_e` is arbitrarily set to :math:`\left( \sum_{i \in {\cal S}} \sigma^2_i \right) / \left(\sum_{i \in {\cal S}} p_i\right)^2`.
If :eq:`singutest2` is true, then :math:`\overline{x^2}` and :math:`\overline{y^2}` are incremented by :math:`\rho_e`.

.. _isoarea_def:

Isophotal areas: :param:`ISOAREA`, :param:`ISOAREAF`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An isophotal area is defined as the number of pixels with values exceeding some threshold above the background. Isophotal areas played an important role in the 80's and the 90's by providing,  at a small computing cost,  morphometric quantities complementary to second-order moments. |SExtractor| computes two isophotal areas inside the detection footprint:

* :param:`ISOAREAF` applies to the (filtered) detection image, above the detection threshold.
* :param:`ISOAREA` applies to the measurement image, with the additional constraint that the background-subtracted value of measurement image pixels must exceed the threshold set with the ``ANALYSIS_THRESH`` configuration parameter.

Photometry
----------

.. important::
  Unless otherwise noted, the pixel values used for computing "isophotal" fluxes and surface brightnesses are taken from the background-subtracted measurement image.

.. _flux_iso_def:

Isophotal flux: :param:`FLUX_ISO`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:param:`FLUX_ISO` is computed simply by integrating the background-subracted pixels values :math:`p_i` from the *measurement* image within the detection footprint
 
.. math::
 :label: fluxiso

   {\tt FLUX\_ISO} = F_{\rm iso} = \sum_{i \in {\cal S}} p_i.

|SExtractor| also provides an estimation of the uncertainty :param:`FLUXERR_ISO`, a magnitude :param:`MAG_ISO` and a magnitude error estimate :param:`MAGERR_ISO`: see :ref:`fluxes_and_magnitudes`.

.. _mag_isocor_def:

Corrected isophotal magnitude: :param:`MAG_ISOCOR`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Clearly the :param:`MAG_ISOCOR` correction works best with stars; and although it gives reasonably accurate results with most disk galaxies, it breaks down for ellipticals because of the broader wings in the profiles.

.. _flux_max_def:

Peak value: :param:`FLUX_MAX`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:param:`FLUX_MAX` is the peak pixel value (above the local background) in the measurement image. Note that it may not correspond to the pixel with coordinates given by :param:`XPEAK` and :param:`YPEAK` (see :ref:`Position of the peak <pospeak_def>`) if ``FILTER`` is set to ``Y`` or if the measurement image differs from the detection image.

.. _class_star_def:

:param:`CLASS_STAR` classifier
------------------------------

.. note::
  The :param:`CLASS_STAR` classifier has been superseded by the :param:`SPREAD_MODEL` estimator (see :ref:`spread_model_def`), which offers better performance by making explicit use of the full, variable |PSF| model.

A good discrimination between stars and galaxies is essential for both galactic and extragalactic statistical studies.
The common assumption is that galaxy images look more extended or fuzzier than those of stars (or |QSO|\ s).
|SExtractor| provides the :param:`CLASS_STAR` catalog parameter for separating both types of sources.
The  :param:`CLASS_STAR` classifier relies on a `multilayer feed-forward neural network <https://en.wikipedia.org/wiki/Multilayer_perceptron>`_ trained using `supervised learning <https://en.wikipedia.org/wiki/Supervised_learning>`_ to estimate the *a posteriori* probability :cite:`Richard1991,Saerens2002` of a |SExtractor| detection to be a point source or an extended object.
Below is a shortened description of the estimator, see :cite:`1996AAS_117_393B` for more details.

Inputs and outputs
~~~~~~~~~~~~~~~~~~

The neural network is a multilayer Perceptron with a single fully connected, hidden layers.
Of all neural networks it is probably the best-studied, and it has been intensively applied with success for many classification tasks.

The classifier (:numref:`fig_classstarnn`) has 10 inputs:

* 8 isophotal areas :math:`A_0..A_7`, measured at isophotes exponentially spaced between the analysis threshold (which may be modified with the ``ANALYSIS_THRESH`` configuration parameter) and the object's peak pixel value
* The object's peak pixel value above the local background :math:`I_{\mathrm max}`
* A |seeing| input, which acts as a tuning button.

The output layer contains only one neuron, as "star" and "galaxy" are two classes mutually exclusive.
The output value is a "stellarity index", which for images that reasonably match those of the training sample is an estimation of the *a posteriori* probability for the classified object to be a point-source.
Hence a :param:`CLASS_STAR` close to 0 means that the object is very likely a galaxy, and 1 that it is a star.
In practice, real data always differ at least slightly from the training sample, and the :param:`CLASS_STAR` output is often a poor approximation of the expected *a posteriori* probabilities.
Nevertheless, a :param:`CLASS_STAR` value closer to 0 or 1 normally indicates a higher confidence in the classification, and the balance between sample completeness and purity may still be adjusted by tweaking the decision threshold .

.. _fig_classstarnn:

.. figure:: figures/classstarnn.*
  :figwidth: 100%
  :align: center

  Architecture of |SExtractor|'s :param:`CLASS_STAR` classifier

The |seeing| input must be set by the user with the ``SEEING_FWHM`` configuration parameter.
If ``SEEING_FWHM`` is set to 0,  it is automatically measured on the |PSF| model which must be provided (using the ``PSF_NAME`` configuration parameter).

If no |PSF| model is available, the ``SEEING_FWHM`` configuration parameter must be adjusted by the user to match the actual average |PSF| |FWHM| on the image.
The accuracy with which ``SEEING_FWHM`` must be set for optimal results ranges from :math:`\pm 20\%` for bright sources to about :math:`\pm 5\%` for the faintest (:numref:`fig_classstar_seeing`). ``SEEING_FWHM`` is expressed in arcseconds.
The ``PIXEL_SCALE`` configuration parameter must therefore also be set by the user if |WCS| information is missing from the |FITS| image header.
There are several ways to measure, directly or indirectly, the size of point sources in |SExtractor|; they may lead to slightly discordant results, depending on the exact shape of the |PSF|.
The measurement :param:`FWHM_IMAGE` (although not the most reliable as an image quality estimator) sets the reference when it comes to setting ``SEEING_FWHM``.

One may check that the ``SEEING_FWHM`` is set correctly by making sure that the typical :param:`CLASS_STAR` value of unclassifiable sources at the faint end of the catalog hovers around the 0.5 mark.

.. _fig_classstar_seeing:

.. figure:: figures/classstar_seeing.*
  :figwidth: 100%
  :align: center

  Architecture of |SExtractor|'s :param:`CLASS_STAR` classifier


Training
~~~~~~~~

This section gives some insight on how the :param:`CLASS_STAR` classifier has been trained.
The main issue with supervised machine learning is the labeling of the large training sample.
Hopefully a big percentage of contemporary astronomical images share a common set of generic features, which can be simulated with sufficient realism to create a large training sample together with the ground truth (labels).
The :param:`CLASS_STAR` classifier was trained on such a sample of artificial images.

Six hundred :math:`512\times512` simulation images containing stars and galaxies were generated to train the network using an early prototype of the |SkyMaker|_ package :cite:`2009MmSAI_80_422B`.
They were done in the blue band, where galaxies present very diversified aspects.
The two parameters governing the shape of the |PSF| (|seeing| |FWHM| and Moffat :math:`\beta` parameter :cite:`1969AA_3_455M`) were chosen randomly with :math:`0.025\le` FWHM :math:`\le 5.5` and :math:`2\le\beta\le4`. Note that the `Moffat function <https://en.wikipedia.org/wiki/Moffat_distribution>`_ used in the simulation is a poor approximation to diffraction-limited images, hence the :param:`CLASS_STAR` classifier is not optimized for space-based images.
The pixel scale was always taken less than :math:`\approx 0.7` FWHM to warrant correct sampling of the image.
Bright galaxies are simply too rare in the sky to constitute a significant training sample on such a small number of simulations.
So, keeping a constant comoving number density, we increased artificially the number of nearby galaxies by making the volume element proportional to :math:`zdz`.
Stars were given a number-magnitude distribution identical to that of galaxies.
**Therefore any pattern presented to the network had a 50% chance to correspond to a star or a galaxy, irrespective of magnitude** [#faint]_.
Crowding in the simulated images was higher than what one sees on real images of high galactic latitude fields, allowing for the presence of many “difficult cases” (close double stars, truncated profiles, etc...) that the neural network classifier had to deal with.

|SExtractor| was run on each image with 8 different extraction thresholds. This led to a catalog with about :math:`10^6` entries with the 10 input parameters plus the class label. Back-propagation learning took about 15 min on a SUN SPARC20 workstation. The final set of synaptic weights was saved to the file :file:`default.nnw` , ready to be used in “feed-forward” mode during source extraction.

.. [#faint]
   Faint galaxies have less chance being detected than faint stars, but it has little effect because the ones that are lost at a given magnitude are predominantly the most extended and consequently the easiest to classify.

