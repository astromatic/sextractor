.. File Model.rst

.. include:: global.rst

Model fitting
=============

Fitting procedure
-----------------

SExtractor can fit models to the images of detected objects since version 2.8. The fit is performed by minimizing the loss function

.. math::
  :label: loss_func

  \lambda(\boldsymbol{q}) = \sum_i \left(g\left(\frac{p_i - \tilde{m}_i(\boldsymbol{q})}{\sigma_i}\right)\right)^2 + \sum_j \left(\frac{f_j(Q_j) - \mu_j}{s_j}\right)^2

with respect to components of the model parameter vector :math:`\boldsymbol{q}`. :math:`\boldsymbol{q}` comprises parameters describing the shape of the model and the model pixel coordinates :math:`\boldsymbol{x}`.

Modified least squares
~~~~~~~~~~~~~~~~~~~~~~

The first term in :eq:`loss_func` is a modified `weighted sum of squares <http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>`_ that aims at minimizing the residuals of the fit. :math:`p_i`, :math:`\tilde{m}_i(\boldsymbol{q})` and :math:`\sigma_i` are respectively the pixel value above the background, the value of the resampled model, and the pixel value uncertainty at image pixel :math:`i`.
:math:`g(u)` is a derivable monotonous function that reduces the influence of large deviations from the model, such as the contamination by neighbors (:numref:`fig_robustgalfit`):

.. math::
  :label: modified_lsq

  g(u) = \left\{
    \begin{array}{rl}
       u_0 \log \left(1 + \frac{u}{u_0}\right) & \mbox{if } u \ge 0,\\
      -u_0 \log \left(1 - \frac{u}{u_0}\right) & \mbox{otherwise}.\\
    \end{array}
  \right.

:math:`u_0` sets the level below which :math:`g(u)\approx u`.
In practice, choosing :math:`u_0 = \kappa \sigma_i` with :math:`\kappa = 10` makes the first term in :eq:`loss_func` behave like a traditional weighted sum of squares for residuals close to the noise level.

.. caution::
  The cost function in :eq:`loss_func` is optimized for Gaussian noise and makes model-fitting in |SExtractor| appropriate only for image noise with a |pdf| symmetrical around the mean.

.. _fig_robustgalfit:

.. figure:: figures/robustgalfit.*
   :figwidth: 100%
   :align: center

   Effect of the modified least squares loss function on fitting a model to a galaxy with a bright neighbor. *Left*: the original image; *Middle*: residuals of the model fitting with a regular least squares (:math:`\kappa = +\infty`); *Right*: modified least squares with :math:`\kappa = 10`.

The vector :math:`\tilde{\boldsymbol{m}}(\boldsymbol{q})` is obtained by convolving the high resolution model :math:`\boldsymbol{m}(\boldsymbol{q})` with the local PSF model :math:`\boldsymbol{\phi}` and applying a resampling operator :math:`\mathbf{R}(\boldsymbol{x})` to generate the final model raster at position :math:`\boldsymbol{x}` at the nominal image resolution:

.. math::
  :label: model_convolution

  \tilde{\boldsymbol{m}}(\boldsymbol{q}) = \mathbf{R}(\boldsymbol{x}) (\boldsymbol{m}(\boldsymbol{q})*\boldsymbol{\phi}).

:math:`\mathbf{R}(\boldsymbol{x})` depends on the pixel coordinates :math:`\boldsymbol{x}` of the model centroid:

.. math::
  :label: model_resampling

  \mathbf{R}_{ij}(\boldsymbol{x}) =  h\left(\boldsymbol{x}_j - \eta.(\boldsymbol{x}_i - \boldsymbol{x})\right),

where :math:`h` is a 2-dimensional interpolant (interpolating function), :math:`\boldsymbol{x}_i` is the coordinate vector of image pixel :math:`i`, :math:`\boldsymbol{x}_j` the coordinate vector of model sample :math:`j`, and :math:`\eta` is the image-to-model sampling step ratio (sampling factor) which is by default defined by the PSF model sampling.
We adopt a Lánczos-4 function :cite:`duchon1979` as interpolant.

Change of variables
~~~~~~~~~~~~~~~~~~~

Many model parameters are valid only over a restricted domain.
Fluxes, for instance, cannot be negative. 
In order to avoid invalid values and also to facilitate convergence, a change of variables is applied individually to each model parameter:

.. math::
  :label: change_of_variables

  q_j = f_j(a_j, b_j, Q_j).

The "model" variable :math:`q_j` is bounded by the lower limit :math:`a_j` and the upper limit :math:`b_j` by construction.
The "engine" variable :math:`Q_j` can take any value, and is actually the parameter that is being adjusted in the fit, although it does not have any physical meaning.
In |SExtractor| three different types of transforms :math:`f_j()` are applied, depending on the parameter (:numref:`change_of_variable_table`).

.. _change_of_variable_table:

.. list-table:: Types of changes of variables applied to model parameters
  :header-rows: 1

  * - Type
    - Model :math:`\stackrel{f^{-1}}{\to}` Engine
    - Engine :math:`\stackrel{f}{\to}` Model
    - Examples
  * - Unbounded (linear)
    - :math:`Q_j = q_j`
    - :math:`q_j = Q_j`
    - | :param:`SPHEROID_POSANGLE`
      | :param:`DISK_POSANGLE`
  * - Bounded linear
    - :math:`Q_j = \ln \frac{q_j - a_j}{b_j - q_j}`
    - :math:`q_j = \frac{b_j - a_j}{1 + \exp -Q_j} + a_j`
    - | :param:`XMODEL_IMAGE`
      | :param:`SPHEROID_SERSICN`
  * - Bounded logarithmic
    - :math:`Q_j = \ln \frac{\ln q_j - \ln a_j}{\ln b_j - \ln q_j}`
    - :math:`q_j = a_j \frac{\ln b_j - \ln a_j}{1 + \exp -Q_j}`
    - | :param:`FLUX_SPHEROID`
      | :param:`DISK_ASPECT`

In practice, this approach works well, and was found to be much more reliable than a box constrained algorithm :cite:`Kanzow2004375`.

Regularization
~~~~~~~~~~~~~~

Although minimizing the (modified) weighted sum of least squares gives a solution that fits best the data, it does not necessarily correspond to the most probable solution given what we know about celestial objects.
The discrepancy is particularly significant in very faint (|SNR| :math:`\le 20`) and barely resolved galaxies, for which there is a tendency to overestimate the elongation, known as the "noise bias" in the weak-lensing community :cite:`2004MNRAS_353_529H,2012MNRAS_424_2757M,2012MNRAS_425_1951R,2012MNRAS_427_2711K`.
To mitigate this issue, |SExtractor| implements a simple `Tikhonov regularization <https://en.wikipedia.org/wiki/Tikhonov_regularization>`_ scheme on selected engine parameters, in the form of an additional penalty term in :eq:`loss_func`.
This term acts as a Gaussian prior on the selected *engine* parameters. However for the associated *model* parameters, the change of variables can make the (improper) prior far from Gaussian.
Currently the only regularized parameter is :param:`SPHEROID_ASPECT_IMAGE` (and its derivatives :param:`SPHEROID_ASPECT_WORLD`, :param:`ELLIP1MODEL_IMAGE`, etc.), for which :math:`\mu_{\tt SPHEROID\_ASPECT} = 0` and :math:`s_{\tt SPHEROID\_ASPECT} = 1`
(:numref:`fig_aspectprior`).

.. _fig_aspectprior:

.. figure:: figures/aspectprior.*
   :figwidth: 100%
   :align: center

   Effect of the Gaussian prior on the :param:`SPHEROID_ASPECT_IMAGE` model parameter. *Left:* change of variables between the model (in abscissa) and the engine (in ordinate) parameters. *Right*: equivalent (improper) prior applied to :param:`SPHEROID_ASPECT_IMAGE` for :math:`\mu_{\tt SPHEROID\_ASPECT} = 0` and :math:`s_{\tt SPHEROID\_ASPECT} = 1` in equation :eq:`loss_func`.

.. _model_minimization_def:

Minimization
~~~~~~~~~~~~

Minimization of the loss function :math:`\lambda(\boldsymbol{q})` is carried out using the `Levenberg-Marquardt algorithm <http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm>`_, and more specifically the |LevMar|_ implementation :cite:`lourakis04LM`.
The library approximates the Jacobian matrix of the model from finite differences using Broyden's :cite:`Broyden1965ACo` rank one updates.
The fit is done inside a disk which diameter is scaled to include the isophotal footprint of the object, plus the FWHM of the PSF, plus a 20 % margin.

The number of iterations is returned in the :param:`NITER_MODEL` measurement parameter.
It is generally a few tens.
The final value of the modified chi square term in :eq:`loss_func`, divided by the number of degrees of freedom, is returned in :param:`CHI2_MODEL`.

.. _flags_model_def:

The :param:`FLAGS_MODEL` catalog parameter flags various issues which may happen during the fitting process (see the :ref:`flagging` section for details on how flags are managed in |SExtractor|):

.. _flags_model_table:

.. csv-table:: :param:`FLAGS_MODEL` flag description
  :header: "Value", "Meaning"
  :widths: 3 60

  1, "the unconvolved, supersampled model raster exceeds 512×512 pixels and had to be resized"
  2, "the convolved, resampled model raster exceeds 512×512 pixels and had to be resized"
  4, "not enough pixels are available for model fitting on the measurement image (less pixels than fit parameters)"
  8, "at least one of the fitted parameters hits the lower bound"
  16, "at least one of the fitted parameters hits the upper bound"

:math:`1\,\sigma` error estimates are provided for most measurement parameters; they are obtained from the full covariance matrix of the fit, which is itself computed by inverting the approximate `Hessian matrix <https://en.wikipedia.org/wiki/Hessian_matrix>`_ of :math:`\lambda(\boldsymbol{q})` at the solution.

.. _models_def:

Models
------

Models contain one or more components, which share their central coordinates. For instance, a galaxy model may be composed of a spheroid (bulge) and a disk components. Both components are concentric but they may have different scales, aspect ratios and position angles. Adding a component is done simply by invoking one of its measurement parameters in the parameter file, e.g., :param:`DISK_SCALE_IMAGE`.

The present version of |SExtractor| supports the following models

- :param:`BACKOFFSET`: flat background offset

  Relevant measurement parameters: :param:`FLUX_BACKOFFSET`, :param:`FLUXERR_BACKOFFSET`

.. math::
  :label: backoffset_model

  m_{\tt BACKOFFSET}(r) = m_0


- :param:`POINT_SOURCE`: point source

  Relevant measurement parameters: :param:`FLUX_POINTSOURCE`, :param:`FLUXERR_POINTSOURCE`, :param:`MAG_POINTSOURCE`, :param:`MAGERR_POINTSOURCE`, :param:`FLUXRATIO_POINTSOURCE`, :param:`FLUXRATIOERR_POINTSOURCE`

.. math::
  :label: pointsource_model

  m_{\tt POINTSOURCE}(r) = m_0 \delta(r)

- :param:`DISK`: exponential disk

  Relevant measurement  parameters:
  :param:`FLUX_DISK`, :param:`FLUXERR_DISK`, :param:`MAG_DISK`, :param:`MAGERR_DISK`,
  :param:`FLUXRATIO_DISK`, :param:`FLUXRATIOERR_DISK`,
  :param:`FLUX_MAX_DISK`, :param:`MU_MAX_DISK`,
  :param:`FLUX_EFF_DISK`, :param:`MU_EFF_DISK`,
  :param:`FLUX_MEAN_DISK`, :param:`MU_MEAN_DISK`,
  :param:`DISK_SCALE_IMAGE`, :param:`DISK_SCALEERR_IMAGE`,
  :param:`DISK_SCALE_WORLD`, :param:`DISK_SCALEERR_WORLD`,
  :param:`DISK_ASPECT_IMAGE`, :param:`DISK_ASPECTERR_IMAGE`,
  :param:`DISK_ASPECT_WORLD`, :param:`DISK_ASPECTERR_WORLD`,
  :param:`DISK_INCLINATION`, :param:`DISK_INCLINATIONERR`,
  :param:`DISK_THETA_IMAGE`, :param:`DISK_THETAERR_IMAGE`,
  :param:`DISK_THETA_WORLD`, :param:`DISK_THETAERR_WORLD`,
  :param:`DISK_THETA_SKY`, :param:`DISK_THETA_J2000`, :param:`DISK_THETA_B1950`

.. math::
  :label: disk_model

  m_{\tt DISK}(r) = m_0 \exp \left( - {r\over h}\right) 

- :param:`SPHEROID`: Sérsic (:math:`R^{1/n}`) spheroid

  :param:`FLUX_SPHEROID`, :param:`FLUXERR_SPHEROID`, :param:`MAG_SPHEROID`, :param:`MAGERR_SPHEROID`,
  :param:`FLUXRATIO_SPHEROID`, :param:`FLUXRATIOERR_SPHEROID`,
  :param:`FLUX_MAX_SPHEROID`, :param:`MU_MAX_SPHEROID`,
  :param:`FLUX_EFF_SPHEROID`, :param:`MU_EFF_SPHEROID`,
  :param:`FLUX_MEAN_SPHEROID`, :param:`MU_MEAN_SPHEROID`,
  :param:`SPHEROID_SCALE_IMAGE`, :param:`SPHEROID_SCALEERR_IMAGE`,
  :param:`SPHEROID_SCALE_WORLD`, :param:`SPHEROID_SCALEERR_WORLD`,
  :param:`SPHEROID_ASPECT_IMAGE`, :param:`SPHEROID_ASPECTERR_IMAGE`,
  :param:`SPHEROID_ASPECT_WORLD`, :param:`SPHEROID_ASPECTERR_WORLD`,
  :param:`SPHEROID_INCLINATION`, :param:`SPHEROID_INCLINATIONERR`,
  :param:`SPHEROID_THETA_IMAGE`, :param:`SPHEROID_THETAERR_IMAGE`,
  :param:`SPHEROID_THETA_WORLD`, :param:`SPHEROID_THETAERR_WORLD`,
  :param:`SPHEROID_THETA_SKY`, :param:`SPHEROID_THETA_J2000`, :param:`SPHEROID_THETA_B1950`
  :param:`SPHEROID_SERSICN`, :param:`SPHEROID_SERSICNERR`

.. math::
  :label: spheroid_model

  m_{\tt SPHEROID}(r) = m_0 \exp \left(- b(n)\,\left({R\over R_e}\right)^{1/n}\right),

where, for the :cite:`1968adga_book_S` model, :math:`b(n)` is the solution to

.. math::
  :label: bofn

  2 \gamma[2\,n,b(n)] = \Gamma(2\,n)

An accurate approximation for the solution for :math:`b(n)` of :eq:`bofn` is :cite:`1999AA_352_447C`:

.. math::

  b(n) = 2\,n - {1\over3} + {4\over 405\,n} + {46\over 25515\,n^2} + {131\over 1148175\,n^3}

Experience shows that the de Vaucouleurs spheroid + exponential disk
combination provides fairly accurate and robust fits for moderately
resolved faint galaxies. An adjustable Sérsic index may offer lower
residuals on spheroids and/or well-resolved galaxies, but makes the fit
less robust and more sensitive to PSF model errors.

The Sérsic profile is very cuspy in the center for
:math:`n>2`. To avoid huge wings in the FFTs when convolving the profile
with the PSF, the profile is split between a 3rd order polynomial,
analytically fit to match, in intensity and its 1st and 2nd spatial
derivatives, the Sérsic profile at :math:`R=4\,\rm pixels`,
:math:`I(r) = I_0 + (r/a)^3`, which has zero first and 2nd derivative at
the center, i.e. a homogeneous core on one hand, and a residual with
finite extent on the other.

For the fit of the spheroid component, the apparent ellipticity allowed
is taken in the range :math:`[0.5, 2]` . This obviously forbids very
flat spheroids to avoid confusion with a flattened disk. By allowing
ellipticities greater than unity, SExtractor avoids dichotomies of
position angle when the ellipticity is very low. The Sérsic index is
allowed values between 1 and 10.

.. _spread_model_def:

Model-based star-galaxy separation: :param:`SPREAD_MODEL`
---------------------------------------------------------

The :param:`SPREAD_MODEL` estimator has been developed as a star/galaxy classifier for the DESDM pipeline :cite:`2012SPIE_8451E_0DM`, and has also been used in other surveys :cite:`2012ApJ_757_83D,2013AA_554A_101B`.
:param:`SPREAD_MODEL` indicates which of the best fitting local PSF model resampled at the current position :math:`\tilde{\boldsymbol{\phi}}` (representing a point source) or a slightly ``fuzzier'' resampled model :math:`\tilde{\boldsymbol{G}}` (representing a galaxy) matches best the image data.
:math:`\tilde{\boldsymbol{G}}` is obtained by convolving the local PSF model with a circular exponential model with scalelength = 1/16 |FWHM|, and resampling the result at the current position on the pixel grid. :param:`SPREAD_MODEL` is normalized to allow comparing sources with different PSFs throughout the field:

.. math::
  :label: spread_model

  {\tt SPREAD\_MODEL} = \frac{\tilde{\boldsymbol{G}}^\mathsf{T} {\bf W}\,\boldsymbol{p}}{\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf W}\,\boldsymbol{p}}
    - \frac{\tilde{\boldsymbol{G}}^\mathsf{T} {\bf W}\,\tilde{\boldsymbol{\phi}}}{\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf W}\,\tilde{\boldsymbol{\phi}}},

where :math:`\boldsymbol{p}` is the image vector centered on the source.
:math:`{\bf W}` is a weight matrix constant along the diagonal except for bad pixels where the weight is 0.

.. warning::
  The definition of :param:`SPREAD_MODEL` above differs from the one given in previous papers, which was incorrect, although in practice both estimators give very similar results.

By construction, :param:`SPREAD_MODEL` is close to zero for point sources, positive for extended sources (galaxies), and negative for detections smaller than the |PSF|, such as cosmic ray hits.
On images taken with linear detectors, the average :param:`SPREAD_MODEL` of point sources should not depend on flux nor |SNR|.
This property may be used to identify bad exposures or |PSF| modeling issues (:numref:`fig_spread`).
More importantly, this makes :param:`SPREAD_MODEL` a very convenient estimator for star-galaxy classification, using a positive threshold to identify extended sources.

.. _fig_spread:

.. figure:: figures/spread.png
   :figwidth: 100%
   :align: center

   |SNR| *vs* :param:`SPREAD_MODEL` for three exposures from :cite:`2013AA_554A_101B`.
   *Left plot*: "good" exposure; extended sources (galaxies and nebulous features) are on the right handside of the stellar locus, and electronic glitches create a small cloud of points on the left handside.
   *Middle*: exposure with an unusual amount of electronic glitches.
   *Right*: exposure with tracking/guiding issues; the |PSF| is too complex for individual sources to be identified as a single objects.


The |pdf| of :param:`SPREAD_MODEL` is close to Gaussian for isolated point sources at a given |SNR|; it gets larger for the fainter sources because of image noise.
In order to maintain a certain level of purity or completeness across the whole magnitude range, it is therefore necessary to take into account the uncertainty on :param:`SPREAD_MODEL`, which can be estimated by propagating the uncertainties on individual pixel values:

.. math::
  :label: spreaderr_model

  \begin{eqnarray}
  {\tt SPREADERR\_MODEL} & = & \frac{1}{(\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf W}\,\boldsymbol{p})^2} \left((\tilde{\boldsymbol{G}}^\mathsf{T} {\bf V}\,\tilde{\boldsymbol{G}})\,(\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf W}\,\boldsymbol{p})^2\right.\nonumber \\
  & & + (\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf V}\,\tilde{\boldsymbol{\phi}})\,(\tilde{\boldsymbol{G}}^\mathsf{T} {\bf W}\,\boldsymbol{p})^2\nonumber \\
  & & \left. - 2\,(\tilde{\boldsymbol{G}}^\mathsf{T} {\bf V}\,\tilde{\boldsymbol{\phi}})\,(\tilde{\boldsymbol{G}}^\mathsf{T} {\bf W}\,\boldsymbol{p})\,(\tilde{\boldsymbol{\phi}}^\mathsf{T} {\bf W}\,\boldsymbol{p}) \right)^{1/2},
  \end{eqnarray}

where :math:`{\bf V}` is the noise covariance matrix, which we assume to be diagonal.
In practice, one may for instance adopt a threshold for star-galaxy separation which is a combination of a fixed and a variable components, such as :math:`\sqrt{\epsilon^2 + (\kappa \times {\tt SPREADERR\_MODEL})^2}`, with :math:`\epsilon \approx 5.10^{-3}` and :math:`\kappa \approx 4`.

..
   Models are measured according to the following table.

   \begin{aligned}
   \hbox{{\tt FLUX\_BACKOFFSET} or {\tt FLUXERR\_BACKOFFSET}} &\to& \hbox{background}
   \nonumber \\
   \hbox{{\tt DISK\_xxx}} &\to& \hbox{exponential disk} \nonumber \\
   \hbox{{\tt SPHEROID\_SERSICN} or {\tt SPHEROID\_SERSICNERR}} &\to&
   \hbox{S\'ersic} \nonumber \\
   \hbox{{\tt SPHEROID\_xxx} without {\tt SPEHEROID\_SERSICN[ERR]}} &\to&
   \hbox{de Vaucouleurs (}n=4 \hbox{ S\'ersic)} \nonumber \\ 
   \hbox{{\tt MODEL\_xxx} only} &\to& \hbox{S\'ersic [???]} \nonumber \\
   \hbox{{\tt SPHEROID\_xxx} and {\tt DISK\_xxx}}&\to& \hbox{S\'ersic spheroid + exponential disk [???]}   \nonumber \end{aligned}

  Table [modeltriggers] should be interpreted as meaning that if one of
  the parameters given in the parameter file (e.g. default.param) includes
  the string on the left of the arrow, the model to the right of the arrow
  is triggered. For example, when including parameters that contain the
  string ‘MODEL’, both galaxies and stars are fit with convolutions of
  Sérsic models with the PSF. If no SPHEROID\_xxx or DISK\_xxx parameter
  is present, but the model-fitting process is nevertheless triggered by
  the presence of other measurement parameters or relevant
  CHECKIMAGE\_TYPEs , a single component with Sérsic profile and
  adjustable Sérsic index :math:`n` is fitted.

  The number of parameters that are fit are 2 for the global center, 4 per
  model for the scale, normalization, aspect ratio and position angle,
  plus the index for the Sérsic model. For example, fitting a Sérsic +
  exponential disk involves a fitting 11 parameters.

  The measurement parameters related to model-fitting follow the usual
  SExtractor rules:

  Flux measurements are available in ADUs (FLUX\_xxx parameters) or
  magnitudes (MAG\_xxx parameters), Coordinates and radii are available in
  pixels or celestial units (provided that the FITS image header contains
  the appropriate WCS information).

  xxxMODEL\_yyy measurement parameters deal with the global fitted model,
  i.e. the sum of all components (e.g. chi-square per d.o.f. CHI2\_MODEL,
  PSF-corrected ellipticities E1/2MODEL\_IMAGE, EPS1/2 MODEL\_IMAGE).

  :math:`1\,\sigma` error estimates xxxERR\_yyy are provided for most
  measurement parameters; they are obtained by marginalizing the full
  covariance matrix of the fit.

  Since the model fitting involves convolution with the PSF, it is
  imperative to launch PSFEx before launching SExtractor. In practice, the
  sequence of operations is:

  #. Run SExtractor to prepare PSFEx;

  #. Run PSFEx to prepare model fits in SExtractor;

  #. Run SExtractor with model fit parameters.
