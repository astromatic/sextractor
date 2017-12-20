.. File Model.rst

.. include:: global.rst

Model fitting
=============

Fitting procedure
-----------------

SExtractor can fit models to the images of detected objects since version 2.8. The fit is performed by minimizing the loss function

.. math::
  :label: loss_func

  \lambda(\boldsymbol{q}) = \sum_i \left(g\left(\frac{p_i - \tilde{m}_i(\boldsymbol{q})}{\sigma_i}\right)\right)^2 + \sum_j \frac{q_j - \mu_j}{}

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

Regularization
~~~~~~~~~~~~~~

.. _model_minimization_def:

Minimization
~~~~~~~~~~~~

Minimization of the loss function :math:`\lambda(\boldsymbol{q})` is carried out using the `Levenberg-Marquardt algorithm <http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm>`_, and more specifically the |LevMar|_ implementation :cite:`lourakis04LM`.
The fit is done inside a disk which diameter is scaled to include the isophotal footprint of the object, plus the FWHM of the PSF, plus a 20 % margin.
The number of iterations is returned in the :param:`NITER_MODEL` measurement parameter.
It is generally a few tens.
The final value of the modified chi square term in :eq:`loss_func`, divided by the number of degrees of freedom, is returned in :param:`CHI2_MODEL`.
The :param:`FLAGS_MODEL` parameter flags various issues which may happen during the fitting process (see the flags section for details on how flags are managed in |SExtractor|):

.. csv-table:: :param:`FLAGS_MODEL` flag description
  :header: "Value", "Meaning"
  :widths: 3 60

  1, "the unconvolved, supersampled model raster exceeds 512×512 pixels and had to be resized"
  2, "the convolved, resampled model raster exceeds 512×512 pixels and had to be resized"
  4, "not enough pixels are available for model fitting on the measurement image (less pixels than fit parameters)"
  8, "at least one of the fitted parameters hits the lower bound"
  16, "at least one of the fitted parameters hits the upper bound"

:math:`1\,\sigma` error estimates are provided for most measurement parameters; they are obtained by marginalizing the full covariance matrix of the fit.

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
