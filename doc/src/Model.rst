.. File Model.rst

.. include:: global.rst

Model fitting
=============

SExtractor can fit models to the images of detected objects since version 2.8. The fit is performed by minimizing the loss function

.. math::
  :label: loss_func

  \lambda(\boldsymbol{q}) = \sum_i \left(g\left(\frac{p_i - \hat{m}_i(\boldsymbol{q})}{\sigma_i}\right)\right)^2 + \sum_j \frac{q_j - \mu_j}{}

with respect to components of the model parameter vector :math:`\boldsymbol{q}`. :math:`\boldsymbol{q}` comprises parameters describing the shape of the model and the model pixel coordinates :math:`\boldsymbol{x}`.

The first term in :eq:`loss_func` is a modified `weighted sum of squares <http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>`_ that aims at minimizing the residuals of the fit. :math:`p_i`, :math:`\hat{m}_i(\boldsymbol{q})` and :math:`\sigma_i` are respectively the pixel value above the background, the value of the resampled model, and the pixel value uncertainty at image pixel :math:`i`.
:math:`g(u)` is a derivable monotonous function that reduces the influence of large deviations from the model (e.g., contamination by neighbors):

.. math::
  :label: loss_func

  g(u) = \left\{
    \begin{array}{rl}
       u_0 \log \left(1 + \frac{u}{u_0}\right) & \mbox{if } u \ge 0,\\
      -u_0 \log \left(1 - \frac{u}{u_0}\right) & \mbox{otherwise.}\\
    \end{array}
  \right.

The vector :math:`\hat{\boldsymbol{m}}(\boldsymbol{q})` is obtained by convolving the high resolution model :math:`\boldsymbol{m}(\boldsymbol{q})` with the local PSF model :math:`\boldsymbol{\phi}` and applying a resampling operator :math:`\mathbf{R}(\boldsymbol{x})` to generate the final model raster at position :math:`\boldsymbol{x}` at the nominal image resolution:

.. math::
  :label: model_resampling

  \hat{\boldsymbol{m}}(\boldsymbol{q}) = \mathbf{R}(\boldsymbol{x}) (\boldsymbol{m}(\boldsymbol{q})*\boldsymbol{\phi}).

Levenberg-Marquardt minimization, inside a disk which diameter is scaled to include the isophotal footprint plus a 20 % margin, plus the size of the PSF model image.

The models that can be fit are:

- Exponential disk

  .. math::

      \Sigma_{\tt ExpDisk}(R) = \Sigma(0) \exp \left (- {R\over h}\right ) 
      \label{expdisk}

- Sérsic (:math:`R^{1/n}`) spheroid (bulg)

  .. math::

      \Sigma_{\tt Sersic}(R) = \Sigma(0) \exp \left [- b(n)\,\left({R\over
          R_e}\right)^{1/n}\right ] \ ,
      \label{sersic}

where, for the :raw-latex:`\cite{sersic:1968}` model, :math:`b(n)` is the solution of

  .. math::

      2 \gamma[2\,n,b(n)] = \Gamma(2\,n)
      \label{bofn}

  An accurate approximation for the solution for :math:`b(n)` of equation (bofn) is :raw-latex:`\citep{ciotti:bertin:1999}`

   .. math::

      b(n) = 2\,n - {1\over3} + {4\over 405\,n} + {46\over 25515\,n^2} + {131\over
        1148175\,n^3}

-  :raw-latex:`\cite{devaucouleurs48}` spheroid (bulge, eq. [[sersic]], with :math:`n=4`)

-  Exponential disk + Sérsic (:math:`R^{1/n}`) spheroid (bulge)

-  Point source

-  Background (constant)

For these models, SExtractor can compute fluxes and magnitudes, as well
as sizes (disk scale length for the disks and effective — projected
half-light — radii for the spheroids), characteristic surface
magnitudes, and Sérsic index, as well as their uncertainties.

The models are concentric (they assume the same center) and are all
convolved with the PSF, given by the .psf file, which must be determined
by first running PSFEx (see below).

Unfortunately, the Sérsic profile is very cuspy in the center for
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

Models are measured according to the following table.

.. math::

   \begin{aligned}
   \hbox{{\tt FLUX\_BACKOFFSET} or {\tt FLUXERR\_BACKOFFSET}} &\to& \hbox{background}
   \nonumber \\
   \hbox{{\tt DISK\_xxx}} &\to& \hbox{exponential disk} \nonumber \\
   \hbox{{\tt SPHEROID\_SERSICN} or {\tt SPHEROID\_SERSICNERR}} &\to&
   \hbox{S\'ersic} \nonumber \\
   \hbox{{\tt SPHEROID\_xxx} without {\tt SPEHEROID\_SERSICN[ERR]}} &\to&
   \hbox{de Vaucouleurs (}n=4 \hbox{ S\'ersic)} \nonumber \\ 
   \hbox{{\tt MODEL\_xxx} only} &\to& \hbox{S\'ersic [???]} \nonumber \\
   \hbox{{\tt SPHEROID\_xxx} and {\tt DISK\_xxx}}&\to& \hbox{S\'ersic spheroid + exponential disk [???]} \nonumber \end{aligned}

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

Experience shows that the de Vaucouleurs spheroid + exponential disk
combination provides fairly accurate and robust fits for moderately
resolved faint galaxies. An adjustable Sérsic index may offer lower
residuals on spheroids and/or well-resolved galaxies, but makes the fit
less robust and more sensitive to PSF model errors. One might think of
adding some mechanism to lock or unlock the Sérsic index automatically
in future versions of SExtractor.

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
