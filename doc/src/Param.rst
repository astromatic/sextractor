.. File Param.rst

The catalog parameter file
============================

In addition to the configuration file detailed above, |SExtractor| requires a
file containing the list of parameters that will be listed in the output
catalog for every detection. This allows the software to compute only
catalog parameters that are needed. The name of this
catalog-parameter file is traditionally suffixed with ``.param``, and must
be specified using the ``PARAMETERS_NAME`` config parameter. The full set
of parameters can be queried with the command

.. code-block:: console

 $ sex -dp

Format
------

The format of the catalog parameter list is ASCII, and there must be
*a single keyword per line*. Presently two kinds of keywords are
recognized by |SExtractor|: scalars and vectors. Scalars, like ``X_IMAGE``,
produce single numbers in the output catalog. Vectors, like ``MAG_APER(4)``
or ``VIGNET(15,15)``, produce arrays of numbers. The ordering of measurements
in the output catalog is identical to that of the keywords in the parameter
list. Comments are allowed, they must begin with a ``#``.

Variants
--------

For many catalog parameters, especially those related to flux,
position, or shape, several variants of the same measurement are
available:

Fluxes
~~~~~~

Fluxes may be expressed in linear (ADU) units or as Pogson
:cite:`1856MNRAS_17_12P` magnitudes. Flux measurements in ADUs
are prefixed with ``FLUX_``, for example: ``FLUX_AUTO``, ``FLUX_ISO``, etc.
Magnitudes are prefixed with ``MAG_`` e.g., ``MAG_AUTO``, ``MAG_ISO``, ... In
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

Flux uncertainties follow a scheme similar to that of fluxes. Flux uncertainties are
prefixed with ``FLUXERR_``, as in ``FLUXERR_AUTO`` or ``FLUXERR_ISO``. Magnitude
uncertainties start with ``MAGERR_``, for instance: ``MAGERR_AUTO``,
``MAGERR_ISO``,... Magnitude uncertainties :math:`\sigma_m` are derived
from the estimated 1-\ :math:`\sigma` flux error :math:`\sigma_f`:

.. math::

   \sigma_m = \left\{\begin{array}{ll}
   (2.5/\ln 10) (\sigma_f/f)\ &\mbox{if } f > 0\\
   99.0 &\mbox{otherwise}.
   \end{array}\right.

.. _coord_suffix:

Positions and shapes
~~~~~~~~~~~~~~~~~~~~

Positions, distances and position angles are computed in pixel coordinates. They may be expressed in image pixels, world coordinates, or in celestial coordinates, depending on the suffix:

_IMAGE
  Measurements are given in pixel coordinates, in units of pixels. For example: ``Y_IMAGE``, ``ERRAWIN_IMAGE``, ``THETA_IMAGE`` etc. Following the FITS convention, in |SExtractor| the center of the first image pixel has coordinates (1.0,1.0). Position angles are counted from the *x* axis (axis 1), positive towards the *y* axis (axis 2)

_WORLD
  Measurements are given in so-called “world coordinates”, converted from pixel coordinates using the local Jacobian of the transformation between both systems. This requires World Coordinate System (|WCS|_) metadata :cite:`2002AA_395_1061G` to be present in the FITS image header(s). Position angles are counted from the first world axis, positive towards the second world axis.

_SKY, _J2000, _B1950
  Measurements are given in celestial (equatorial) coordinates, converted from pixel coordinates using the local Jacobian of the transformation between both systems. Positions and distances are in units of degrees. This requires celestial |WCS| metadata :cite:`2002AA_395_1077C` to be present in the FITS image header(s). _SKY measurements are given in the native world coordinate system. _J2000 and _B1950 measurements are automatically converted from the native |WCS|, taking into account the change of reference frame. In all cases, positions angles are counted East-of-North.

.. include:: keys.rst

