.. File PositionWin.rst

.. include:: global.rst

Windowed positional parameters
==============================

Measurements performed through a *window* function (an *envelope*) do not have many of the drawbacks of isophotal measurements.
|SExtractor| implements “windowed” versions for most of the measurements described in the :ref:`previous section<pos_iso_def>`:

.. note::
  Unless otherwise noted, all parameter names given below are only prefixes. They must be followed by _IMAGE if the results shall be expressed in pixel coordinates or :param:`_WORLD`, :param:`_SKY`, :param:`_J2000` or :param:`_B1950` for |WCS|_ coordinates (see :ref:`coord_suffix`).

+-------------------------------------------------+----------------------------------------------------------+
| Isophotal parameters                            | Equivalent windowed parameters                           |
+=================================================+==========================================================+
| :param:`X`, :param:`Y`                          | :param:`XWIN`, :param:`YWIN`                             |
+-------------------------------------------------+----------------------------------------------------------+
| :param:`ERRA`, :param:`ERRB`, :param:`ERRTHETA` | :param:`ERRAWIN`, :param:`ERRBWIN`, :param:`ERRTHETAWIN` |
+-------------------------------------------------+----------------------------------------------------------+
| :param:`A`, :param:`B`, :param:`THETA`          | :param:`AWIN`, :param:`BWIN`, :param:`THETAWIN`          |
+-------------------------------------------------+----------------------------------------------------------+
| :param:`X2`, :param:`Y2`, :param:`XY`           | :param:`X2WIN`, :param:`Y2WIN`, :param:`XYWIN`           |
+-------------------------------------------------+----------------------------------------------------------+
| :param:`CXX`, :param:`CYY`, :param:`CXY`        | :param:`CXXWIN`, :param:`CYYWIN`, :param:`CXYWIN`        |
+-------------------------------------------------+----------------------------------------------------------+

The computations involved are roughly the same except that the pixel values are integrated within a circular Gaussian window as opposed to the object’s isophotal footprint.
The Gaussian window is scaled to each object; the Gaussian FWHM is the diameter of the disk that contains half of the object flux (:math:`d_{50}`).
Note that in double-image mode (§[chap:using]) the window is scaled based on the *measurement* image.

.. _pos_win_def:

Windowed centroid: :param:`XWIN`, :param:`YWIN`
-----------------------------------------------

This is an iterative process.
Computation starts by initializing the windowed centroid coordinates :math:`\overline{x_{\tt WIN}}^{(0)}` and :math:`\overline{y_{\tt WIN}}^{(0)}` to their basic :math:`\overline{x}` and :math:`\overline{y}` isophotal equivalents, respectively.
Then at each iteration :math:`t`, :math:`\overline{x_{\tt WIN}}` and :math:`\overline{y_{\tt WIN}}` are refined using:

.. math::
  :label: xywin

   \begin{aligned}
   {\tt XWIN}^{(t+1)} & = & \overline{x_{\tt WIN}}^{(t+1)}
   = \overline{x_{\tt WIN}}^{(t)} + 2\,\frac{\sum_{r_i^{(t)} < r_{\rm max}}
   w_i^{(t)} I_i \ (x_i  - \overline{x_{\tt WIN}}^{(t)})}
   {\sum_{r_i^{(t)} < r_{\rm max}} w_i^{(t)} I_i},\\
   {\tt YWIN}^{(t+1)} & = & \overline{y_{\tt WIN}}^{(t+1)}
   = \overline{y_{\tt WIN}}^{(t)} + 2\,\frac{\sum_{r_i^{(t)} < r_{\rm max}}
   w_i^{(t)} I_i\ (y_i - \overline{y_{\tt WIN}}^{(t)})}
   {\sum_{r_i^{(t)} < r_{\rm max}} w_i^{(t)} I_i},
   \end{aligned}

where

.. math::
  :label: wi_t

  w_i^{(t)} = \exp \left(-\frac{r_i^{(t)^2}}{2s_{\tt WIN}^2} \right),

with

.. math::
  :label: ri_t

   r_i^{(t)} = \sqrt{\left(x_i - \overline{x_{\tt WIN}}^{(t)}\right)^2 + \left(y_i
   - \overline{y_{\tt WIN}}^{(t)}\right)^2}

and :math:`s_{\tt WIN} = d_{50} / \sqrt{8 \ln 2}`.
Process stops when the change in position between two iterations is less than :math:`2\times10^{-4}` pixel, a condition which is achieved in about 3 to 5 iterations in practice.

Although they are slower, it is recommended to use whenever possible windowed position parameters instead of their isophotal equivalents; the measurements they provide are generally much more accurate (:numref:`fig_xwinprec`).
The centroiding accuracy of :param:`XWIN_IMAGE` and :param:`YWIN_IMAGE` is actually very close to that of PSF-fitting on focused and properly sampled star images. Windowed measurements can also be applied to galaxies.
It has been verified that for isolated objects with Gaussian-like profiles, their accuracy is close to the theoretical limit set by image noise [#win_accuracy]_.

.. _fig_xwinprec:

.. figure:: figures/xwinprec.*
   :figwidth: 100%
   :align: center

   Comparison between isophotal and windowed centroid measurement residuals on simulated, background noise-limited images.
   *Left*: histogram of the difference between :param:`X_IMAGE` and the true centroid in x.
   *Right*: histogram of the difference between :param:`XWIN_IMAGE` and the true centroid in x.

.. _moments_win_def:

Windowed 2nd order moments: :param:`X2`, :param:`Y2`, :param:`XY`
-----------------------------------------------------------------

Windowed second-order moments are computed on the image data once the :ref:`centering process <pos_win_def>` has converged:

.. math::
  :label: x2y2win

   \begin{aligned}
   {\tt X2WIN} & = \overline{x_{\tt WIN}^2}
   = & \frac{\sum_{r_i < r_{\rm max}} w_i I_i (x_i - \overline{x_{\tt WIN}})^2}
   {\sum_{r_i < r_{\rm max}} w_i I_i},\\
   {\tt Y2WIN} & = \overline{y_{\tt WIN}^2}
   = & \frac{\sum_{r_i < r_{\rm max}} w_i I_i (y_i - \overline{y_{\tt WIN}})^2}
   {\sum_{r_i < r_{\rm max}} w_i I_i},\\
   {\tt XYWIN} & = \overline{xy_{\tt WIN}}
   = & \frac{\sum_{r_i < r_{\rm max}} w_i I_i (x_i - \overline{x_{\tt WIN}})
   (y_i - \overline{y_{\tt WIN}})}
   {\sum_{r_i < r_{\rm max}} w_i I_i}.\end{aligned}

Windowed second-order moments are typically twice smaller than their isophotal equivalent.

.. _shape_win_def:

Windowed shape parameters: :param:`AWIN`, :param:`BWIN`, :param:`THETAWIN`
--------------------------------------------------------------------------

They are computed from the windowed 2nd order moments exactly the same
way as in :eq:`theta0_3` and :eq:`aimage_2` from the :ref:`isophotal shape parameter<shape_iso_def>` section:

.. math::
  :label: shapewin

   \begin{aligned}
   {\tt AWIN}^2 & = & \frac{\overline{x_{\tt WIN}^2}+\overline{y_{\tt WIN}^2}}{2}
       + \sqrt{\left(\frac{\overline{x_{\tt WIN}^2}-\overline{y_{\tt WIN}^2}}{2}\right)^2 + \overline{xy_{\tt WIN}}^2},\\
   {\tt BWIN}^2 & = & \frac{\overline{x_{\tt WIN}^2}+\overline{y_{\tt WIN}^2}}{2}
       - \sqrt{\left(\frac{\overline{x_{\tt WIN}^2}-\overline{y_{\tt WIN}^2}}{2}\right)^2 + \overline{xy_{\tt WIN}}^2},\\
   \tan (2\,{\tt THETAWIN}) & = & 2 \frac{\overline{xy_{\tt WIN}}}{\overline{x_{\tt WIN}^2} - \overline{y_{\tt WIN}^2}}.
   \end{aligned}


.. _ellipse_win_def:

Windowed ellipse parameters: :param:`CXXWIN`, :param:`CYYWIN`, :param:`CXYWIN`
------------------------------------------------------------------------------

Ellipse parameters are computed from the windowed 2nd order moments exactly the same way as in :eq:`ellipse_2` describing the :ref:`isophotal ellipse parameters<ellipse_iso_def>`:

.. math::
  :label: ellipsewin_2

   \begin{aligned}
   {\tt CXXWIN} & = & \frac{\cos^2 {\tt THETAWIN}}{{\tt AWIN}^2} + \frac{\sin^2
   {\tt THETAWIN}}{{\tt BWIN}^2} =
   \frac{\overline{y_{\tt WIN}^2}}{\overline{x_{\tt WIN}^2} \overline{y_{\tt WIN}^2} - \overline{xy_{\tt WIN}}^2}\\
   {\tt CYYWIN} & = & \frac{\sin^2 {\tt THETAWIN}}{{\tt
   AWIN}^2} + \frac{\cos^2 {\tt THETAWIN}}{{\tt BWIN}^2} =
   \frac{\overline{x_{\tt WIN}^2}}{\overline{x_{\tt WIN}^2} \overline{y_{\tt WIN}^2} - \overline{xy_{\tt WIN}}^2}\\
   {\tt CXYWIN} & = & 2 \,\cos {\tt THETAWIN}\,\sin {\tt
   THETAWIN} \left( \frac{1}{{\tt AWIN}^2} - \frac{1}{{\tt BWIN}^2}\right) = -2\,
   \frac{\overline{xy_{\tt WIN}}}{\overline{x_{\tt WIN}^2} \overline{y_{\tt WIN}^2} - \overline{xy_{\tt WIN}}^2}.
   \end{aligned}

.. _poserr_win_def:

Windowed position uncertainties: :param:`ERRX2WIN`, :param:`ERRY2WIN`, :param:`ERRXYWIN`, :param:`ERRAWIN`, :param:`ERRBWIN`, :param:`ERRTHETAWIN`, :param:`ERRCXXWIN`, :param:`ERRCYYWIN`, :param:`ERRCXYWIN`
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Windowed position uncertainties are computed on the image data once the centering process of the :ref:`windowed centroid <pos_win_def>` has converged.
Assuming that noise is uncorrelated among pixels, standard error propagation applied to :eq:`xywin` writes:

.. math::
  :label: errwin

   \begin{aligned}
   {\tt ERRX2WIN} & = {\rm var}(\overline{x_{\tt WIN}})
   = & 4\,\frac{\sum_{r_i < r_{\rm max}} w_i^2 \sigma^2_i (x_i-\overline{x_{\tt WIN}})^2}
   {\left(\sum_{r_i < r_{\rm max}} w_i I_i\right)^2},\\
   {\tt ERRY2WIN} & = {\rm var}(\overline{y_{\tt WIN}})
   = & 4\,\frac{\sum_{r_i < r_{\rm max}} w_i^2 \sigma^2_i (y_i-\overline{y_{\tt WIN}})^2}
   {\left(\sum_{r_i < r_{\rm max}} w_i I_i\right)^2},\\
   {\tt ERRXYWIN} & = {\rm cov}(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})
   = & 4\,\frac{\sum_{r_i < r_{\rm max}}
   w_i^2 \sigma^2_i (x_i-\overline{x_{\tt WIN}})(y_i-\overline{y_{\tt WIN}})}
   {\left(\sum_{r_i < r_{\rm max}} w_i I_i\right)^2}.
   \end{aligned}

Semi-major axis :param:`ERRAWIN`, semi-minor axis :param:`ERRBWIN`, and position angle :param:`ERRTHETAWIN` of the :math:`1\sigma` position error ellipse are computed from the covariance matrix elements :math:`{\rm var}(\overline{x_{\tt WIN}})`, :math:`{\rm var}(\overline{y_{\tt WIN}})`, :math:`{\rm cov}(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})`, similarly to the :ref:`isophotal error ellipse <poserr_iso_def>`:

.. math::
  :label: errabthetawin

   \begin{aligned}
   {\tt ERRAWIN}^2 & = & \frac{{\rm var}(\overline{x_{\tt WIN}})+{\rm var}(\overline{y_{\tt WIN}})}{2}
       + \sqrt{\left(\frac{{\rm var}(\overline{x_{\tt WIN}})-{\rm var}(\overline{y_{\tt WIN}})}{2}\right)^2
       + {\rm cov}^2(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})},\\
   {\tt ERRBWIN}^2 & = & \frac{{\rm var}(\overline{x_{\tt WIN}})+{\rm var}(\overline{y_{\tt WIN}})}{2}
       - \sqrt{\left(\frac{{\rm var}(\overline{x_{\tt WIN}})-{\rm var}(\overline{y_{\tt WIN}})}{2}\right)^2
       + {\rm cov}^2(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})},\\
   \tan (2{\tt ERRTHETAWIN}) & = & 2 \,\frac{{\rm cov}(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})}
                       {{\rm var}(\overline{x_{\tt WIN}}) - {\rm var}(\overline{y_{\tt WIN}})}.
    \end{aligned}

The error ellipse parameters are:

.. math::
  :label: errellipsewin

   \begin{aligned}
   {\tt ERRCXXWIN} & = & \frac{\cos^2 {\tt ERRTHETAWIN}}{{\tt ERRAWIN}^2} +
   \frac{\sin^2 {\tt ERRTHETAWIN}}{{\tt ERRBWIN}^2} = \frac{{\rm
   var}(\overline{y_{\tt WIN}})}{{\rm var}(\overline{x_{\tt WIN}}) {\rm var}(\overline{y_{\tt WIN}}) -
   {\rm cov}^2(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})},\\
   {\tt ERRCYYWIN} & = & \frac{\sin^2 {\tt ERRTHETAWIN}}{{\tt ERRAWIN}^2} +
   \frac{\cos^2 {\tt ERRTHETAWIN}}{{\tt ERRBWIN}^2} =
   \frac{{\rm var}(\overline{x_{\tt WIN}})}{{\rm var}(\overline{x_{\tt WIN}}) {\rm var}(\overline{y_{\tt WIN}}) -
   {\rm cov}^2(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})},\\
   {\tt ERRCXYWIN} & = & 2 \cos {\tt
   ERRTHETAWIN}\sin {\tt ERRTHETAWIN} \left( \frac{1}{{\tt ERRAWIN}^2} -
   \frac{1}{{\tt ERRBWIN}^2}\right)\\ & = & -2 \frac{{\rm
   cov}(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})}{{\rm var}(\overline{x_{\tt WIN}}) {\rm var}(\overline{y_{\tt WIN}}) -
   {\rm cov}^2(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})}.
   \end{aligned}

.. _flags_win_def:

Windowed measurement flags: :param:`FLAGS_WIN`
----------------------------------------------

The :param:`FLAGS_WIN` catalog parameter flags various issues which may happen with windowed measurements (see the :ref:`flagging` section for details on how flags are managed in |SExtractor|):

.. _flags_win_table:

.. csv-table:: :param:`FLAGS_WIN` description
  :header: "Value", "Meaning"
  :widths: 3 60

  1, ":ref:`Windowed second-order moments <moments_win_def>` are inconsistent (:math:`\overline{x_{\tt WIN}^2}\overline{y_{\tt WIN}^2} \le \overline{x_{\tt WIN}y_{\tt WIN}}^2`)"
  2, ":ref:`Windowed second-order moments <moments_win_def>` are less than or equal to 0"
  4, "Windowed flux is less than or equal to 0"

.. [#win_accuracy] See http://www.astromatic.net/forum/showthread.php?tid=581 .

