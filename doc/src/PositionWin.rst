.. File PositionWin.rst

Windowed positional parameters
==============================

Measurements performed through a *window* function (an *envelope*) do not have many of the drawbacks of isophotal measurements. |SExtractor| implements “windowed” versions for most
of the measurements described in the :ref:`previous section<position_iso>`:

.. note::
  Unless otherwise noted, all parameter names given below are only prefixes. They must be followed by _IMAGE if the results shall be expressed in pixel coordinates or _WORLD, _SKY, _J2000 or _B1950 for |WCS|_ coordinates (see :ref:`coord_suffix`).

+----------------------+-------------------------------+
| Isophotal parameters | Equivalent windowed parameters|
+======================+===============================+
| X, Y                 | XWIN, YWIN                    |
+----------------------+-------------------------------+
| ERRA, ERRB, ERRTHETA | ERRAWIN, ERRBWIN, ERRTHETAWIN |
+----------------------+-------------------------------+
| A, B, THETA          | AWIN, BWIN, THETAWIN          |
+----------------------+-------------------------------+
| X2, Y2, XY           | X2WIN, Y2WIN, XYWIN           |
+----------------------+-------------------------------+
| CXX, CYY, CXY        | CXXWIN, CYYWIN, CXYWIN        |
+----------------------+-------------------------------+

The computations involved are roughly the same except that the pixel
values are integrated within a circular Gaussian window as opposed to
the object’s isophotal footprint. The Gaussian window is scaled to each
object; its FWHM is the diameter of the disk that contains half of the
object flux (:math:`d_{50}`). Note that in double-image mode
(§[chap:using]) the window is scaled based on the *measurement* image.

.. _xywin:

Windowed centroid: XWIN, YWIN
-----------------------------

This is an iterative process. The computation starts by initializing the
windowed centroid coordinates :math:`\overline{x_{\tt WIN}}^{(0)}` and
:math:`\overline{y_{\tt WIN}}^{(0)}` to their basic :math:`\overline{x}`
and :math:`\overline{y}` isophotal equivalents, respectively. Then at
each iteration :math:`t`, :math:`\overline{x_{\tt WIN}}` and
:math:`\overline{y_{\tt WIN}}` are refined using:

.. math::
  :label: xywin

   \begin{aligned}
   {\tt XWIN}^{(t+1)} & = & \overline{x_{\tt WIN}}^{(t+1)}
   = \overline{x_{\tt WIN}}^{(t)} + 2\,\frac{\sum_{r_i^{(t)} < r_{\rm max}}
   w_i^{(t)} I_i \ (x_i  - \overline{x_{\tt WIN}}^{(t)})}
   {\sum_{r_i^{(t)} < r_{\rm max}} w_i^{(t)} I_i},\\
   \label{eq:ywin}
   {\tt YWIN}^{(t+1)} & = & \overline{y_{\tt WIN}}^{(t+1)}
   = \overline{y_{\tt WIN}}^{(t)} + 2\,\frac{\sum_{r_i^{(t)} < r_{\rm max}}
   w_i^{(t)} I_i\ (y_i - \overline{y_{\tt WIN}}^{(t)})}
   {\sum_{r_i^{(t)} < r_{\rm max}} w_i^{(t)} I_i},\end{aligned}

where

.. math::
  :label: wi_t

  w_i^{(t)} = \exp \left(-\frac{r_i^{(t)^2}}{2s_{\tt WIN}^2} \right),

with

.. math::
  :label: ri_t

   r_i^{(t)} = \sqrt{\left(x_i - \overline{x_{\tt WIN}}^{(t)}\right)^2 + \left(y_i
   - \overline{y_{\tt WIN}}^{(t)}\right)^2}

and :math:`s_{\tt WIN} = d_{50} / \sqrt{8 \ln 2}`. The process stops
when the change in position between two iterations is less than
:math:`2\times10^{-4}` pixel, a condition which is generally achieved in
about 3 to 5 iterations.

Although the iterative nature of the processing slows down the
processing , it is recommended to use whenever possible windowed
parameters instead of their isophotal equivalents, since the
measurements they provide are much more precise (:numref:`fig_xwinprec`).
The precision in centroiding offered by XWIN_IMAGE and YWIN_IMAGE is
actually very close to that of PSF-fitting on focused and properly
sampled star images, and can also be applied to galaxies. It has been
verified that for isolated, Gaussian-like PSFs, its accuracy is close to
the theoretical limit set by image noise [#win_accuracy]_.

.. _fig_xwinprec:

.. figure:: figures/xwinprec.*
   :figwidth: 100%
   :align: center

   Comparison between isophotal and windowed centroid measurement residuals on
   simulated, background noise-limited images.
   *Left*: histogram of the difference between X_IMAGE and the true centroid
   in x.
   *Right*: histogram of the difference between XWIN_IMAGE and the true centroid
   in x.

Windowed 2nd order moments: X2, Y2, XY
--------------------------------------

Windowed second-order moments are computed on the image data once the :ref:`centering process <xywin>` has converged:

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

Windowed second-order moments are typically twice smaller than their
isophotal equivalent.

Windowed shape parameters: AWIN, BWIN, THETAWIN
-----------------------------------------------

They are computed from the windowed 2nd order moments exactly the same
way as in :eq:`theta0_3` and :eq:`aimage_2` from the :ref:`isophotal shape parameter<shape_iso>` section:

.. math::
  :label: shapewin

   \begin{aligned}
   {\tt AWIN}^2 & = & \frac{\overline{x_{\tt WIN}^2}+\overline{y_{\tt WIN}^2}}{2}
       + \sqrt{\left(\frac{\overline{x_{\tt WIN}^2}-\overline{y_{\tt WIN}^2}}{2}\right)^2 + \overline{xy_{\tt WIN}}^2},\\
   {\tt BWIN}^2 & = & \frac{\overline{x_{\tt WIN}^2}+\overline{y_{\tt WIN}^2}}{2}
       - \sqrt{\left(\frac{\overline{x_{\tt WIN}^2}-\overline{y_{\tt WIN}^2}}{2}\right)^2 + \overline{xy_{\tt WIN}}^2},\\
   \tan (2\,{\tt THETAWIN}) & = & 2 \frac{\overline{xy_{\tt WIN}}}{\overline{x_{\tt WIN}^2} - \overline{y_{\tt WIN}^2}}.
   \end{aligned}


Windowed ellipse parameters: CXXWIN, CYYWIN, CXYWIN
---------------------------------------------------

They are computed from the windowed 2nd order moments exactly the same
way as in :eq:`ellipse_2` describing the :ref:`isophotal ellipse parameters<ellipse_iso>`:

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
   \frac{\overline{xy_{\tt WIN}}}{\overline{x_{\tt WIN}^2} \overline{y_{\tt WIN}^2} - \overline{xy_{\tt WIN}}^2}\end{aligned}

Windowed position uncertainties: ERRX2WIN, ERRY2WIN, ERRXYWIN, ERRAWIN, ERRBWIN, ERRTHETAWIN, ERRCXXWIN, ERRCYYWIN, ERRCXYWIN
-----------------------------------------------------------------------------------------------------------------------------

Windowed position uncertainties are computed on the image data once the
centering process of the :ref:`windowed centroid <xywin>` has converged. Assuming that
noise is uncorrelated among pixels, standard error propagation applied
to :eq:`xywin` writes:

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
   {\left(\sum_{r_i < r_{\rm max}} w_i I_i\right)^2}.\end{aligned}

The semi-major axis ERRAWIN, semi-minor axis ERRBWIN, and position angle
ERRTHETAWIN of the :math:`1\sigma` position error ellipse are computed
from the covariance matrix elements
:math:`{\rm var}(\overline{x_{\tt WIN}})`,
:math:`{\rm var}(\overline{y_{\tt WIN}})`,
:math:`{\rm cov}(\overline{x_{\tt WIN}},\overline{y_{\tt WIN}})`,
similarly to the :ref:`isophotal error ellipse <poserr>`:

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

And the error ellipse parameters are:

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



.. [#win_accuracy] See http://www.astromatic.net/forum/showthread.php?tid=581 .

.. include:: keys.rst

