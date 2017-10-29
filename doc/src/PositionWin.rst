.. File PositionWin.rst

Windowed positional parameters
==============================

Parameters measured within an object’s isophotal limit are sensitive to
two main factors: 1) changes in the detection threshold, which create a
variable bias and 2) irregularities in the object’s isophotal
boundaries, which act as additional “noise” in the measurements.

Measurements performed through a *window* function (an *envelope*) do
not have such drawbacks. |SExtractor| implements “windowed” versions for most
of the measurements described in [chap:isoparam]:

+----------------------------------------+-------------------------------------------------+
| Isophotal parameters                   | Equivalent windowed parameters                  |
+========================================+=================================================+
| X_IMAGE, Y_IMAGE                       | XWIN_IMAGE, YWIN_IMAGE                          |
+----------------------------------------+-------------------------------------------------+
| ERRA_IMAGE, ERRB_IMAGE, ERRTHETA_IMAGE | ERRAWIN_IMAGE, ERRBWIN_IMAGE, ERRTHETAWIN_IMAGE |
+----------------------------------------+-------------------------------------------------+
| A_IMAGE, B_IMAGE, THETA_IMAGE          | AWIN_IMAGE, BWIN_IMAGE, THETAWIN_IMAGE          |
+----------------------------------------+-------------------------------------------------+
| X2_IMAGE, Y2_IMAGE, XY_IMAGE           | X2WIN_IMAGE, Y2WIN_IMAGE, XYWIN_IMAGE           |
+----------------------------------------+-------------------------------------------------+
| CXX_IMAGE, CYY_IMAGE, CXY_IMAGE        | CXXWIN_IMAGE, CYYWIN_IMAGE, CXYWIN_IMAGE        |
+----------------------------------------+-------------------------------------------------+

The computations involved are roughly the same except that the pixel
values are integrated within a circular Gaussian window as opposed to
the object’s isophotal footprint. The Gaussian window is scaled to each
object; its FWHM is the diameter of the disk that contains half of the
object flux (:math:`d_{50}`). Note that in double-image mode
(§[chap:using]) the window is scaled based on the *measurement* image.

Windowed centroid: XWIN, YWIN
-----------------------------

This is an iterative process. The computation starts by initializing the
windowed centroid coordinates :math:`\overline{x_{\tt WIN}}^{(0)}` and
:math:`\overline{y_{\tt WIN}}^{(0)}` to their basic :math:`\overline{x}`
and :math:`\overline{y}` isophotal equivalents, respectively. Then at
each iteration :math:`t`, :math:`\overline{x_{\tt WIN}}` and
:math:`\overline{y_{\tt WIN}}` are refined using:

.. math::

   \begin{aligned}
   \label{eq:xwin}
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

.. math:: w_i^{(t)} = \exp \left(-\frac{r_i^{(t)^2}}{2s_{\tt WIN}^2} \right),

with

.. math::

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
the theoretical limit set by image noise [1]_.

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

Windowed second-order moments are computed on the image data once the
centering process from §[chap:wincent] has converged:

.. math::

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

Windowed ellipse parameters: CXXWIN, CYYWIN, CXYWIN
---------------------------------------------------

They are computed from the windowed 2nd order moments exactly the same
way as in §[chap:cxx].

Windowed position uncertainties: ERRX2WIN, ERRY2WIN, ERRXYWIN, ERRAWIN, ERRBWIN, ERRTHETAWIN, ERRCXXWIN, ERRCYYWIN, ERRCXYWIN
-----------------------------------------------------------------------------------------------------------------------------

Windowed position uncertainties are computed on the image data once the
centering process from §[chap:wincent] has converged. Assuming that
noise is uncorrelated among pixels, standard error propagation applied
to ([eq:xwin]) and ([eq:xwin]) gives us:

.. math::

   \begin{aligned}
   {\tt ERRX2WIN} & = {\rm var}(\overline{x_{\tt WIN}})
   = & 4\,\frac{\sum_{r_i < r_{\rm max}} w_i^2 \sigma^2_i (x_i-\overline{x})^2}
   {\left(\sum_{r_i < r_{\rm max}} w_i I_i\right)^2},\\
   {\tt ERRY2WIN} & = {\rm var}(\overline{y_{\tt WIN}})
   = & 4\,\frac{\sum_{r_i < r_{\rm max}} w_i^2 \sigma^2_i (y_i-\overline{y})^2}
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
exactly as in §[chap:poserr]: see eqs. ([eq:erra]), ([eq:errb]),
([eq:errtheta]), ([eq:errcxx]), ([eq:errcyy]) and ([eq:errcxy]).

.. [1]
   see http://www.astromatic.net/forum/showthread.php?tid=581

.. include:: keys.rst

