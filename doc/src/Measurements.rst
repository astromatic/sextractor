.. File Measurements.rst

.. include:: global.rst

Measurements
============

Once sources have been detected and deblended, they enter the measurement phase.
|SExtractor| performs three categories of measurements: isophotal, full, and model-fitting.

.. _isophotal_measurements:

Isophotal
  Measurements are made on the isophotal object footprints, which are defined on the filtered detection image. Only pixels with values above the detection threshold (set with ``DETECT_THRESH``) are considered [#thresh]_, which makes the analysis extremely fast, but obviously strongly dependent on the threshold itself. This is an issue particularly when the amplitude of the bakground noise varies over the image. Many of the isophotal measurements (e.g., :param:`X_IMAGE`, :param:`Y_IMAGE`, :param:`FLUX_ISO`) are necessary for the internal operations of |SExtractor| and are therefore executed even if they are not requested.

Full
  Measurements have access to all pixels of the image. These measurements are generally more sophisticated, less affected by variable biases induced by the detection threshold, and still reasonably fast. They are done at a later stage of the processing, after CLEANing and MASKing.

Model-fitting
  Measurements require PSF models [#psf_models]_. They are often the most accurate and can recover the flux of saturated objects. They are also much slower, allowing typically only a few tens of objects to be processed every second.

.. toctree::

  Position
  PositionWin
  Photom
  Model

.. [#thresh] For some isophotal measurements pixel values also have to exceed the local analysis threshold set with ``ANALYSIS_THRESH``.
.. [#psf_models] PSF models be computed using the |PSFEx|_ package.


