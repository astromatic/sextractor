.. File Input.rst

.. include:: global.rst

.. _input:

Input files
-----------

|SExtractor| accepts images stored in |FITS|_ :cite:`1981AAS_44_363W`.
Both "Basic FITS" (one single header and one single body) and |MEF|_
files are recognized.
Binary |SExtractor| catalogs produced from |MEF| images are |MEF| files themselves.
If the catalog output format is set to ASCII, all catalogs from the individual extensions are concatenated in one big file; the :param:`EXT_NUMBER` catalog parameter can be used to tell which extension the detection belongs to.

.. important::
  Contrary to most other astronomy packages, |SExtractor| does not rely on the `FITSIO <https://heasarc.gsfc.nasa.gov/fitsio/>`_ library and instead uses its own library for managing |FITS| files. As a consequence, some features of |FITS| such as image compression/tiling are not supported at this time.

For images with :math:`{\rm NAXIS} > 2`, only the first data-plane is loaded.
In |SExtractor|, as in all similar programs, |FITS| axis #1 is traditionally
referred to as the *x* axis, and |FITS| axis #2 as the *y* axis.

Double image mode
~~~~~~~~~~~~~~~~~

When the pixel data for a given field are available using the same pixel grid in several photometric channels, it is often best to measure object characteristics, like magnitudes, exactly at the same positions and through the same apertures for the different channels.
This to derive precise color indices for example.
|SExtractor| makes this possible by providing a special mode dubbed "double image mode" where detections are made on one image while measurements are carried out on another. Both images must have the exact same dimensions, obviously.
By repeatedly running |SExtractor| with various "measurement images" while keeping the same "detection image", one ends up with a set of catalogs having the same sources measured through different channels.
The detection image will generally be chosen in the band where the data are the deepest.
Alternatively, using a ":math:`\chi2` image" :cite:`1999AJ_117_68S`  [#swarpchi2]_ as a detection image, will allow most of the sources present in at least one channel to be detected and measured.

Double image mode is automatically engaged when providing |SExtractor| with two images instead of one:

.. code-block:: console

  $ sex detection.fits,measurement.fits

A space may be used instead of a coma. In the example above, :file:`detection.fits` is used as a detection image, while measurements are carried out on :file:`measurement.fits`.

.. [#swarpchi2]
   :math:`\chi2` images can easily be generated using the |SWarp|_ package :cite:`2002ASPC_281_228B`.


