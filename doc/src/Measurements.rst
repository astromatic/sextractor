Measurements
============

Once sources have been detected and deblended, they enter the
measurement phase. |SExtractor| performs two categories of measurements.
Measurements from the first category are made on the isophotal object
profiles. Only pixels above the detection threshold are considered. Many
of these isophotal measurements (like ``X_IMAGE``, ``Y_IMAGE``, etc.) are
necessary for the internal operations of |SExtractor| and are therefore
executed even if they are not requested. Measurements from the second
category have access to all pixels of the image. These measurements are
generally more sophisticated and are done at a later stage of the
processing (after CLEANing and MASKing).

.. toctree::
   :numbered:
   :maxdepth: 2

  Position
  PositionWin

.. include:: keys.rst

