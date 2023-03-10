.. File Flagging.rst

.. include:: global.rst

.. _flagging:

Flagging
========

Both *internal* and *external* flags are available for each detection as catalog parameters:

* Several sets of internal flags are produced by the various processes within |SExtractor|:

  * :ref:`Extraction flags <flags_def>` (:param:`FLAGS`)
  * :ref:`Weighting flags <flags_weight_def>` (:param:`FLAGS_WEIGHT`)
  * :ref:`Windowed measurement flags <flags_win_def>` (:param:`FLAGS_WIN`)
  * :ref:`Model-fitting flags <flags_model_def>` (:param:`FLAGS_MODEL`);

* :ref:`External flags <imaflags_iso_def>` are derived from *flag maps*. Flag maps are images in integer format having the same dimensions as the science images, with pixel values that can be used to flag some pixels (for instance, "bad" or noisy pixels). Different combinations of flags can be applied within the isophotal area that defines each object, to produce a unique value that will be written to the catalog.

Each catalog parameter comprises several flag bits as a sum of powers of 2 coded in decimal. For example, a saturated detection close to an image boundary will have :param:`FLAGS` = 4+8 = 12 (see below).

.. _flags_def:

Extraction flags: :param:`FLAGS`
--------------------------------

:param:`FLAGS` contains 8 flag bits with basic warnings about the source extraction process, in order of increasing concern.

.. _flags_table:

.. csv-table:: :param:`FLAGS` description
  :header: "Value", "Meaning"
  :widths: 3 60

  1, ":ref:`aperture photometry <photometry>` is likely to be biased by neighboring sources or by more than 10% of bad pixels in any aperture"
  2, "the object has been deblended"
  4, "at least one object pixel is saturated"
  8, "the isophotal footprint of the detected object is truncated (too close to an image boundary)"
  16, "at least one photometric aperture is incomplete or corrupted (hitting buffer or memory limits)"
  32, "the isophotal footprint is incomplete or corrupted (hitting buffer or memory limits)"
  64, "a memory overflow occurred during deblending"
  128, "a memory overflow occurred during extraction"

.. _imaflags_iso_def:

External flags: :param:`IMAFLAGS_ISO`
-------------------------------------

The processing of external flags is triggered whenever :param:`IMAFLAGS_ISO` or :param:`NIMAFLAGS_ISO` are present in the catalog parameter file.
The file names of the images containing flag maps in |FITS| format are set with the ``FLAG_IMAGE`` configuration keyword.
Flag maps must be stored as 2-dimensional arrays of 8, 16 or 32 bits integers.
The integers may be signed or unsigned, although for Boolean operations the sign bit will be interpreted as bit 7, 15 or 31, depending on integer size.
Obviously all flag maps must have the same dimensions as the image used for detection. 
Flag maps can be created using, e.g., the |WeightWatcher|_ software :cite:`2008ASPC_394_619M`.

The values of flag map pixels that overlap the isophotal area of a given detected object are combined and stored in the catalog as the long integer :param:`IMAFLAGS_ISO`.
The ``FLAG_TYPE`` configuration keyword sets the type of combination applied individually to each flag map.
5 types of combination are currently supported:

- ``OR``: arithmetic (bit-to-bit) **OR** of flag map pixels.
- ``AND``: arithmetic (bit-to-bit) **AND** of non-zero flag map pixels.
- ``MIN``: minimum of the (signed) flag map pixel values.
- ``MAX``: maximum of the (signed) flag map pixel values.
- ``MOST``: most frequent non-zero flag map pixel value.

The :param:`NIMAFLAGS_ISO` catalog parameter contains the number(s) of relevant flag map pixels: the number of non-zero flag map pixels if ``FLAG_TYPE`` is set to ``OR`` or ``AND``, or the number of pixels with value :param:`IMAFLAGS_ISO` if ``FLAG_TYPE`` is set to ``MIN`` , ``MAX``, or ``MOST``.
