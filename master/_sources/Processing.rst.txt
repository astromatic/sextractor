.. File Processing.rst

.. include:: global.rst

.. _processing:

Processing
==========

The complete analysis of an image is fully automated (:numref:`fig_sexlayout`). Two passes are made through the data.
During the first pass, a `model of the sky background <background_mode>`_ is built, and several global statistics are computed.
During the second pass, image pixels are background-subtracted, filtered and segmented on-the-fly.
Detections are then deblended, pruned (“CLEANed”), and enter the measurement phase.
Finally, the measured quantities are written to the output catalog, after cross-matching with an optional input list.

.. _fig_sexlayout:

.. figure:: figures/sexlayout.*
   :align: center

   Layout of the main |SExtractor| procedures. *Dashed arrows* represent
   optional inputs. 

The following sections describe each of these operations in more detail.

.. toctree::

   Background
   Weighting
   Flagging

