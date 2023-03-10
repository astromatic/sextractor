.. File ClassStar.rst

.. include:: global.rst

.. _class_star_def:

:param:`CLASS_STAR` classifier
==============================

.. note::
  The :param:`CLASS_STAR` classifier has been superseded by the :param:`SPREAD_MODEL` estimator (see :ref:`spread_model_def`), which offers better performance by making explicit use of the full, variable |PSF| model.

A good discrimination between stars and galaxies is essential for both galactic and extragalactic statistical studies.
The common assumption is that galaxy images look more extended or fuzzier than those of stars (or |QSO|\ s).
|SExtractor| provides the :param:`CLASS_STAR` catalog parameter for separating both types of sources.
The  :param:`CLASS_STAR` classifier relies on a `multilayer feed-forward neural network <https://en.wikipedia.org/wiki/Multilayer_perceptron>`_ trained using `supervised learning <https://en.wikipedia.org/wiki/Supervised_learning>`_ to estimate the *a posteriori* probability :cite:`Richard1991,Saerens2002` of a |SExtractor| detection to be a point source or an extended object.
Below is a shortened description of the estimator, see :cite:`1996AAS_117_393B` for more details.

Inputs and outputs
------------------

The neural network is a multilayer Perceptron with a single fully connected, hidden layers.
Of all neural networks it is probably the best-studied, and it has been intensively applied with success for many classification tasks.

The classifier (:numref:`fig_classstarnn`) has 10 inputs:

* 8 isophotal areas :math:`A_0..A_7`, measured at isophotes exponentially spaced between the analysis threshold (which may be modified with the ``ANALYSIS_THRESH`` configuration parameter) and the object's peak pixel value
* The object's peak pixel value above the local background :math:`I_{\mathrm max}`
* A |seeing| input, which acts as a tuning button.

The output layer contains only one neuron, as "star" and "galaxy" are two classes mutually exclusive.
The output value is a "stellarity index", which for images that reasonably match those of the training sample is an estimation of the *a posteriori* probability for the classified object to be a point-source.
Hence a :param:`CLASS_STAR` close to 0 means that the object is very likely a galaxy, and 1 that it is a star.
In practice, real data always differ at least slightly from the training sample, and the :param:`CLASS_STAR` output is often a poor approximation of the expected *a posteriori* probabilities.
Nevertheless, a :param:`CLASS_STAR` value closer to 0 or 1 normally indicates a higher confidence in the classification, and the balance between sample completeness and purity may still be adjusted by tweaking the decision threshold .

.. _fig_classstarnn:

.. figure:: figures/classstarnn.*
  :figwidth: 100%
  :align: center

  Architecture of |SExtractor|'s :param:`CLASS_STAR` classifier

The |seeing| input must be set by the user with the ``SEEING_FWHM`` configuration parameter.
If ``SEEING_FWHM`` is set to 0,  it is automatically measured on the |PSF| model which must be provided (using the ``PSF_NAME`` configuration parameter).

If no |PSF| model is available, the ``SEEING_FWHM`` configuration parameter must be adjusted by the user to match the actual average |PSF| |FWHM| on the image.
The accuracy with which ``SEEING_FWHM`` must be set for optimal results ranges from :math:`\pm 20\%` for bright sources to about :math:`\pm 5\%` for the faintest (:numref:`fig_classstar_seeing`). ``SEEING_FWHM`` is expressed in arcseconds.
The ``PIXEL_SCALE`` configuration parameter must therefore also be set by the user if |WCS| information is missing from the |FITS| image header.
There are several ways to measure, directly or indirectly, the size of point sources in |SExtractor|; they may lead to slightly discordant results, depending on the exact shape of the |PSF|.
The measurement :param:`FWHM_IMAGE` (although not the most reliable as an image quality estimator) sets the reference when it comes to setting ``SEEING_FWHM``.

One may check that the ``SEEING_FWHM`` is set correctly by making sure that the typical :param:`CLASS_STAR` value of unclassifiable sources at the faint end of the catalog hovers around the 0.5 mark.

.. _fig_classstar_seeing:

.. figure:: figures/classstar_seeing.*
  :figwidth: 100%
  :align: center

  Architecture of |SExtractor|'s :param:`CLASS_STAR` classifier


Training
--------

This section gives some insight on how the :param:`CLASS_STAR` classifier has been trained.
The main issue with supervised machine learning is the labeling of the large training sample.
Hopefully a big percentage of contemporary astronomical images share a common set of generic features, which can be simulated with sufficient realism to create a large training sample together with the ground truth (labels).
The :param:`CLASS_STAR` classifier was trained on such a sample of artificial images.

Six hundred :math:`512\times512` simulation images containing stars and galaxies were generated to train the network using an early prototype of the |SkyMaker|_ package :cite:`2009MmSAI_80_422B`.
They were done in the blue band, where galaxies present very diversified aspects.
The two parameters governing the shape of the |PSF| (|seeing| |FWHM| and Moffat :math:`\beta` parameter :cite:`1969AA_3_455M`) were chosen randomly with :math:`0.025\le` FWHM :math:`\le 5.5` and :math:`2\le\beta\le4`. Note that the `Moffat function <https://en.wikipedia.org/wiki/Moffat_distribution>`_ used in the simulation is a poor approximation to diffraction-limited images, hence the :param:`CLASS_STAR` classifier is not optimized for space-based images.
The pixel scale was always taken less than :math:`\approx 0.7` FWHM to warrant correct sampling of the image.
Bright galaxies are simply too rare in the sky to constitute a significant training sample on such a small number of simulations.
So, keeping a constant comoving number density, we increased artificially the number of nearby galaxies by making the volume element proportional to :math:`zdz`.
Stars were given a number-magnitude distribution identical to that of galaxies.
**Therefore any pattern presented to the network had a 50% chance to correspond to a star or a galaxy, irrespective of magnitude** [#faint]_.
Crowding in the simulated images was higher than what one sees on real images of high galactic latitude fields, allowing for the presence of many “difficult cases” (close double stars, truncated profiles, etc...) that the neural network classifier had to deal with.

|SExtractor| was run on each image with 8 different extraction thresholds. This led to a catalog with about :math:`10^6` entries with the 10 input parameters plus the class label. Back-propagation learning took about 15 min on a SUN SPARC20 workstation. The final set of synaptic weights was saved to the file :file:`default.nnw` , ready to be used in “feed-forward” mode during source extraction.

.. [#faint]
   Faint galaxies have less chance being detected than faint stars, but it has little effect because the ones that are lost at a given magnitude are predominantly the most extended and consequently the easiest to classify.
