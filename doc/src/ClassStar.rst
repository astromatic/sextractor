.. File ClassStar.rst

.. include:: global.rst

Star-galaxy classifier
======================

A good discrimination between stars and galaxies is essential for both galactic and extragalactic statistical studies.
The common assumption is that galaxy images look more extended or fuzzier than those of stars (or |QSO|\ s).
The :param:`CLASS_STAR` classifier relies on a `multilayer feed-forward neural network <https://en.wikipedia.org/wiki/Multilayer_perceptron>`_ trained using `supervised learning <https://en.wikipedia.org/wiki/Supervised_learning>`_ to estimate the *a posteriori* probability :cite:`Richard1991,Saerens2002` of a |SExtractor| detection to be a point source or an extended object.
Below is a shortened description of the estimator, see :cite:`1996AAS_117_393B` for more details.

Inputs and outputs
------------------

The neural network is a multilayer Perceptron with a single fully connected, hidden layers.
Of all neural networks it is probably the best-studied, and it has been intensively applied with success for many classification tasks.

The classifier (:numref:`fig_classstarnn`) has 10 inputs:

* 8 isophotal areas :math:`A_0..A_7`, measured at isophotes exponentially spaced between the analysis threshold (which may be modified with the ``ANALYSIS_THRESH`` configuration parameter) and the object's peak pixel value
* The object's peak pixel value above the local background :math:`I_{\mathrm max}`
* A |seeing| input, which must be set by the user with the ``SEEING_FWHM`` configuration parameter.

The output layer contains only one neuron, as "star" and "galaxy" are two classes mutually exclusive.
The output value is a "stellarity index", which for images that reasonably match those of the training sample is an estimation of the *a posteriori* probability for the classified object to be a point-source.
Hence a :param:`CLASS_STAR` close to 0 means that the object is very likely a galaxy, and 1 that it is a star,

.. _fig_classstarnn:

.. figure:: figures/classstarnn.*
  :figwidth: 100%
  :align: center

  Architecture of |SExtractor|'s :param:`CLASS_STAR` classifier

The ``SEEING_FWHM`` configuration parameter must be adjusted by the user


Training
--------

The main issue with supervised machine learning is the labeling of the large training sample.
Hopefully a big percentage of contemporary astronomical images share a common set of generic features, which can be simulated with sufficient realism to create a large training sample together with the ground truth (labels).
The :param:`CLASS_STAR` classifier was trained on such a sample of artificial images.

Six hundred :math:`512\times512` simulation images containing stars and galaxies were generated to train the network using an early prototype of the |SkyMaker|_ package :cite:`2009MmSAI_80_422B`.
They were done in the blue band, where galaxies present very diversified aspects.
The two parameters governing the shape of the |PSF| (|seeing| |FWHM| and Moffat :math:`\beta` parameter :cite:`1969AA_3_455M`) were chosen randomly with :math:`0.025\le` FWHM :math:`\le 5.5` and :math:`2\le\beta\le4`.
The pixel scale was always taken less than :math:`\approx 0.7` FWHM to warrant correct sampling of the image.
Bright galaxies are simply too rare in the sky to constitute a significant training sample on such a small number of simulations.
So, keeping a constant comoving number density, we increased artificially the number of nearby galaxies by making the volume element proportional to :math:`zdz`.
Stars were given a number-magnitude distribution identical to that of galaxies.
**Therefore any pattern presented to the network has 50% chance to correspond to a star or a galaxy, irrespective of the magnitude** [3]_.
Crowding in the simulated images is higher than what one sees on real images of high galactic latitude fields, allowing for the presence of many “difficult cases” (close double stars, truncated profiles, etc...) that the neural network classifier had to deal with.

|SExtractor| was run on each image with 8 different extraction thresholds. This led to a catalog with about :math:`10^6` entries with the 10 input parameters plus the class label. Back-propagation learning took about 15 min on a SUN SPARC20 workstation. The final set of synaptic weights was saved to the file :file:`default.nnw` , ready to be used in “feed-forward” mode during source extraction.


It has been proven, both theoretically and experimentally (Richard &
Lippmann :raw-latex:`\cite{richard:lippmann}` and references therein),
that the cost function of a multi-layered neural net as ours is
minimized when outputs estimate Bayesian *a posteriori* probabilities.
Although the present classifier does not output exact Bayesian
probabilities (because real data unavoidably differ a bit from simulated
ones), it provides at least a confidence estimate in the classification.
This is obvious on Fig. [fig:classvsmag], where the stellarity index
tends to “degenerate” to intermediate values as magnitude increases and
classification becomes harder. Besides, qualitative estimates of
confidence done by eye are in good agreement with those of the
classifier.

Thus, the problem of contamination mentioned in Sect. [training],
arising when the classification becomes difficult and the two classes
are in unbalanced proportions, can be straightforwardly solved by
adjusting the decision boundary between “stars” and “galaxies” within
the stellarity index range ([0,1]).


In conclusion, if one excludes heavily saturated stars, the
classification of “bright” objects appears to be reliable at least at a
conservative :math:`\approx`\ 95% level, on both simulated and real
images with a good sampling. Figure [fig:ccdfield] shows an example of a
CCD image and its interpretation by SExtractor.

.. [1]
   Doing so, we neglect of course the change in appearance of galaxies
   from UV to NIR, from large to narrow-band filters; but the effect of
   this change is not larger than the one due to the natural spread of
   galaxy types in the B band.

.. [2]
   Optical artifacts like spikes or scattered light rays would certainly
   be worth identifying, as they are found to pollute all survey
   catalogs, especially around bright stars; but the classification
   parameters used here are simply inadequate.

.. [3]
   Faint galaxies have less chance being detected than faint stars, but
   it would have little effect since the ones that are lost at a given
   magnitude are predominantly the most extended and consequently the
   easiest to classify.
