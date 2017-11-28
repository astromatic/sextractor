.. File Introduction.rst

.. include:: global.rst

Introduction
============

|SExtractor|_ (Source-Extractor) is a program that builds a catalog of objects from an astronomical image. It is particularly oriented towards the reduction of large scale galaxy-survey data, but it also performs well on moderately crowded star fields. Its main features are:

* Support for multi-extension FITS (|MEF|_)
* Speed: up to about 50 Mpixel/s or 10,000 sources/s with a 3 GHz processor
* Ability to work with very large images (up to 2Gx2G pixels on 64 bit machines)
  thanks to buffered image access
* Real-time filtering of images to improve detectability
* Robust deblending of overlapping extended objects
* Flexible catalog output of desired parameters only
* Pixel-to-pixel photometry in dual-image mode
* Fast and accurate Point-Spread-Function and galaxy model fitting.
* Handling of weight maps and flag maps.
* Optimum handling of images with variable SNR.
* Built-in catalog cross-identification.
* Special mode for photographic scans.
* |XML|_ |VOTable|_-compliant catalog output.
* |XSLT|_ filter sheet provided for convenient access to metadata from a regular web browser.

