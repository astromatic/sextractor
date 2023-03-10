.. File Using.rst

.. include:: global.rst

.. _using-sextractor:

Using SExtractor
================

|SExtractor| is run from the shell with the following syntax:

.. code-block:: console

  $ sex Image1 [Image2] -c configuration-file [-Parameter1 Value1 -Parameter2 Value2 ...]

The parts enclosed within brackets are optional.
Any `-Parameter Value` statement in the command-line overrides the corresponding definition in the configuration file or any default value (see :ref:`configuration section<config_file>`).

.. toctree::

   Input
   Output
   Config
   Param

