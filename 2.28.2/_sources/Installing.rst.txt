.. File Installing.rst

.. include:: global.rst

***********************
Installing the software
***********************

Hardware requirements
=====================

|SExtractor| runs in (ANSI) text-mode from a shell. A graphical environment
is not necessary to operate the software.

Memory requirements depend on the size of the images to be processed.
Processing a single image should typically require about 100MB of memory.
For large images (hundreds of Mpixels or more), or in double-image / weighted mode, |SExtractor|'s  memory footprint should be around 500MB, and up to 2GB in the worst cases.
Swap-space can be put to contribution, although a strong performance hit is to be expected.

Obtaining |SExtractor|
----------------------

For Linux users, the simplest way to have |SExtractor| up and running is to install the standard binary package the comes with your Linux distribution.
Run, e.g., ``apt-get sextractor`` (on Debian) or ``dnf sextractor`` (Fedora) as root and |SExtractor|, as well as all its dependencies, will automatically be installed.
If you decided to install the package this way you may skip the following and move straight to the :ref:`next section <Using Sextractor>`.

However if |SExtractor| is not available in your distribution, or to obtain the most recent version, the |SExtractor| source package can be downloaded from `the official GitHub repository <https://github.com/astromatic/sextractor>`_ .
One may choose `one of the stable releases <https://github.com/astromatic/sextractor/releases>`_, or for the fearless, `a copy of the current master development branch <https://github.com/astromatic/sextractor/archive/master.zip>`_.

Software requirements
---------------------

|SExtractor| has been developed on `GNU/Linux <http://en.wikipedia.org/wiki/Linux>`_ machines and should compile on any `POSIX <http://en.wikipedia.org/wiki/POSIX>`_-compliant system (this includes |OSX|_ and `Cygwin <http://www.cygwin.com>`_ on |Windows|_, at the price of some difficulties with the configuration), provided that the *development* packages of the following libraries have been installed:

* |ATLAS|_ V3.6 and above [#atlas_install]_,
* |FFTw|_ V3.0 and above [#fftw_install]_, 

On Fedora/Redhat distributions for instance, the development packages above are available as ``atlas-devel`` and ``fftw-devel``.
Note that |ATLAS| and |FFTw| are not necessary if |SExtractor| is linked with |Intel|'s |MKL|_ library.

Installation
------------

To install from the |GitHub| source package, you must first uncompress the archive:

.. code-block:: console

  $ unzip sextractor-<version>.zip

A new directory called :file:`sextractor-<version>` should now appear at the current location on your disk.
Enter the directory and generate the files required by the `autotools <http://en.wikipedia.org/wiki/GNU_Build_System>`_, which the package relies on:

.. code-block:: console

  $ cd sextractor-<version>
  $ sh autogen.sh

A :program:`configure` script is created.
This script has many options, which may be listed with the ``--help`` option:

.. code-block:: console

  $ ./configure --help

No options are required for compiling with the default GNU C compiler (:program:`gcc`) if all the required libraries are installed at their default locations:

.. code-block:: console

  $ ./configure

Compared to :program:`gcc` and the librairies above, the combination of the |Intel| compiler (:program:`icc`) and the |MKL|_ libraries can give the |SExtractor| executable a strong boost in performance, thanks to better vectorized code.
If :program:`icc` and the |MKL| are installed on your system [#geticc]_ , you can take advantage of them using

.. code-block:: console

  $ ./configure --enable-mkl

Additionally, if the |SExtractor| binary is to be run on a different machine that does not have :program:`icc` and the |MKL| installed (e.g., a cluster computing node), you must configure a partially statically linked executable using

.. code-block:: console

  $ ./configure --enable-mkl --enable-auto-flags --enable-best-link

In all cases, |SExtractor| can now be compiled with

.. code-block:: console

  $ make -j

An :file:`src/sex` executable is created. For system-wide installation, run the usual

.. code-block:: console

  $ sudo make install

You may now check that the software is properly installed by simply typing in your shell:

.. code-block:: console

  $ sex

which will return the version number and other basic information (note that some shells require the :program:`rehash` command to be run before making a freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#atlas_install] Use the ``--with-atlas`` and/or ``--with-atlas-incdir`` options of the |SExtractor| :command:`configure` script to specify the |ATLAS| library and include paths if |ATLAS| files are  installed at unusual locations.
.. [#fftw_install] Make sure that |FFTW| has been compiled with :command:`configure` options ``--enable-threads --enable-float``.
.. [#geticc] The Linux versions of the |Intel| compiler and |MKL| are `available for free to academic researchers, students, educators and open source contributors <http://software.intel.com/qualify-for-free-software>`_.

