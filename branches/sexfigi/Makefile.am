# Main Makefile for SExtractor
# Copyright (C) 2002-2007 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src tests
EXTRA_DIST		= config doc AUTHORS BUGS ChangeLog COPYING HISTORY \
			  INSTALL README THANKS \
			  acx_prog_cc_optim.m4
RPM_ROOTDIR		= /usr/src/RPM
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O -g -funroll-loops -fomit-frame-pointer -Wall \
		-march=opteron" TPXLIBS="-lfftw3 -lm" \
		rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec

	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O -g -pthread -funroll-loops \
		-fomit-frame-pointer -Wall -march=opteron" \
		TPXLIBS="-lpthread -lfftw3 -lm" \
		rpmbuild -ba --target=x86_64 --clean  $(PACKAGE_NAME)-mp.spec

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O -g -funroll-loops -fomit-frame-pointer -Wall \
		-m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" \
		TPXLIBS="-L/usr/lib/ -lm -L/usr/local/lib/ -lfftw3" \
		rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec

	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O -g -pthread -funroll-loops \
		-fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse \
		-mtune=athlon" \
		TPXLIBS="-L/usr/lib/ -lpthread -lm -L/usr/local/lib/ -lfftw3" \
		rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME)-mp.spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -g -axWP -ip -unroll" TPXLIBS="-lfftw3 -lm" \
		rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec

	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -g -axWP -ip -unroll -pthread" \
	TPXLIBS="-lpthread -lfftw3 -lm" \
	rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME)-mp.spec

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -g -axKWNPB -ip -unroll" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lfftw3 -lm" \
	rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec

	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -g -axKWNPB -ip -unroll -pthread" \
	TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lpthread -lfftw3f -lm" \
	rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME)-mp.spec

athlon:
	$(MAKE) CFLAGS="-O -g -pthread -funroll-loops -fomit-frame-pointer \
		-Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" \
		LIBS="-L/usr/lib/ -L/usr/local/lib/ -lfftw3 -lm"

icc:
	$(MAKE) CC="icc" CFLAGS="-O3 -g -axWP -ip -unroll -pthread"
