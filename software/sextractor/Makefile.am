# Main Makefile for SExtractor
# Copyright (C) 2002-2007 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
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
		-march=opteron" rpmbuild -ba --target=x86_64 \
		--clean $(PACKAGE_NAME).spec                                                             

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O -g -funroll-loops -fomit-frame-pointer -Wall -m32 \
		-march=i686 -msse -mfpmath=sse -mtune=athlon" rpmbuild \
		-ba --target=i686 --clean $(PACKAGE_NAME).spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -g -axWP -ip -unroll" rpmbuild -ba --target=x86_64 \
		--clean $(PACKAGE_NAME).spec                                                             

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -g -axKWNPB -ip -unroll" rpmbuild \
		-ba --target=i686 --clean $(PACKAGE_NAME).spec

icc:
	$(MAKE) CC="icc" CFLAGS="-O3 -g -axWP -ip -unroll" LIBS="-lm"

athlon:
	$(MAKE) CFLAGS="-O -g -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" LIBS="-lm"

debug:
	$(MAKE) CFLAGS="-O3 -funroll-loops -fomit-frame-pointer -Wall -m32 -g"
