#
#				Makefile.am
#
# Main Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SExtractor
#
#	Copyright:		(C) 2002-2017 IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SExtractor is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SExtractor is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		19/06/2017
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACLOCAL_AMFLAGS		= -I m4
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= src man doc tests debian
dist_pkgdata_DATA	= config/* xsl/sextractor.xsl
EXTRA_DIST		= autogen.sh AUTHORS BUGS ChangeLog COPYRIGHT \
			  HISTORY INSTALL LICENSE README.md THANKS
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
DEB_SRCDIR		= ~/DEBIAN/SOURCES

dist-hook:
	rm -rf `find $(distdir) -type d \( -name .svn -o -name .git \)`

rpm:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-best:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_BEST="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

deb:	dist
	mkdir -p $(DEB_SRCDIR)
	tar -xvf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz -C $(DEB_SRCDIR) 
	cd $(DEB_SRCDIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)
	debuild -j4 -b -us -uc

deb-best:	dist
	mkdir -p $(DEB_SRCDIR)
	tar -xvf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz -C $(DEB_SRCDIR) 
	cd $(DEB_SRCDIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)
	debuild --set-envvar=USE_BEST=1 --preserve-env --preserve-envvar=PATH -j4 -b -us -uc

