/*=============================================================================
*
*   WCSLIB - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995-1999, Mark Calabretta
*
*   This library is free software; you can redistribute it and/or modify it
*   under the terms of the GNU Library General Public License as published
*   by the Free Software Foundation; either version 2 of the License, or (at
*   your option) any later version.
*
*   This library is distributed in the hope that it will be useful, but
*   WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library
*   General Public License for more details.
*
*   You should have received a copy of the GNU Library General Public License
*   along with this library; if not, write to the Free Software Foundation,
*   Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*   Correspondence concerning WCSLIB may be directed to:
*      Internet email: mcalabre@atnf.csiro.au
*      Postal address: Dr. Mark Calabretta,
*                      Australia Telescope National Facility,
*                      P.O. Box 76,
*                      Epping, NSW, 2121,
*                      AUSTRALIA
*
*   Author: Mark Calabretta, Australia Telescope National Facility
*   IRAF's TNX added by E.Bertin 2000/03/28
*   $Id: proj.h,v 1.1.1.1 2002/03/15 16:33:26 bertin Exp $
*===========================================================================*/

#ifndef WCSLIB_PROJ
#define WCSLIB_PROJ

#ifdef __cplusplus
extern "C" {
#endif

struct prjprm {
   int flag;
   int n;
   double r0;
   double p[200];
   double w[10];
   struct tnxaxis	*tnx_latcor;
   struct tnxaxis	*tnx_lngcor;
   struct poly		*inv_x;
   struct poly		*inv_y;
};

#if __STDC__ || defined(__cplusplus)
   int azpset(struct prjprm *);
   int azpfwd(const double, const double, struct prjprm *, double *, double *);
   int azprev(const double, const double, struct prjprm *, double *, double *);
   int tanset(struct prjprm *);
   int tanfwd(const double, const double, struct prjprm *, double *, double *);
   int tanrev(const double, const double, struct prjprm *, double *, double *);
   int sinset(struct prjprm *);
   int sinfwd(const double, const double, struct prjprm *, double *, double *);
   int sinrev(const double, const double, struct prjprm *, double *, double *);
   int stgset(struct prjprm *);
   int stgfwd(const double, const double, struct prjprm *, double *, double *);
   int stgrev(const double, const double, struct prjprm *, double *, double *);
   int arcset(struct prjprm *);
   int arcfwd(const double, const double, struct prjprm *, double *, double *);
   int arcrev(const double, const double, struct prjprm *, double *, double *);
   int zpnset(struct prjprm *);
   int zpnfwd(const double, const double, struct prjprm *, double *, double *);
   int zpnrev(const double, const double, struct prjprm *, double *, double *);
   int zeaset(struct prjprm *);
   int zeafwd(const double, const double, struct prjprm *, double *, double *);
   int zearev(const double, const double, struct prjprm *, double *, double *);
   int airset(struct prjprm *);
   int airfwd(const double, const double, struct prjprm *, double *, double *);
   int airrev(const double, const double, struct prjprm *, double *, double *);
   int cypset(struct prjprm *);
   int cypfwd(const double, const double, struct prjprm *, double *, double *);
   int cyprev(const double, const double, struct prjprm *, double *, double *);
   int carset(struct prjprm *);
   int carfwd(const double, const double, struct prjprm *, double *, double *);
   int carrev(const double, const double, struct prjprm *, double *, double *);
   int merset(struct prjprm *);
   int merfwd(const double, const double, struct prjprm *, double *, double *);
   int merrev(const double, const double, struct prjprm *, double *, double *);
   int ceaset(struct prjprm *);
   int ceafwd(const double, const double, struct prjprm *, double *, double *);
   int cearev(const double, const double, struct prjprm *, double *, double *);
   int copset(struct prjprm *);
   int copfwd(const double, const double, struct prjprm *, double *, double *);
   int coprev(const double, const double, struct prjprm *, double *, double *);
   int codset(struct prjprm *);
   int codfwd(const double, const double, struct prjprm *, double *, double *);
   int codrev(const double, const double, struct prjprm *, double *, double *);
   int coeset(struct prjprm *);
   int coefwd(const double, const double, struct prjprm *, double *, double *);
   int coerev(const double, const double, struct prjprm *, double *, double *);
   int cooset(struct prjprm *);
   int coofwd(const double, const double, struct prjprm *, double *, double *);
   int coorev(const double, const double, struct prjprm *, double *, double *);
   int bonset(struct prjprm *);
   int bonfwd(const double, const double, struct prjprm *, double *, double *);
   int bonrev(const double, const double, struct prjprm *, double *, double *);
   int pcoset(struct prjprm *);
   int pcofwd(const double, const double, struct prjprm *, double *, double *);
   int pcorev(const double, const double, struct prjprm *, double *, double *);
   int glsset(struct prjprm *);
   int glsfwd(const double, const double, struct prjprm *, double *, double *);
   int glsrev(const double, const double, struct prjprm *, double *, double *);
   int parset(struct prjprm *);
   int parfwd(const double, const double, struct prjprm *, double *, double *);
   int parrev(const double, const double, struct prjprm *, double *, double *);
   int aitset(struct prjprm *);
   int aitfwd(const double, const double, struct prjprm *, double *, double *);
   int aitrev(const double, const double, struct prjprm *, double *, double *);
   int molset(struct prjprm *);
   int molfwd(const double, const double, struct prjprm *, double *, double *);
   int molrev(const double, const double, struct prjprm *, double *, double *);
   int cscset(struct prjprm *);
   int cscfwd(const double, const double, struct prjprm *, double *, double *);
   int cscrev(const double, const double, struct prjprm *, double *, double *);
   int qscset(struct prjprm *);
   int qscfwd(const double, const double, struct prjprm *, double *, double *);
   int qscrev(const double, const double, struct prjprm *, double *, double *);
   int tscset(struct prjprm *);
   int tscfwd(const double, const double, struct prjprm *, double *, double *);
   int tscrev(const double, const double, struct prjprm *, double *, double *);
   int tnxset(struct prjprm *);
   int tnxfwd(const double, const double, struct prjprm *, double *, double *);
   int tnxrev(const double, const double, struct prjprm *, double *, double *);
   int raw_to_pv(struct prjprm *, double, double, double *, double *);
#else
   int azpset(), azpfwd(), azprev();
   int tanset(), tanfwd(), tanrev();
   int sinset(), sinfwd(), sinrev();
   int stgset(), stgfwd(), stgrev();
   int arcset(), arcfwd(), arcrev();
   int zpnset(), zpnfwd(), zpnrev();
   int zeaset(), zeafwd(), zearev();
   int airset(), airfwd(), airrev();
   int cypset(), cypfwd(), cyprev();
   int carset(), carfwd(), carrev();
   int merset(), merfwd(), merrev();
   int ceaset(), ceafwd(), cearev();
   int copset(), copfwd(), coprev();
   int codset(), codfwd(), codrev();
   int coeset(), coefwd(), coerev();
   int cooset(), coofwd(), coorev();
   int bonset(), bonfwd(), bonrev();
   int pcoset(), pcofwd(), pcorev();
   int glsset(), glsfwd(), glsrev();
   int parset(), parfwd(), parrev();
   int aitset(), aitfwd(), aitrev();
   int molset(), molfwd(), molrev();
   int cscset(), cscfwd(), cscrev();
   int qscset(), qscfwd(), qscrev();
   int tscset(), tscfwd(), tscrev();
   int tnxset(), tnxfwd(), tnxrev();
#endif
/*
extern const char *prjset_errmsg[];
extern const char *prjfwd_errmsg[];
extern const char *prjrev_errmsg[];
*/
#define PRJSET 137

#ifdef __cplusplus
};
#endif

#endif /* WCSLIB_PROJ */
