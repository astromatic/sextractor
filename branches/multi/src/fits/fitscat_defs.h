/*
*				fitscat_defs.h
*
* Internal definitions for the LDACTools FITS library.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1995-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		15/02/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"2.0"
#define		DATE		"2003-x-x"
#define		HAVE_SYS_MMAN_H	1
#endif

/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"LDACTools"
#define		COPYRIGHT	"Emmanuel BERTIN (bertin@iap.fr)"
#define		INSTITUTE	"IAP/Leiden"


/*----------------------------- Internal constants --------------------------*/

#define	MAXCHAR		512		/* max. number of characters */
#define	OUTPUT		stdout		/* where all msgs are sent */
#define	KBYTE		1024		/* 1 kbyte! */
#define	MBYTE		(1024*KBYTE)	/* 1 Mbyte! */
#define	GBYTE		(1024*MBYTE)	/* 1 Gbyte! */
#define	DATA_BUFSIZE	(4*MBYTE)	/* data buffer size for I/O's */
#define	BODY_DEFRAM	(256*MBYTE)	/* a fair number by 1999 standards */
#define	BODY_DEFVRAM	(1.9*GBYTE)	/* a fair number by 1999 standards */
#define	BODY_DEFSWAPDIR	"/tmp"		/* OK at least for Unix systems */

#define	BIG		1e+30		/* a huge number */
#define	TINY		(1.0/BIG)	/* a tiny number */
#ifndef PI
#define	PI		3.14159265359	/* never met before? */
#endif

/* NOTES:
We must have:		MAXCHARS >= 16
			DATA_BUFSIZE >= 2 although DATA_BUFSIZE >= 100000
					  is better!!
*/

/*--------------------- in case of missing constants ------------------------*/

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef	EXIT_SUCCESS
#define	EXIT_SUCCESS		0
#endif
#ifndef	EXIT_FAILURE
#define	EXIT_FAILURE		-1
#endif

/*--------------------------------- typedefs --------------------------------*/
typedef	unsigned char	BYTE;			/* a byte */
typedef	int		LONG;			/* for DEC-Alpha... */
	
/*----------------------------- Internal constants --------------------------*/
char		gstr[MAXCHAR];

/*----------------------------- External constants --------------------------*/

extern int	bswapflag;		/* != 0 if bytes are swapped/IEEE */

/*------------------------------- Other Macros -----------------------------*/

#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif

#define QFREAD(ptr, size, file, fname) \
		{if (fread(ptr, (size_t)(size), (size_t)1, file)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname);;}

#define QFWRITE(ptr, size, file, fname) \
		{if (fwrite(ptr, (size_t)(size), (size_t)1, file)!=1) \
		   error(EXIT_FAILURE, "*Error* while writing ", fname);;}

#define	QFSEEK(file, offset, pos, fname) \
		{if (FSEEKO(file, offset, pos)) \
		   error(EXIT_FAILURE,"*Error*: File positioning failed in ", \
			fname);;}

#define	QFTELL(file, pos, fname) \
		{if ((pos=FTELLO(file))==-1) \
		   error(EXIT_FAILURE,"*Error*: File position unknown in ", \
			fname);;}


#define	QFREE(x)	{free(x); x = NULL;}

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))


#define	QPRINTF		if (qflag) fprintf

#define	QFPRINTF(w,x)	{if (qflag) \
				fprintf(w, "\33[1M> %s\n\33[1A",x);;}


#define	QGETKEY(tab, key, keyname, dest) \
	{if (!(key = name_to_key(tab, keyname))) \
	   error(EXIT_FAILURE, "*Error*: No such parameter in catalog: ", \
			keyname); \
	 dest = key->ptr;}

#define MIN(a,b) (a<b?a:b)
