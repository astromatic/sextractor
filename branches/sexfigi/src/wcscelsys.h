/*
 				wcscelsys.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	LDACTools+
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for fitswcs.c
*
*	Last modify:	14/01/2001
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*-------------------------------- constants --------------------------------*/

/* equatorial coordinates of origin and pole and rotation sign of equatorial,*/
/* galactic, and ecliptic reference frames, from Allen Ast. Quant., 4th ed. */
char	celsysname[][2][8] = {  {"RA--", "DEC-"},
				{"GLON", "GLAT"},
				{"ELON", "ELAT"},
				{""}};
double	celsysorig[][2] = {	{0.0, 0.0},
				{266.40499625, -28.93617242},
				{0.0, 0.0} },
	celsyspole[][2] = {	{0.0, 90.0},
				{192.85948123, 27.12825120},
				{273.85261111, 66.99111111}},
/* Note: the code to handle the rotation sign is not yet implemented!!! */
	celsyssign[]	= {	 1.0,
				 1.0,
				 1.0};

