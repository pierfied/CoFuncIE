/*
  This header file contains the structures used to
  contain galaxy information. 

  The struct "galaxy"
  contains the raw data; that is RA, DEC, a redshift.
  Later, the uncertainty in the redshift may be added to 
  this as well.

  The struct "cartesianGalaxy"
  contains a galaxy's cartesian coordinates,
  namely x, y, and z.
*/

#ifndef GALSTRUCTS_H
#define GALSTRUCTS_H

typedef struct{
	double ra;
	double dec;
	double z;
  double z_photo;
	double z_err;
} galaxy;

typedef struct {
	double x;
	double y;
	double z;
  double r_photo;
  double r_err;
} cartesianGalaxy;

#endif

