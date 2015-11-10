/*
  This file contains the assumed cosmology we will use.
  It is in a separated file for easy access.

  This header file contains the prototypes for
  cosmology.c
  See accompanying .c file for more information.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef COSMOLOGY_H

#define cosmo_M 0.3
#define cosmo_k 0.0
#define cosmo_lambda 0.7
#define cosmo_h 0.7

double comovingDistance(double z);
void setupRedshiftInterp(double* *r, double* *z, int* numInterpPts);
double interpRedshift(double rQuery, double* r, double* z, int numInterpPts);
double interpDist(double zQuery, double* r, double* z, int numInterpPts);

#endif
