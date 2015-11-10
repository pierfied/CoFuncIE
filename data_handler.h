/*
  This header file contains the prototypes for
  data_handler.c. 
  See accompanying .c file for more information.
*/

#include <stdio.h>
#include <math.h>
#include "cosmology.h"
#include "galaxy_structs.h"

#ifndef HANDLER_H

galaxy* readData(int* numGals);
cartesianGalaxy* convertToCartesian(galaxy* gals, int numGals);
galaxy* convertFromCartesian(cartesianGalaxy* cartGals, galaxy* origGals,int numGals);
galaxy *trimGalaxyList(galaxy *gals, int *numGals);

//TODO: These functions should be moved to the cosmology.c program
double comovingDistance(double z);
void setupRedshiftInterp(double* *r, double* *z, int* numInterpPts);
double interpRedshift(double rQuery, double* r, double* z, int numInterpPts);
double interpDist(double zQuery, double* r, double* z, int numInterpPts);

#endif
