#include <gsl/gsl_interp.h>
#include "cosmology.h"
#include "density_map.h"

#ifndef RADIAL_LIKE_H
#define RADIAL_LIKE_H

double *fineMap(double *map, int numVoxelsPerDim, int resoultion);
double trilInterp(double *map, int numVoxelsPerDim, int resoultion,
	int x, int y, int z);
double drawRSample(galaxy *gal, mapData md, double *fMap, int resolution, int n);
double drawFromCDF(double cdf[], double x[], int numPts);
void drawGalRSamps(galaxy *gals, int numGals, mapData md);

double drawFromCDFPDF(double cdf[], double pdf[], double curR, double x[], int numPts);

#endif