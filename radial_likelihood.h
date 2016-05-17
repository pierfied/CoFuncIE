#ifndef RADIAL_LIKE_H
#define RADIAL_LIKE_H

double *fineMap(double *map, int numVoxelsPerDim, int resoultion);
double *deltaMap(double *map, int numVoxelsPerDim, int dim);
double trilInterp(double *map, int numVoxelsPerDim, int resoultion,
	int x, int y, int z);

#endif