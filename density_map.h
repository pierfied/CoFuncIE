#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data_handler.h"

#ifndef DENSITY_MAP_H
#define DENSITY_MAP_H

double *generateMap(galaxy *gals, int numGals, int numVoxelsPerDim, int xStart,
	int yStart, int zStart, int boxLength, int **voxels);

typedef struct{
	double *map;

	int xStart;
	int yStart;
	int zStart;

	int numVoxelsX;
	int numVoxelsY;
	int numVoxelsZ;

	int boxLengthX;
	int boxLengthY;
	int boxLengthZ;
} mapData;

#endif
