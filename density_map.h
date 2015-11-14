#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data_handler.h"

#ifndef DENSITY_MAP_H
#define DENSITY_MAP_H

double *generateMap(galaxy *gals, int numGals, int numVoxelsPerDim, int xStart,
	int yStart, int zStart, int boxLength);

#endif
