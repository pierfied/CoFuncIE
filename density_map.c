#include "density_map.h"

double *generateMap(galaxy *gals, int numGals, int numVoxelsPerDim, int xStart,
	int yStart, int zStart, int boxLength){

	// Allocate the voxel array.
	double *voxels = calloc(pow(numVoxelsPerDim, 3), sizeof(double));

	// Convert the galaxy list to cartesian coordinates.
	cartesianGalaxy *cartGals = convertToCartesian(gals, numGals);

	// Loop through each galaxy and add it to the corresponding voxel.
	double voxelLength = (double)boxLength/numVoxelsPerDim;
	int ind;
	#pragma omp parallel for
	for(ind = 0; ind < numGals; ind++){
		// Adjust the location of the galaxy.
		double x = (cartGals+ind)->x - xStart;
		double y = (cartGals+ind)->y - yStart;
		double z = (cartGals+ind)->z - zStart;

		// Calculate the index of the galaxy.
		int i = x/voxelLength;
		int j = y/voxelLength;
		int k = z/voxelLength;
		int index = i + numVoxelsPerDim*j + pow(numVoxelsPerDim, 2)*k;

		// Increment the voxel at this location.
		#pragma omp critical
			(*(voxels + index))++;
	}

	// Loop through each voxel and calculate the contrasts.
	double p0 = numGals/pow(numVoxelsPerDim, 3);
	#pragma omp parallel for
	for(ind = 0; ind < (int)pow(numVoxelsPerDim, 3); ind++){
		*(voxels+ind) = (*(voxels+ind) - p0)/p0;
	}

	return voxels;
}