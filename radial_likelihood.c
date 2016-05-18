#include "radial_likelihood.h"

double *fineMap(double *map, int numVoxelsPerDim, int resoultion){
	// Calculate the number of voxels per dimension in the new resolution.
	int newNumVPD = numVoxelsPerDim * resoultion;

	// Allocate the fine map.
	double *fMap = malloc(sizeof(double) * pow(newNumVPD,3));

	// Loop through each x,y,z index in the fine-map.
	int i;
	#pragma omp parallel for
	for(i = 0; i < newNumVPD; i++){
		int j;
		for(j = 0; j < newNumVPD; j++){
			int k;
			for(k = 0; k < newNumVPD; k++){
				// Calculate the index location in the fine map.
				int index = i + j*newNumVPD + k*pow(newNumVPD,2);

				// Calculate the interpolation bounds.
				int minInterpBound = resoultion/2.0;
				int maxInterpBound = newNumVPD - resoultion/2.0;

				// Check that the points of interest are within the
				// interpolation bounds. If so interpolate, otherwise
				// use the value from the coarse map voxel.
				if(i > minInterpBound && j > minInterpBound
					&& k > minInterpBound && i < maxInterpBound
					&& j < maxInterpBound && k < maxInterpBound){

					// Perform a shift as if the interpolation points were
					// at the center of each voxel.
					int x = i - resoultion/2.0;
					int y = j - resoultion/2.0;
					int z = k - resoultion/2.0;

					// Set the fine map at the current location equal to the
					// interpolated value.
					fMap[index] = trilInterp(map, numVoxelsPerDim, resoultion,
						x, y, z);
				}else{
					// Calculate the voxel location in the coarse map.
					int coarseX = i/resoultion;
					int coarseY = j/resoultion;
					int coarseZ = k/resoultion;

					// Calculate the coarse map index loccation.
					int coarseInd = coarseX + coarseY*numVoxelsPerDim
						+ coarseZ*pow(numVoxelsPerDim,2);

					// Set the fine map value equal to the coarse map value.
					fMap[index] = map[coarseInd];
				}
			}
		}
	}

	return fMap;
}

double trilInterp(double *map, int numVoxelsPerDim, int resoultion,
	int x, int y, int z){

	// Calculate the voxel location in the coarse map.
	int coarseX = x/resoultion;
	int coarseY = y/resoultion;
	int coarseZ = z/resoultion;

	// Calculate the relative differences of the query point.
	double xd = x/(double)resoultion - coarseX;
	double yd = y/(double)resoultion - coarseY;
	double zd = z/(double)resoultion - coarseZ;

	// Calculate the index locations in the coarse map of all eight
	// interpolation points of interest.
	int ind000 = coarseX + coarseY*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind001 = coarseX + coarseY*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind010 = coarseX + (coarseY + 1)*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind011 = coarseX + (coarseY + 1)*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind100 = (coarseX + 1) + coarseY*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind101 = (coarseX + 1) + coarseY*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind110 = (coarseX + 1) + (coarseY + 1)*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind111 = (coarseX + 1) + (coarseY + 1)*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);

	// Interpolate along the x axis reducing to four interpolation points.
	double c00 = map[ind000]*(1-xd) + map[ind100]*xd;
	double c01 = map[ind001]*(1-xd) + map[ind101]*xd;
	double c10 = map[ind010]*(1-xd) + map[ind110]*xd;
	double c11 = map[ind011]*(1-xd) + map[ind111]*xd;

	// Interpoalte along the y axis reducing to two interpolation points.
	double c0 = c00*(1-yd) + c10*yd;
	double c1 = c01*(1-yd) + c11*yd;

	// Interpolate along the z axis producing the interpolated value.
	double c = c0*(1-zd) + c1*zd;

	return c;
}