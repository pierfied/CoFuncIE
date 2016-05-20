#include "radial_likelihood.h"

#define PI 3.14159265359

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

double drawRSample(galaxy *gal, mapData md){
	// Create a fine version of the density contrast map using trilinear
	// interpolation.
	int resolution = 10;
	double *fMap = fineMap(md.map, md.numVoxelsX, resolution);

	// Convert the galaxy of interest into cartesian coordinates and calculate
	// the radial distance.
	cartesianGalaxy *cGal = convertToCartesian(gal, 1);
	double dist = sqrt(cGal->x*cGal->x + cGal->y*cGal->y + cGal->z*cGal->z);

	// Get radial distance error from the photo-z.
	double sigR = cGal->r_err;

	// Calculate the unit vector components of the galaxy position.
	double xUnit = cGal->x/dist;
	double yUnit = cGal->y/dist;
	double zUnit = cGal->z/dist;

	// Get the offsets of the box.
	int xStart = md.xStart;
	int yStart = md.yStart;
	int zStart = md.zStart;

	// Get the number of voxels in each direction.
	int numVoxelsX = md.numVoxelsX;
	int numVoxelsY = md.numVoxelsY;
	int numVoxelsZ = md.numVoxelsZ;

	// Get the length of the box in each direction.
	int boxLengthX = md.boxLengthX;
	int boxLengthY = md.boxLengthY;
	int boxLengthZ = md.boxLengthZ;

	// Calculate the length of each voxel (should be the same FIX!!!!).
	double voxelLengthX = boxLengthX/(double)(numVoxelsX*resolution);
	double voxelLengthY = boxLengthY/(double)(numVoxelsY*resolution);
	double voxelLengthZ = boxLengthZ/(double)(numVoxelsZ*resolution);

	// Calculate the number of voxels in each direction in the fine map.
	double fineNVPDX = numVoxelsX*resolution;
	double fineNVPDY = numVoxelsY*resolution;
	double fineNVPDZ = numVoxelsZ*resolution;

	// Set the maximum radial distance to calculate out to.
	int maxR = 1500;

	// Allocate the pdf and cdf arrays.
	double pdf[maxR];
	double cdf[maxR];

	FILE *fp;
	fp = fopen("rlike.txt","w");

	// Initialize the normalization constant for the pdf.
	double normConst = 0;

	// Loop over all radial distances of interest.
	int r;
	for(r = 0; r < maxR; r++){
		// Caluclate the x, y, and z values for the given radius.
		double x = r * xUnit;
		double y = r * yUnit;
		double z = r * zUnit;

		// Calculate the gaussian term from the photo-z.
		double gaussR = 1/(sigR * sqrt(2*PI)) * exp(-0.5*pow((r-dist)/sigR,2));

		// Check whether the point of interest is within the box or not.
		if(x < xStart || x > xStart + boxLengthX || y < yStart ||
			y > yStart + boxLengthY || z < zStart || 
			z > zStart + boxLengthZ){

			// Assume that the delta map is zero outside the box. Therefore,
			// the pdf is only the gaussian term.
			pdf[r] = gaussR;
		}else{
			// Calculate the index location in the fine map.
			int a = (x - xStart)/voxelLengthX;
			int b = (y - yStart)/voxelLengthY;
			int c = (z - zStart)/voxelLengthZ;
			int index = a + fineNVPDY*b + pow(fineNVPDZ, 2)*c;

			// Calculate the map term (1+delta)
			double mapTerm = 1 + fMap[index];

			// Calculate the pdf at the given radius.
			pdf[r] = gaussR * mapTerm;
		}

		// Add to the normalization constant sum.
		normConst += pdf[r];
	}

	// Declare the radius array.
	double y[maxR];
	double prevCDF = 0;
	for(r = 0; r < maxR; r++){
		// Normalize the current pdf value and set the cdf at the current r.
		pdf[r] /= normConst;
		cdf[r] = prevCDF + pdf[r];

		// Store the cdf at the current r.
		prevCDF = cdf[r];

		// Simply store the current r in the radius array (used for
		// drawing from the cdf later).
		y[r] = r;

		fprintf(fp,"%d,%f,%f\n",r,pdf[r],cdf[r]);
	}

	printf("Photo-z R: %f\n", dist);

	int i;
	for(i = 0; i < 10; i++){
		double rSamp = drawFromCDF(cdf,y,maxR);
		printf("Random R Sample: %f\n", rSamp);
	}

	exit(0);

	return 0;
}

double drawFromCDF(double cdf[], double x[], int numPts){
	// Generate a random number between zero and one.
	double rnd = ((double)rand()/(double)RAND_MAX);;

	// Loop through the CDF values.
	int i;
	for(i = 0; i < numPts-1; i++){
		// Check if the random value is between the current cdf value
		// and the next.
		if(cdf[i] < rnd && rnd < cdf[i+1]){
			// Perform linear interpolation between the points and return
			// the corresponding x value.
			return x[i] + (x[i+1]-x[i])*(rnd-cdf[i])/(cdf[i+1]-cdf[i]);
		}
	}

	// If no value was found within the bounds, use the final point.
	return x[numPts-1];
}