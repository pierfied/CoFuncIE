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
	int resolution = 10;
	double *fMap = fineMap(md.map, md.numVoxelsX, resolution);

	cartesianGalaxy *cGal = convertToCartesian(gal, 1);
	double dist = sqrt(cGal->x*cGal->x + cGal->y*cGal->y + cGal->z*cGal->z);

	double sigR = cGal->r_err;

	double xUnit = cGal->x/dist;
	double yUnit = cGal->y/dist;
	double zUnit = cGal->z/dist;

	int xStart = md.xStart;
	int yStart = md.yStart;
	int zStart = md.zStart;

	int numVoxelsX = md.numVoxelsX;
	int numVoxelsY = md.numVoxelsY;
	int numVoxelsZ = md.numVoxelsZ;

	int boxLengthX = md.boxLengthX;
	int boxLengthY = md.boxLengthY;
	int boxLengthZ = md.boxLengthZ;

	double voxelLengthX = boxLengthX/(double)(numVoxelsX*resolution);
	double voxelLengthY = boxLengthY/(double)(numVoxelsY*resolution);
	double voxelLengthZ = boxLengthZ/(double)(numVoxelsZ*resolution);

	double fineNVPDX = numVoxelsX*resolution;
	double fineNVPDY = numVoxelsY*resolution;
	double fineNVPDZ = numVoxelsZ*resolution;

	int maxR = 1500;

	double pdf[maxR];
	double cdf[maxR];

	FILE *fp;
	fp = fopen("rlike.txt","w");

	double normConst = 0;

	int r;
	for(r = 0; r < maxR; r++){
		double x = r * xUnit;
		double y = r * yUnit;
		double z = r * zUnit;

		double gaussR = 1/(sigR * sqrt(2*PI)) * exp(-0.5*pow((r-dist)/sigR,2));

		if(x < xStart || x > xStart + boxLengthX || y < yStart ||
			y > yStart + boxLengthY || z < zStart || 
			z > zStart + boxLengthZ){

			pdf[r] = gaussR;
		}else{
			int a = (x - xStart)/voxelLengthX;
			int b = (y - yStart)/voxelLengthY;
			int c = (z - zStart)/voxelLengthZ;
			int index = a + fineNVPDY*b + pow(fineNVPDZ, 2)*c;

			double mapTerm = 1 + fMap[index];

			pdf[r] = gaussR * mapTerm;
		}

		normConst += pdf[r];
	}

	double y[maxR];
	double prevCDF = 0;
	for(r = 0; r < maxR; r++){
		pdf[r] /= normConst;
		cdf[r] = prevCDF + pdf[r];

		prevCDF = cdf[r];

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
	double rnd = ((double)rand()/(double)RAND_MAX);;

	int i;
	for(i = 0; i < numPts-1; i++){
		if(cdf[i] < rnd && rnd < cdf[i+1]){
			return x[i] + (x[i+1]-x[i])*(rnd-cdf[i])/(cdf[i+1]-cdf[i]);
		}
	}

	return x[numPts-1];
}