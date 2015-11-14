#include <stdio.h>
#include "data_handler.h"
#include "galaxy_structs.h"
#include "density_map.h"
#include "map_likelihood.h"

int main(){
	int numGals;
	galaxy *gals = readData(&numGals);

	printf("%d\n", numGals);

	cartesianGalaxy *cartGals = convertToCartesian(gals, numGals);

	int i;/*
	for(i = 0; i < 10; i++){
		printf("x: %e\ty: %e\tz: %e\n", cartGals->x, cartGals->y, cartGals->z);
		cartGals++;
	}

	double maxZ = (gals++)->z_red;
	double minZ = maxZ;
	double avgZ = maxZ;
	for(i = 1; i < numGals; i++){
		if(gals->z_red < minZ) minZ = gals->z_red;
		if(gals->z_red > maxZ) maxZ = gals->z_red;
		avgZ += gals->z_red;
		gals++;
	}

	printf("Max: %e\nMin: %e\nAvg: %e\n", maxZ, minZ, avgZ/numGals);*/

	//galaxy *newGals = convertFromCartesian(cartGals,gals,numGals);

	double *r, *z;
	int numInterpPts;
	setupRedshiftInterp(&r, &z, &numInterpPts);
	printf("Query = %f\n", interpRedshift(500, r, z, numInterpPts));
	printf("Query = %f\n", interpDist(0.033, r, z, numInterpPts));

	int xStart = 600;
	int yStart = 0;
	int zStart = 0;
	int numVoxelsPerDim = 12;
	int boxLength = 400;

	galaxy *trimmedGals = trimGalaxyList(gals, &numGals, xStart, yStart,
		zStart, boxLength);

	double *voxels = generateMap(trimmedGals, numGals, numVoxelsPerDim,
		xStart, yStart, zStart, boxLength);

	for(i = 0; i < numVoxelsPerDim; i++){
		printf("%f\n", *(voxels+i));
	}

	// Declare and initialize the xi samples.
	int numSamps = 5;
	double rSamp[] = {0, 250, 500, 750, 1000};
	double xiSamp[] = {1, 0.7, 0.4, 0.1, 0.01};

	gsl_spline *spline = initCorrSpline(numSamps, rSamp, xiSamp);

	double *cov = generateCov(numVoxelsPerDim, boxLength, spline);

	FILE *fp;
	fp = fopen("Catalogues/cov.txt", "w");
	int n = pow(numVoxelsPerDim, 3);
	for(i = 0; i < n; i++){
		int j;
		for(j = 0; j < n; j++){
			int index = i + n*j;
			fprintf(fp, "%e ", *(cov + index));
		}
		fprintf(fp, "\n");
	}

	double *invCov = invertCov(cov, numVoxelsPerDim);

	fp = fopen("Catalogues/inv_cov.txt", "w");
	n = pow(numVoxelsPerDim, 3);
	for(i = 0; i < n; i++){
		int j;
		for(j = 0; j < n; j++){
			int index = i + n*j;
			fprintf(fp, "%e ", *(cov + index));
		}
		fprintf(fp, "\n");
	}


	/*
	FILE *fp;
	fp = fopen("Catalogues/trimmedGals.csv", "w");

	for(i = 0; i < numGals; i++){
		fprintf(fp, "%e, %e, %e\n", trimmedGals->ra, trimmedGals->dec,
			trimmedGals->z_red);
		trimmedGals++;
	}
	fclose(fp);
	*/
}
