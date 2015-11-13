/*
  This file contains functions used to read in galaxy data 
  from a catalogue, trim the catalogue, and convert it to 
  cartesian coordinates.
*/

#include "data_handler.h"

#define PI 3.14159265359

galaxy *readData(int *numGals){
	// Open the file containing the location of the dataset.
	char *fname = "data_loc";
	FILE *fp;
	fp = fopen(fname, "r");

	// Read in the location of the dataset.
	char loc[100];
	fscanf(fp, "%s", loc);
	fclose(fp);

	// Open the file containing the dataset and get the number of galxies.
	fp = fopen(loc, "r");
	*numGals = 0;
	while(fgets(loc, sizeof(loc), fp) != NULL){
		(*numGals)++;
	}
	rewind(fp);

	// Declare the galaxy list.
	galaxy *gals = malloc(*numGals * sizeof(galaxy));

	// Read in the galaxies.
	int i;
	galaxy *curGal = gals;
	for(i = 0; i < *numGals; i++){
		fscanf(fp, "%le, %le, %le, %le", &(curGal->ra), &(curGal->dec),
			&(curGal->z_red), &(curGal->z_err));
		curGal++;
	}
	fclose(fp);

	// Return the list of galaxies.
	return gals;
}

cartesianGalaxy *convertToCartesian(galaxy *gals, int numGals){
	// Declare the cartesian coordinate galaxy list.
	cartesianGalaxy *cartGals = malloc(numGals * sizeof(cartesianGalaxy));

	// Setup the interpolator.
	double *rs, *zs;
	int numInterpPts;
	setupRedshiftInterp(&rs, &zs, &numInterpPts);

	// Loop through each galaxy in parallel and convert to cartesian
	// coordinates.
	int i;
	//omp_set_num_threads(32);
	#pragma omp parallel for
	for(i = 0; i < numGals; i++){
		// Convert ra and dec.
		double phi = (gals+i)->ra * (PI/180);
		double theta = PI/2.0 - (gals+i)->dec * (PI/180);

		// Calculate the line of sight distance.
		double r = interpDist((gals+i)->z_red, rs, zs, numInterpPts);

		// Set the x, y, and z coordinates.
		(cartGals+i)->x = r * sin(theta) * cos(phi);
		(cartGals+i)->y = r * sin(theta) * sin(phi);
		(cartGals+i)->z = r * cos(theta);
	}

	return cartGals;
}

galaxy *convertFromCartesian(cartesianGalaxy *cartGals, galaxy *origGals,
	int numGals){

	// Declare the regular coordinate galaxy list.
	galaxy *gals = malloc(numGals * sizeof(galaxy));

	// Setup the interpolator.
	double *rs, *zs;
	int numInterpPts;
	setupRedshiftInterp(&rs, &zs, &numInterpPts);

	// Loop through each galaxy and perform root finding to calculate the
	// redshift, using bisection.
	int i;
	#pragma omp parallel for
	for(i = 0; i < numGals; i++){
		double x = (cartGals+i)->x;
		double y = (cartGals+i)->y;
		double z = (cartGals+i)->z;

		// Calculate the radial distance.
		double r = pow(x*x+y*y+z*z,0.5);

		// Set the values for ra, dec, and z.
		(gals+i)->ra = (origGals+i)->ra;
		(gals+i)->dec = (origGals+i)->dec;
		(gals+i)->z_red = interpRedshift(r, rs, zs, numInterpPts);
		(gals+i)->z_err = (origGals+i)->z_err;
	}

	return gals;
}

galaxy *trimGalaxyList(galaxy *gals, int *numGals, int xStart, int yStart,
	int zStart, int boxLength){

	// Convert the galaxy list to cartesian coordinates.
	cartesianGalaxy *cartGals = convertToCartesian(gals, *numGals);

	// Loop through each galaxy and check that it is within the outer bounds.
	// If so mark it as 1, otherwise 0.
	int *index = malloc(*numGals * sizeof(int));
	int i;
	int sum = 0;
	for(i = 0; i < *numGals; i++){
		if(((cartGals+i)->x >= xStart &&
				(cartGals+i)->x <= (xStart + boxLength)) &&
			((cartGals+i)->y >= yStart &&
				(cartGals+i)->y <= (yStart + boxLength)) &&
			((cartGals+i)->z >= zStart &&
				(cartGals+i)->z <= (zStart + boxLength))){
			
			*(index+i) = 1;
			sum++;
		}else{
			*(index+i) = 0;
		}
	}

	// Create a new list of galaxies from the trimmed subset.
	galaxy *trimmedGals = malloc(sum * sizeof(galaxy));
	galaxy *curGal = trimmedGals;
	for(i = 0; i < *numGals; i++){
		if(*(index+i)){
			curGal->ra = (gals+i)->ra;
			curGal->dec = (gals+i)->dec;
			curGal->z_red = (gals+i)->z_red;
			curGal->z_err = (gals+i)->z_err;
			curGal++;
		}
	}

	// Free up the old cartesian galaxy list.
	free(cartGals);

	// Update the number of galaxies.
	*numGals = sum;

	return trimmedGals;
}

/*
galaxy *trimGalaxyList(galaxy *gals, int *numGals){
	// Read in the specifications for trimming the galaxies.
	FILE *fp;
	fp = fopen("bounds.conf", "r");
	int xOutMin, xInMin, xInMax, xOutMax;
	int yOutMin, yInMin, yInMax, yOutMax;
	int zOutMin, zInMin, zInMax, zOutMax;
	int voxelLength;
	fscanf(fp, "X = %d %d %d %d\n", &xOutMin, &xInMin, &xInMax, &xOutMax);
	fscanf(fp, "Y = %d %d %d %d\n", &yOutMin, &yInMin, &yInMax, &yOutMax);
	fscanf(fp, "Z = %d %d %d %d\n", &zOutMin, &zInMin, &zInMax, &zOutMax);
	fscanf(fp, "Voxel Length = %d", &voxelLength);
	fclose(fp);

	// Convert the galaxy list to cartesian coordinates.
	cartesianGalaxy *cartGals = convertToCartesian(gals, *numGals);

	// Loop through each galaxy and check that it is within the outer bounds.
	// If so mark it as 1, otherwise 0.
	int *index = malloc(*numGals * sizeof(int));
	int i;
	int sum = 0;
	for(i = 0; i < *numGals; i++){
		if(((cartGals+i)->x >= xOutMin && (cartGals+i)->x <= xOutMax) &&
			((cartGals+i)->y >= yOutMin && (cartGals+i)->y <= yOutMax) &&
			((cartGals+i)->z >= zOutMin && (cartGals+i)->z <= zOutMax)){
			
			*(index+i) = 1;
			sum++;
		}else{
			*(index+i) = 0;
		}
	}

	// Create a new list of galaxies from the trimmed subset.
	galaxy *trimmedGals = malloc(sum * sizeof(galaxy));
	galaxy *curGal = trimmedGals;
	for(i = 0; i < *numGals; i++){
		if(*(index+i)){
			curGal->ra = (gals+i)->ra;
			curGal->dec = (gals+i)->dec;
			curGal->z_red = (gals+i)->z_red;
			curGal->z_err = (gals+i)->z_err;
			curGal++;
		}
	}

	// Free up the old cartesian galaxy list.
	free(cartGals);

	// Update the number of galaxies.
	*numGals = sum;

	return trimmedGals;
}
*/