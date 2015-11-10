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
		fscanf(fp, "%le, %le, %le", &(curGal->ra), &(curGal->dec),
			&(curGal->z_red));
		curGal++;
	}
	fclose(fp);

	// Return the list of galaxies.
	return gals;
}

double comovingDistance(double z){
	// Define various cosmological constants.
	double M = cosmo_M;
	double k = cosmo_k;
	double lambda = cosmo_lambda;
	double h = cosmo_h;
	double Dh = 3000/h; // [Mpc]

	// Perform the integration.
	int numSteps = 5000;
	double stepsize = z/numSteps;
	double curZ;
	double fx, fxdx;
	double sum = 0;
	for(curZ = 0; curZ < z; curZ += stepsize){
	  fx = 1./sqrt(M*pow(1+curZ,3)
		       + k*pow(1+curZ,2) + lambda);
	  fxdx = 1./sqrt(M*pow(1+curZ+stepsize,3) 
			 + k*pow(1+curZ+stepsize,2) + lambda);
	  sum += stepsize/2.*(fx+fxdx);
	}

	// Error is proportional to this factor
	// We can choose to print this out if we want
	double error = -Dh*z*z*z/(12*numSteps*numSteps);
	//printf("Error in comoving distance integration = %f\n",error);

	return Dh * sum;
}

void setupRedshiftInterp(double **r, double **z, int *numInterpPts){
	// Define various parameters.
	*numInterpPts = 101;
	*r = malloc(*numInterpPts * sizeof(double));
	*z = malloc(*numInterpPts * sizeof(double));
	double zMin = 0;
	double zMax = 1;
	double stepsize = (zMax - zMin)/(*numInterpPts - 1);

	// Set the values at r = 0.
	**r = 0;
	**z = 0;

	// Loop through each interpolation point and calculate the value of z.
	int i;
	for(i = 1; i < *numInterpPts; i++){
		*(*z+i) = zMin + i * stepsize;
		*(*r+i) = comovingDistance(*(*z+i));
	}
}

double interpRedshift(double rQuery, double *r, double *z, int numInterpPts){
	// Perform linear interpolation for the redshift at the query radius.
	int i;
	for(i = 0; i < (numInterpPts - 1); i++){
		if(*(r+i+1) > rQuery){
			double zInterp = *(z+i)
				+ (*(z+i+1) - *(z+i))*(rQuery - *(r+i))/(*(r+i+1) - *(r+i));
			return zInterp;
		}
	}

	return -1;
}

double interpDist(double zQuery, double *r, double *z, int numInterpPts){
	// Perform linear interpolation for the redshift at the query radius.
	int i;
	for(i = 0; i < (numInterpPts - 1); i++){
		if(*(z+i+1) > zQuery){
			double zInterp = *(r+i)
				+ (*(r+i+1) - *(r+i))*(zQuery - *(z+i))/(*(z+i+1) - *(z+i));
			return zInterp;
		}
	}

	return -1;
}

cartesianGalaxy *convertToCartesian(galaxy *gals, int numGals){
	// Declare the cartesian coordinate galaxy list.
	cartesianGalaxy *cartGals = malloc(numGals * sizeof(cartesianGalaxy));

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
		double r = comovingDistance((gals+i)->z_red);

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

		double threshold = 1e-6;
		double a = 0.01;
		double b = 1;
		double z_red, mid, fa, fmid;

		// Perform the bisection root finding.
		while(1){
			mid = (a + b)/2.0;
			fmid = comovingDistance(mid) - r;

			// Check for convergence.
			if(fabs(fmid) < threshold || (b-a)/2.0 < threshold){
				z_red = mid;
				break;
			}

			// Determine the new interval.
			fa = comovingDistance(a) - r;
			if((fa > 0 && fmid > 0) || (fa < 0 && fmid < 0)){
				a = mid;
			}else{
				b = mid;
			}
		}

		// Set the values for ra, dec, and z.
		(gals+i)->ra = (origGals+i)->ra;
		(gals+i)->dec = (origGals+i)->dec;
		(gals+i)->z_red = z_red;
	}

	return gals;
}

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
			curGal++;
		}
	}

	// Free up the old cartesian galaxy list.
	free(cartGals);

	// Update the number of galaxies.
	*numGals = sum;

	return trimmedGals;
}
