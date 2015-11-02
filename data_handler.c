#include <stdio.h>
#include <math.h>

#define PI 3.14159265359

typedef struct{
	double ra;
	double dec;
	double z_red;
} galaxy;

typedef struct {
	double x;
	double y;
	double z;
} cartesianGalaxy;

int getNumGals(){
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
	int numGals = 0;
	while(fgets(loc, sizeof(loc), fp) != NULL){
		numGals++;
	}
	fclose(fp);

	return numGals;
}

galaxy *readData(){
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
	int numGals = 0;
	while(fgets(loc, sizeof(loc), fp) != NULL){
		numGals++;
	}
	rewind(fp);

	// Declare the galaxy list.
	galaxy *gals = malloc(numGals * sizeof(galaxy));

	// Read in the galaxies.
	int i;
	galaxy *curGal = gals;
	for(i = 0; i < numGals; i++){
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
	double M = 0.3;
	double k = 0;
	double lambda = 0.7;
	double h = 0.7;
	double Dh = 3000/h; // [Mpc]

	// Perform the integration.
	int numSteps = 1000;
	double stepsize = z/numSteps;
	double curZ;
	double sum = 0;
	for(curZ = 0; curZ < z; curZ += stepsize){
		sum += stepsize / sqrt(M*pow(1+curZ,3)
			+ k*pow(1+curZ,2) + lambda);
	}

	return Dh * sum;
}

cartesianGalaxy *convertCartesian(galaxy *gals, int numGals){
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
