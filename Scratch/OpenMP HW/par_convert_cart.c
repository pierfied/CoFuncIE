/*

Reads in a set of ra, dec, and z coordinates and converts
them into cartesian coordinates. The format for the input
is ra, dec, and z in that order, separated by spaces.
Each galaxy coordinate data is on its own line.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

typedef struct {
	double ra;
	double dec;
	double z;
} galaxy;

typedef struct {
	double x;
	double y;
	double z;
} cartesianGalaxy;

cartesianGalaxy convertCartesian(galaxy *gal);
double comovingDistance(double z);

int main(){
	// Count the number of lines.
	int numLines = 0;
	char tmp[100];
	while(fgets(tmp, sizeof(tmp), stdin) != NULL){
		numLines++;
	}
	int numGals = numLines;
	rewind(stdin);

	// Create the catalog.
	galaxy *catalog = malloc(numGals*sizeof(galaxy));

	// Read in each variable.
	int i;
	galaxy *curGal = catalog;
	for(i = 0; i < numLines; i++){
		fscanf(stdin, "%le %le %le", &curGal->ra, 
			&curGal->dec, &curGal->z);
		curGal++;
	}

	// Create the cartesian catalog.
	cartesianGalaxy *cartCatalog =
		malloc(numGals*sizeof(cartesianGalaxy));

	// Convert each galaxy position into cartesian
	// coordinates.
	curGal = catalog;
	cartesianGalaxy *curCartGal = cartCatalog;
	i = 0;
	omp_set_num_threads(4);
	while(i < numGals){
		#pragma omp parallel
		{
			galaxy *curThreadGal;
			cartesianGalaxy *curThreadCartGal;

			#pragma omp critical
			{
				curThreadGal = curGal++;
				curThreadCartGal = curCartGal++;
				i++;
			}

			*curThreadCartGal =
				convertCartesian(curThreadGal);
		}
	}

	// Print out the results to a txt file.
	FILE *file;
	file = fopen("cartesian_coords.txt","w+");
	curCartGal = cartCatalog;
	for(i = 0; i < numGals; i++){
		fprintf(file, "%e %e %e\n", curCartGal->x,
			curCartGal->y, curCartGal->z);
		curCartGal++;
	}
	close(file);
}

cartesianGalaxy convertCartesian(galaxy *gal){
	// Convert ra and dec.
	double phi = gal->ra * (PI/180);
	double theta = PI/2.0 - gal->dec * (PI/180);

	// Calculate the line of sight distance.
	double r = comovingDistance(gal->z);

	// Calculate the x, y, and z coordiantes.
	cartesianGalaxy cartGal;
	cartGal.x = r * sin(theta) * cos(phi);
	cartGal.y = r * sin(theta) * sin(phi);
	cartGal.z = r * cos(theta);

	return cartGal;
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