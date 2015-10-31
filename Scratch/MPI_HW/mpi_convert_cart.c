#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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

void runMaster(char *fname, int world_size);
void runSlave(int rank);
double comovingDistance(double z);

void main(int argc, char *argv[]){
	// Initialize MPI.
	MPI_Init(NULL, NULL);

	// Get the rank and size.
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Check if this is the master or a slave node.
	if(rank == 0){
		// Check that there is a filename passed in
		// as input.
		if(argc < 2){
			printf("No filename was passed as input.\n");
			exit(-1);
		}

		runMaster(argv[1], world_size);
	}else{
		runSlave(rank);
	}

	// Close MPI
	MPI_Finalize();
}

void runMaster(char *fname, int world_size){
	// Open the file.
	FILE *file;
	file = fopen(fname, "r");

	// Count the number of lines.
	int numLines = 0;
	char tmp[100];
	while(fgets(tmp, sizeof(tmp), file) != NULL){
		numLines++;
	}
	int numGals = numLines;
	rewind(file);

	printf("Total Gals: %d\n", numGals);

	// Create the ra, dec, and z vars.
	double *ra = malloc(sizeof(double)*numGals);
	double *dec = malloc(sizeof(double)*numGals);
	double *zRed = malloc(sizeof(double)*numGals);

	// Read in the data.
	int i;
	for(i = 0; i < numGals; i++){
		double x;
		fscanf(file, "%le %le %le", (ra + i),
			(dec + i), (zRed + i));
	}
	close(file);

	// Split and send the data to each node.
	int numGalsPerNode =
		(int)ceil((double)numGals/(world_size-1));
	int startIndex = 0;
	for(i = 1; i < world_size; i++){
		int endIndex = startIndex + numGalsPerNode;

		// Check that the end index doesn't go above
		// the number of galaxies.
		if(endIndex >= numGals){
			endIndex = numGals;
		}

		// Calculate the number of galaxies to be
		// sent to the node.
		int numGalsToSend = endIndex - startIndex;

		// Tell the node how many galaxies will be sent.
		MPI_Send(&numGalsToSend, 1, MPI_INT, i, 0,
			MPI_COMM_WORLD);

		// Send the ra, dec, and z of each galaxy to
		// the node.
		MPI_Send((ra + startIndex), numGalsToSend,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		MPI_Send((dec + startIndex), numGalsToSend,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		MPI_Send((zRed + startIndex), numGalsToSend,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

		// Increment the index.
		startIndex = endIndex;
	}

	// Create the cartesian coordinate variables.
	double *x = malloc(sizeof(double)*numGals);
	double *y = malloc(sizeof(double)*numGals);
	double *z = malloc(sizeof(double)*numGals);

	// Receive the converted coordiantes from each
	// slave node.
	startIndex = 0;
	for(i = 1; i < world_size; i++){
		// Get the number of galaxies to be added to
		// the list.
		int numGalsToReceive;
		MPI_Recv(&numGalsToReceive, 1, MPI_INT, i, 0,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// Receive the cartesian coordinates for
		// each galaxy and store them.
		MPI_Recv((x + startIndex), numGalsToReceive,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv((y + startIndex), numGalsToReceive,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv((z + startIndex), numGalsToReceive,
			MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);

		// Update the startIndex.
		startIndex += numGalsToReceive;
	}

	// Print out the results to a txt file.
	file = fopen("cartesian_coords.txt","w+");
	for(i = 0; i < numGals; i++){
		fprintf(file, "%e %e %e\n", *(x + i),
			*(y + i), *(z + i));
	}
	close(file);
}

void runSlave(int rank){
	// Get the number of galxies being sent to this
	// node.
	int numGals;
	MPI_Recv(&numGals, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
		MPI_STATUS_IGNORE);

	printf("Num gals processor %d is receiving: %d\n",
		rank, numGals);

	// Receive ra, dec, and zRed.
	double *ra = malloc(sizeof(double)*numGals);
	double *dec = malloc(sizeof(double)*numGals);
	double *zRed = malloc(sizeof(double)*numGals);
	MPI_Recv(ra, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(dec, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(zRed, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// Create the cartesian variables.
	double *x = malloc(sizeof(double)*numGals);
	double *y = malloc(sizeof(double)*numGals);
	double *z = malloc(sizeof(double)*numGals);

	// Loop over all galaxies and convert to
	// cartesian coordinates.
	int i;
	for(i = 0; i < numGals; i++){
		// Convert ra and dec to phi and theta.
		double phi = *(ra + i) * (PI/180);
		double theta = PI/2.0 - *(dec + i) * (PI/180);

		// Calculate the line of sight distance.
		double r = comovingDistance(*(zRed + i));

		// Calculate the x, y, and z coordinates.
		*(x + i) = r * sin(theta) * cos(phi);
		*(y + i) = r * sin(theta) * sin(phi);
		*(z + i) = r * cos(theta);
	}

	// Remind the master of how many galaxies this
	// slave was given.
	MPI_Send(&numGals, 1, MPI_INT, 0, 0,
		MPI_COMM_WORLD);

	// Send the cartesian coordinates to the master.
	MPI_Send(x, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD);
	MPI_Send(y, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD);
	MPI_Send(z, numGals, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD);
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