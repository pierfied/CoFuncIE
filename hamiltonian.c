#include "hamiltonian.h"

double *modifyMap(double *cov, double *invCov, double *map, int *voxels,
	int numVoxelsPerDim){

	// Generate random momenta and the mass matrix and the inverse.
	double *p, *M, *invM;
	p = generateMomenta(numVoxelsPerDim);
	M = generateMassMatDiag(invCov, voxels, numVoxelsPerDim);
	invM = invertDiagMat(M, numVoxelsPerDim);

	// Generate a random number of steps so that the total is of order
	// one for the given dt.
	srand((unsigned)time(NULL));
	int n = 1e2*(0.8 + 0.4*((double)rand()/(double)RAND_MAX));
	double epsilon = 1e-2;

	// Create a new map, which will be called r (position vector).
	// Copy values of the original map into this one.
	int i;
	int numVoxels = pow(numVoxelsPerDim, 3);
	double *r = malloc(numVoxels * sizeof(double));

	// Create the position values from the denstiy maps.
	for(i = 0; i < numVoxels; i++){
		r[i] = log(1 + map[i]);
	}

	// Calculate the average galaxy count.
	double avgN = 0;
	int j;
	for(j = 0; j < numVoxels; j++){
		avgN += voxels[j];
	}
	avgN /= numVoxels;

	// Set the mass of each voxel equal to the average galaxy count
	// of the entire survey.
	for(i = 0; i < numVoxels; i++){
		p[i] = 0;
		M[i] = avgN;
	}

	// Calculate the mean value of r.
	double mean = -0.5*(*cov);

	// Create the parameter structure.
	potentialParams params;
	params.mean = mean;
	params.invCov = invCov;
	params.voxels = voxels;
	params.numVoxels = numVoxels;
	params.avgN = avgN;

	// Create the function pointer for the hamiltonian.
	double (*force)(double*,int,void*);
	force = &potentialForce;

	// Perform the integration.
	leapfrogIntegrator(r, p, M, numVoxels, n, epsilon, force, &params);

	// Create the new map.
	double *newMap = malloc(numVoxels * sizeof(double));
	for(i = 0; i < numVoxels; i++){
		newMap[i] = exp(r[i]) - 1;
	}

	// Free up unused variables.
	free(r);

	return newMap;
}

double *generateMassMatDiag(double *invCov, int *voxels, int numVoxelsPerDim){
	// Allocate space for the mass matrix.
	int n = pow(numVoxelsPerDim, 3);
	double *M = malloc(n * sizeof(double));

	// Calculate the average galaxy count.
	double avgN = 0;
	int i;
	for(i = 0; i < n; i++){
		avgN += voxels[i];
	}
	avgN /= n;

	// Loop through each diagonal mass entry.
	for(i = 0; i < n; i++){
		// Set the value of the mass matrix.
		M[i] = invCov[i] - ((voxels[i]-avgN) - avgN);
	}

	return M;
}

double *invertDiagMat(double *mat, int numVoxelsPerDim){
	int n = pow(numVoxelsPerDim, 3);

	// Create the inverse mass matrix and copy the values into it.
	double *invMat = malloc(n * sizeof(double));
	int i;
	for(i = 0; i < n; i++){
		invMat[i] = 1/mat[i];
	}

	return invMat;
}

double *generateMomenta(int numVoxelsPerDim){
	// Initialize the random number generator.
	const gsl_rng_type *T;
  	gsl_rng *r;
  	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 1);//time(NULL)); // Set an initial seed of 0.

	// Allocate the space for the momenta.
	int n = pow(numVoxelsPerDim, 3);
	double *p = malloc(n * sizeof(double));

	// Initialize the momentum values with normally distributed values with
	// a mean of zero and unit variance.
	double sigma = 1;
	int i;
	for(i = 0; i < n; i++){
		p[i] = gsl_ran_gaussian(r, sigma);
	}

	gsl_rng_free(r);

	return p;
}

double springForce(double *x, int ind, void *params){
	double k = ((hoParams*)params)->k[ind];
	return -k*x[ind];
}

double potentialForce(double *x, int ind, void *params){
	// Declare necessary variables.
	double mean, avgN, *invCov;
	int *voxels, numVoxels;

	// Set the variables to the values held in the parameter struct.
	mean = ((potentialParams*)params)->mean;
	invCov = ((potentialParams*)params)->invCov;
	voxels = ((potentialParams*)params)->voxels;
	numVoxels = ((potentialParams*)params)->numVoxels;
	avgN = ((potentialParams*)params)->avgN;

	int i;
	// Calculate the first term of the force.
	double force = 0;
	for(i = 0; i < numVoxels; i++){
		// Get the index for the invcov.
		int index = ind + numVoxels*i;

		force += -invCov[index]*(x[i] - mean);
	}

	// Calculate the second term of the force.
	force += voxels[ind] - avgN*exp(x[ind]);

	return force;
}

void leapfrogIntegrator(double *x, double *p, double *M, int numBodies,
	int numSteps, double epsilon, double (*force)(double*,int,void*),
	void *params){

	// Perform the first half-step in momenta.
	int i;
	#pragma omp parallel for
	for(i = 0; i < numBodies; i++){
		p[i] += 0.5 * epsilon * (*force)(x, i, params);
	}
	
	// Loop through each step.
	for(i = 0; i < numSteps; i++){
		// Perform the full-step in position for all bodies.
		int j;
		#pragma omp parallel for
		for(j = 0; j < numBodies; j++){
			x[j] += epsilon * p[j] / M[j];
		}

		// Perform the full-step in momentum for all bodies.
		#pragma omp parallel for
		for(j = 0; j < numBodies; j++){
			p[j] += epsilon * (*force)(x, j, params);
		}
	}
}
