#include "hamiltonian.h"

double *generateMassMatDiag(double *invCov, int *voxels, int numVoxelsPerDim){
	// Allocate space for the mass matrix.
	int n = pow(numVoxelsPerDim, 3);
	double *M = malloc(n * sizeof(double));

	// Loop through each 
	int i;
	for(i = 0; i < n; i++){
		// Calculate the average galaxy count.
		double avgN = 0;
		int j;
		for(j = 0; j < n; j++){
			avgN += voxels[i];
		}
		avgN /= n;

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