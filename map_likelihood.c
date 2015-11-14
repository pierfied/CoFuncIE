#include "map_likelihood.h"



double *invertCov(double *cov, int numVoxelsPerDim){
	int n = pow(numVoxelsPerDim, 3);

	// Create the inverse covariance matrix and copy the values into it.
	double *invCov = malloc(n*n * sizeof(double));
	int i,j;
	#pragma omp parallel for
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			int index = i + n*j;
			*(invCov + index) = *(cov + index);
		}
	}

	// Perform the cholesky factorization.
	int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',n,invCov,n);

	// Perform the cholesky inversion.
	info = LAPACKE_dpotri(LAPACK_ROW_MAJOR,'U',n,invCov,n);

	// Copy the upper values to the lower diagonal.
	int x,y;
	#pragma omp parallel for
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			int a = n*i + j;
			int b = i + n*j;
			*(invCov + b) = *(invCov + a);
		}
	}

	return invCov;
}

double *generateCov(int numVoxelsPerDim, int boxLength,
	gsl_spline *spline){

	// Calculate the voxel length.
	double voxelLength = (double)boxLength/numVoxelsPerDim;

	// Allocate the distsances matrix.
	int n = pow(numVoxelsPerDim,3);
	double *dists = calloc(n*n, sizeof(double));

	// Declare and initialize the distances matrix.
	//gsl_matrix *dists = gsl_matrix_calloc(n,n);

	// Loop through each pair of two voxels and calcualte
	// the distances between them.
	int x,y,z;
	int i,j,k;
	#pragma omp parallel for
	for(x = 0; x < numVoxelsPerDim; x++){
		for(y = 0; y < numVoxelsPerDim; y++){
			for(z = 0; z < numVoxelsPerDim; z++){
				int a = x + y*numVoxelsPerDim +	z*pow(numVoxelsPerDim,2);

				for(i=0; i < numVoxelsPerDim; i++){
					for(j=0; j < numVoxelsPerDim; j++){
						for(k=0; k < numVoxelsPerDim; k++){
							int b = i + j*numVoxelsPerDim +
								k*pow(numVoxelsPerDim,2);

							// Calculate the distance
							// between the two voxels.
							double r = voxelLength *
								sqrt(pow(x-i,2) + pow(y-j,2) + pow(k-z,2));

							// Calculate the index in the distances matrix.
							int index = a + n*b;

							// Set the values of the matrix.
							*(dists + index) = r;
						}
					}
				}
			}
		}
	}

	// Declare and initialize the xi samples.
	int numSamps = 5;
	double rSamp[] = {0, 250, 500, 750, 1000};
	double xiSamp[] = {1, 0.7, 0.4, 0.1, 0.01};

	// Initialize the spline accelerator.
	gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Initialize the covariance matrix.
    double *cov = calloc(n*n, sizeof(double));

    // Loop through each value in the distance matrix
    // and calculate the covariance matrix.
    int index;
    double r, Qij;
    #pragma omp parallel for
    for(i = 0; i < n; i++){
    	for(j = 0; j < n; j++){
    		// Get the value of r at the current location.
    		index = i + n*j;
    		r = *(dists + index);

    		// Calculate and set the value of the covariance
    		// matrix at the current location.
    		*(cov + index) = log(1 + gsl_spline_eval(spline, r, acc));
    	}
    }

    return cov;
}

gsl_spline *initCorrSpline(int numSamps, double rSamp[], double xiSamp[]){
	// Initialize the spline interpolator.
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, numSamps);

    // Fit the spline.
    gsl_spline_init(spline, rSamp, xiSamp, numSamps);

    return spline;
}