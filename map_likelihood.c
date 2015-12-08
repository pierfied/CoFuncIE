#include "map_likelihood.h"

#define PI 3.14159265358979323846

double mapLnLikelihood(double *map, int *voxels, int numVoxelsPerDim,
	int boxLength, gsl_spline *spline){

	// Initialize the covariance matrix.
	int n = pow(numVoxelsPerDim, 3);
	double *cov = generateCov(numVoxelsPerDim, boxLength, spline);
	double *invCov = invertCov(cov, numVoxelsPerDim);

	// Calculate the mean value for the lognormal map.
	double mean = -0.5*(*cov);

	// Calculate the first term.
	double firstTerm = 0;
	int i,j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// Calculate the index.
			int index = i + n*j;

			double yi = log(1 + *(map+i));
			double yj = log(1 + *(map+j));

			double invQij = *(invCov + index);

			// Add to the sum.
			firstTerm += (yi - mean) * invQij * (yj - mean);
		}
	}
	firstTerm *= -0.5;

	// Calculate the second term.
	double secondTerm = 0;
	for(i = 0; i < n; i++){
		secondTerm += log(1/(1+*(map+i)));
	}

	// Calculate the mean galaxy count.
	double avgN = 0;
	for(i = 0; i < n; i++){
		avgN += *(voxels + i);
	}
	avgN /= n;

	// Calculate the third term.
	double thirdTerm = 0;
	for(i = 0; i < n; i++){
		double R = 1;
		int Nk = *(voxels + i);
		double lambdak = avgN * (1 + *(map+i));

		thirdTerm += Nk*log(lambdak) - lambdak - lnfactorial(Nk);
	}

	return firstTerm + secondTerm + thirdTerm;
}

double lnfactorial(int x){
	if(x <= 1){
		return 0;
	}else{
		double lnfac = 0;
		while(x > 1){
			lnfac += log(x--);
		}

		return lnfac;
	}
}

double *invertCov(double *cov, int numVoxelsPerDim){
	int n = pow(numVoxelsPerDim, 3);

	// Create the inverse covariance matrix and copy the values into it.
	double *invCov = malloc(n*n * sizeof(double));
	int i,j;
	#pragma omp parallel for private(j)
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
	#pragma omp parallel for private(j)
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// Check that the value being looked at is within the upper
			// triangular matrix.
			if(i <= j) continue;

			// Copy the value from the upper to the lower triangle.
			int a = i + n*j;
			int b = n*i + j;
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

	// Loop through each pair of two voxels and calcualte
	// the distances between them.
	int x,y,z;
	int i,j,k;
	#pragma omp parallel for private(y,z,i,j,k)
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

	// Initialize the spline accelerator.
	gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Initialize the covariance matrix.
    double *cov = calloc(n*n, sizeof(double));

    // Loop through each value in the distance matrix
    // and calculate the covariance matrix.
    #pragma omp parallel for private(j)
    for(i = 0; i < n; i++){
    	for(j = 0; j < n; j++){
    		// Get the value of r at the current location.
    		int index = i + n*j;
    		double r = *(dists + index);

    		// Calculate and set the value of the covariance
    		// matrix at the current location.
    		//*(cov + index) = log(1 + gsl_spline_eval(spline, r, acc));
    		if(i == j){
    			*(cov + index) = log(2);
    		}else{
    			*(cov + index) = log(1 + pow(r/(10./0.7),-2));
    		}
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
