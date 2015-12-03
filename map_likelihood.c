#include "map_likelihood.h"
#include <stdio.h>

#define PI 3.14159265358979323846

double mapLnLikelihood(double *map, int *voxels, int numVoxelsPerDim,
	int boxLength, gsl_spline *spline){

	// Initialize the covariance matrix.
	int n = pow(numVoxelsPerDim, 3);
	double *cov = generateCov(numVoxelsPerDim, boxLength, spline);
	double *invCov = invertCov(cov, numVoxelsPerDim);

	// Calculate the mean values for the lognormal map.
	double mean = 0;
	int i,j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// Calculate the index.
			int index = i + n*j;

			// Add to the sum.
			mean += *(cov + index);
		}

		mean *= 0.5;
	}

	printf("COV %e\n", *(cov+1111));
	printf("INV %e\n", *(invCov+1111));
	printf("mean: %f\n", mean);

	FILE *fp;
	fp = fopen("tmp2.out", "w");
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			int index = i + n*j;

			fprintf(fp, "%f\n", *(invCov + index));
		}
	}
	fclose(fp);

	// Calculate the first term.
	double firstTerm = 0;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			// Calculate the index.
			int index = i + n*j;

			double yi = log(1 + *(map+i));
			double yj = log(1 + *(map+j));

			if(i == j && yi != yj){
				printf("WARN\t%d\t%f\t%f\n",i , yi, yj);
				exit(0);
			}

			/*double ui = *(means + i);
			double uj = *(means + j);*/
			double invQij = *(invCov + index);

			if(index == n*n-1){
				printf("yi:\t%f\nyj:\t%f\n", yi, yj);
				printf("invQij\t%f\n", invQij);
			}

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
		double a = pow(avgN * (1 + *(map + i)), *(voxels + i));
		double b = -avgN * (1 + *(map + i));

		thirdTerm += log(a) + b - log(factorial(*(map + i)));
	}

	printf("f: %f\ns: %f\nt: %f\n", firstTerm, secondTerm, thirdTerm);

	return firstTerm + secondTerm + thirdTerm;

}

long factorial(long x){
	if(x <= 1){
		return 1;
	}else{
		return x * factorial(x-1);
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
	double *temp = malloc(n*n * sizeof(double));
	#pragma omp parallel for private(j)
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			int a = n*i + j;
			int b = i + n*j;
			*(temp + b) = *(invCov + a);
		}
	}

	return temp;
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
