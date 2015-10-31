#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_spline.h>

typedef struct{
	double x;
	double y;
	double z;
} galaxy;

galaxy *importGalaxies(char *filename);
gsl_matrix *generateCov(int numVoxelsPerDim,
	double voxelLength);

void main(){
	char filename[100] = "sample.txt";

	printf("Here-1\n");

	// Load the galaxy data.
	galaxy *gals = importGalaxies(filename);

	printf("Here0\n");

	// Declare and intialize various values.
	int boxLength = 500; // Mpc
	int numVoxelsPerDim = 16;
	double voxelLength = boxLength/numVoxelsPerDim;

	printf("Here1\n");

	// Generate the covariance matrix.
	generateCov(numVoxelsPerDim, voxelLength);

	printf("Here2\n");
}

gsl_matrix *generateCov(int numVoxelsPerDim,
	double voxelLength){

	int n = pow(numVoxelsPerDim,3);

	// Declare and initialize the distances matrix.
	gsl_matrix *dists = gsl_matrix_calloc(n,n);

	// Loop through each pair of two voxels and calcualte
	// the distances between them.
	int x,y,z;
	int i,j,k;
	for(x = 0; x < numVoxelsPerDim; x++){
		for(y = 0; y < numVoxelsPerDim; y++){
			for(z = 0; z < numVoxelsPerDim; z++){
				int a = x + y*numVoxelsPerDim +
					z*pow(numVoxelsPerDim,2);

				for(i=0; i < numVoxelsPerDim; i++){
					for(j=0; j < numVoxelsPerDim; j++){
						for(k=0; k < numVoxelsPerDim; k++){
							int b = i + j*numVoxelsPerDim +
								k*pow(numVoxelsPerDim,2);

							// Calculate the distance
							// between the two voxels.
							double r = voxelLength *
								sqrt(pow(x-i,2) + 
								pow(y-j,2) + pow(k-z,2));

							// Set the values of the matrix.
							gsl_matrix_set(dists,a,b,r);
						}
					}
				}
			}
		}
	}

	printf("Here1.5\n");

	// Declare and initialize the xi samples.
	int numSamps = 5;
	double rSamp[] = {0, 250, 500, 750, 1000};
	double xiSamp[] = {1, 0.7, 0.4, 0.1, 0.01};

	// Initialize the spline interpolator.
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline =
    	gsl_spline_alloc(gsl_interp_cspline, numSamps);

    // Fit the spline.
    gsl_spline_init(spline, rSamp, xiSamp, numSamps);

    // Initialize the covariance matrix.
    gsl_matrix *cov = gsl_matrix_calloc(n,n);

    // Loop through each value in the distance matrix
    // and calculate the covariance matrix.
    double r, Qij;
    for(i = 0; i < n; i++){
    	for(j = 0; j < n; j++){
    		// Get the value of r at the current location.
    		r = gsl_matrix_get(dists, i, j);

    		// Calculate the value of the covariance
    		// matrix at the current location.
    		Qij = log(1 + gsl_spline_eval(spline,
    			r, acc));

    		// Set the value of the covariance matrix.
    		gsl_matrix_set(cov, i, j, Qij);
    	}
    }

    printf("Here1.75\n");

    // Perform cholesky decomposition then invert.
    gsl_matrix *inv_cov = cov;
    gsl_linalg_cholesky_decomp(cov);
    gsl_linalg_cholesky_invert(cov);

	FILE *fp;
	fp = fopen("out.txt", "w");
	print_matrix(fp, inv_cov);
}

galaxy *importGalaxies(char *filename){
	// Open the file.
	FILE *fp;
	fp = fopen(filename, "r");

	// Get the number of lines in the file.
	char c;
	int numGals = 0;
	while((c = getc(fp)) != EOF){
		if(c == '\n'){
			numGals++;
		}
	}
	rewind(fp);

	// Allocate the catalog of galaxies.
	galaxy *gals = malloc(numGals*sizeof(galaxy));

	// Load each galaxy.
	int i;
	galaxy *curGal = gals;
	for(i = 0; i < numGals; i++){
		fscanf(fp, "%le %le %le", &curGal->x, 
			&curGal->y, &curGal->z);

		curGal++;
	}

	return gals;
}

int print_matrix(FILE *f, const gsl_matrix *m){
        int status, n = 0;

        size_t i,j;
        for (i = 0; i < m->size1; i++) {
                for (j = 0; j < m->size2; j++) {
                        if ((status = fprintf(f, "%e ",
                        	gsl_matrix_get(m, i, j))) < 0)
                                return -1;
                        n += status;
                }

                if ((status = fprintf(f, "\n")) < 0)
                        return -1;
                n += status;
        }

        return n;
}