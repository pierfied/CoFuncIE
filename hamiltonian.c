#include "hamiltonian.h"

double *modifyMap(double *cov, double *invCov, double *map, int *voxels,
	int numVoxelsPerDim){

	// Generate random momenta and the mass matrix and the inverse.
	double *p, *M, *invM;
	p = generateMomenta(numVoxelsPerDim);
	M = generateMassMatDiag(invCov, voxels, numVoxelsPerDim);
	invM = invertDiagMat(M, numVoxelsPerDim);

	// Generate a random number of steps to take and a random stepsize
	// between 0 and 2.
	srand((unsigned)time(NULL));
	int n = 1e4;//*((double)rand()/(double)RAND_MAX);
	double epsilon = 1e-2;//*((double)rand()/(double)RAND_MAX);

	printf("n = %d\ne = %f\n", n, epsilon);

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
	double sum = 0;
	for(i = 0; i < numVoxels; i++){
		newMap[i] = exp(r[i]) - 1;
		sum += newMap[i];
		if(fabs(r[i]) > 10 || r[i]!=r[i]){
			printf("Warning i = %d\tr = %f\tp = %f\tmap = %d\n", i, r[i], p[i], voxels[i]);
		}
	}
	//printf("Average: %f\n", sum/m);

	// Free up unused variables.
	free(r);
	//free(grad);

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
	double mean, avgN, *invCov;
	int *voxels, numVoxels;

	mean = ((potentialParams*)params)->mean;
	invCov = ((potentialParams*)params)->invCov;
	voxels = ((potentialParams*)params)->voxels;
	numVoxels = ((potentialParams*)params)->numVoxels;
	avgN = ((potentialParams*)params)->avgN;

	double firstTerm = 0;
	int i;
	for(i = 0; i < numVoxels; i++){
		int index = ind + numVoxels*i;

		firstTerm += invCov[index]*(x[ind] - mean);
	}
	firstTerm *= -1;

	double secondTerm = voxels[ind] - avgN * exp(x[ind]);

	//*****************************************************
	if(ind == 0){
		double *map = malloc(sizeof(double)*numVoxels);
		for(i = 0; i < numVoxels; i++){
			map[i] = exp(x[i]) - 1;
		}

		int numSamps = 5;
		double rSamp[] = {0, 250, 500, 750, 1000};
		double xiSamp[] = {0.1, 0.25, 0.5, 0.75, 1};

		gsl_spline *spline = initCorrSpline(numSamps, rSamp, xiSamp);

		double lnLikeMap = mapLnLikelihood(map, voxels, 10, 400,
			spline);

		FILE *fp;
		fp = fopen("hamilton2_long.csv","a");

		fprintf(fp, "%f, ", lnLikeMap);

		free(spline);
		free(map);
		fclose(fp);
	}
	if(ind == -1){
		FILE *fp;
		fp = fopen("likelihood.txt","a");

		double *map = malloc(sizeof(double)*numVoxels);
		for(i = 0; i < numVoxels; i++){
			map[i] = exp(x[i]) - 1;
		}

		int numSamps = 5;
		double rSamp[] = {0, 250, 500, 750, 1000};
		double xiSamp[] = {0.1, 0.25, 0.5, 0.75, 1};

		gsl_spline *spline = initCorrSpline(numSamps, rSamp, xiSamp);

		double lnLikeMap = mapLnLikelihood(map, voxels, 10, 400,
			spline);

		fprintf(fp, "%f\n", lnLikeMap);

		fclose(fp);
	}

	return (firstTerm + secondTerm);
}

void leapfrogIntegrator(double *x, double *p, double *M, int numBodies,
	int numSteps, double epsilon, double (*force)(double*,int,void*),
	void *params){

	// Perform the first half-step in momenta.
	int i;
	for(i = 0; i < numBodies; i++){
		p[i] += 0.5 * epsilon * (*force)(x, i, params);
	}
	printf("\n");

	FILE *fp;
	fp = fopen("hamilton2_long.csv", "w");
	fclose(fp);
	
	// Loop through each step.
	for(i = 0; i < numSteps; i++){
		printf("%d\n",i);
		// Perform the full-step in position for all bodies.
		int j;
		for(j = 0; j < numBodies; j++){
			x[j] += epsilon * p[j] / M[j];
		}

		double *halfP = malloc(sizeof(double) * numBodies);

		// Perform the full-step in momentum for all bodies.
		for(j = 0; j < numBodies; j++){
			double tmp = (*force)(x,j,params);

			// Calculate the k
			halfP[j] = p[j] + epsilon*0.5*tmp;

			p[j] += epsilon * tmp;
			//p[j] += epsilon * (*force)(x, j, params);
		}

		// Calculate the kinetic energy.
		fp = fopen("hamilton2_long.csv","a");
		double kineticTerm = 0;
		for(j = 0; j < numBodies; j++){
			kineticTerm += pow(halfP[j],2)/M[j];
		}
		kineticTerm /= 2;
		fprintf(fp,"%f, %f, %f\n", kineticTerm,x[0],p[0]);
		fclose(fp);

		free(halfP);

		if(i % 100 == -1){
			for(j = 0; j < numBodies; j++){
				p[j] = 0;
			}
		}

		/*
		// Print the position and momentum at the current timestep for all M<0.
		fprintf(fp, "%f", i*epsilon);
		for(j = 0; j < numBodies; j++){
			if(1 > 0){
				fprintf(fp, ",%f,%f,%f", M[j], x[j], p[j]);
			}
		}
		fprintf(fp, "\n");
		*/
	}
}