/*
  This file contains functions used to calculate
  comoving distances and interpolate for a z(r)
  dependence.
*/

#include "cosmology.h"

double comovingDistance(double z){
	// Define various cosmological constants.
	double M = cosmo_M;
	double k = cosmo_k;
	double lambda = cosmo_lambda;
	double h = cosmo_h;
	double Dh = 3000/h; // [Mpc]

	// Perform the integration.
	int numSteps = 5000;
	double stepsize = z/numSteps;
	double curZ;
	double fx, fxdx;
	double sum = 0;
	for(curZ = 0; curZ < z; curZ += stepsize){
	  fx = 1./sqrt(M*pow(1+curZ,3)
		       + k*pow(1+curZ,2) + lambda);
	  fxdx = 1./sqrt(M*pow(1+curZ+stepsize,3) 
			 + k*pow(1+curZ+stepsize,2) + lambda);
	  sum += stepsize/2.*(fx+fxdx);
	}

	// Error is proportional to this factor
	// We can choose to print this out if we want
	double error = -Dh*z*z*z/(12*numSteps*numSteps);
	//printf("Error in comoving distance integration = %f\n",error);

	return Dh * sum;
}

void setupRedshiftInterp(double **r, double **z, int *numInterpPts){
	// Define various parameters.
	*numInterpPts = 101;
	*r = malloc(*numInterpPts * sizeof(double));
	*z = malloc(*numInterpPts * sizeof(double));
	double zMin = 0;
	double zMax = 1;
	double stepsize = (zMax - zMin)/(*numInterpPts - 1);

	// Set the values at r = 0.
	**r = 0;
	**z = 0;

	// Loop through each interpolation point and calculate the value of z.
	int i;
	for(i = 1; i < *numInterpPts; i++){
		*(*z+i) = zMin + i * stepsize;
		*(*r+i) = comovingDistance(*(*z+i));
	}
}

double interpRedshift(double rQuery, double *r, double *z, int numInterpPts){
	// Perform linear interpolation for the redshift at the query radius.
	int i;
	for(i = 0; i < (numInterpPts - 1); i++){
		if(*(r+i+1) > rQuery){
			double zInterp = *(z+i)
				+ (*(z+i+1) - *(z+i))*(rQuery - *(r+i))/(*(r+i+1) - *(r+i));
			return zInterp;
		}
	}

	return -1;
}

double interpDist(double zQuery, double *r, double *z, int numInterpPts){
	// Perform linear interpolation for the redshift at the query radius.
	int i;
	for(i = 0; i < (numInterpPts - 1); i++){
		if(*(z+i+1) > zQuery){
			double zInterp = *(r+i)
				+ (*(r+i+1) - *(r+i))*(zQuery - *(z+i))/(*(z+i+1) - *(z+i));
			return zInterp;
		}
	}

	return -1;
}
