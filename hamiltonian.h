#include <mkl_lapacke.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spline.h>
#include "map_likelihood.h"

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

double *invertDiagMat(double *mat, int numVoxelsPerDim);
double *generateMassMatDiag(double *invCov, int *voxels, int numVoxelsPerDim);
double *generateMomenta(int numVoxelsPerDim);
double *modifyMap(double *cov, double *invCov, double *map, int *voxels,
	int numVoxelsPerDim);
double springForce(double *x, int ind, void*);
double potentialForce(double *x, int ind, void *params);
void leapfrogIntegrator(double *x, double *p, double *M, int numBodies,
	int numSteps, double epsilon, double (*force)(double*,int,void*), void *params);
typedef struct{
	double *k;
} hoParams;
typedef struct{
	double mean;
	double *invCov;
	int *voxels;
	int numVoxels;
	double avgN;
}potentialParams;

#endif