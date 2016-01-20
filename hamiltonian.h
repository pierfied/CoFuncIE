#include <mkl_lapacke.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

double *invertDiagMat(double *mat, int numVoxelsPerDim);
double *generateMassMatDiag(double *invCov, int *voxels, int numVoxelsPerDim);
double *generateMomenta(int numVoxelsPerDim);
//double *modifyMap(double *map, int n);

#endif