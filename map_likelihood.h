#include <math.h>
#include <gsl/gsl_spline.h>
#include <mkl_lapacke.h>

#ifndef MAP_LIKELIHOOD_H
#define MAP_LIKELIHOOD_H

double *generateCov(int numVoxelsPerDim, int boxLength,
	gsl_spline *spline);
gsl_spline *initCorrSpline(int numSamps, double rSamp[], double xiSamp[]);
double *invertCov(double *cov, int numVoxelsPerDim);

#endif