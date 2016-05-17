#include "radial_likelihood.h"

double *fineMap(double *map, int numVoxelsPerDim, int resoultion){
	int newNumVPD = numVoxelsPerDim * resoultion;

	double *fMap = malloc(sizeof(double) * pow(newNumVPD,3));

	int i;
	for(i = 0; i < newNumVPD; i++){
		int j;
		for(j = 0; j < newNumVPD; j++){
			int k;
			for(k = 0; k < newNumVPD; k++){
				int index = i + j*newNumVPD + k*pow(newNumVPD,2);

				int minInterpBound = resoultion/2.0;
				int maxInterpBound = newNumVPD - resoultion/2.0;

				if(i > minInterpBound && j > minInterpBound
					&& k > minInterpBound && i < maxInterpBound
					&& j < maxInterpBound && k < maxInterpBound){

					int x = i - resoultion/2.0;
					int y = j - resoultion/2.0;
					int z = k - resoultion/2.0;

					fMap[index] = trilInterp(map, numVoxelsPerDim, resoultion,
						x, y, z);
				}else{
					int coarseX = i/resoultion;
					int coarseY = j/resoultion;
					int coarseZ = k/resoultion;

					int coarseInd = coarseX + coarseY*numVoxelsPerDim
						+ coarseZ*pow(numVoxelsPerDim,2);

					fMap[index] = map[coarseInd];
				}
			}
		}
	}

	return fMap;
}

double trilInterp(double *map, int numVoxelsPerDim, int resoultion,
	int x, int y, int z){

	int coarseX = x/resoultion;
	int coarseY = y/resoultion;
	int coarseZ = z/resoultion;

	double xd = x/(double)resoultion - coarseX;
	double yd = y/(double)resoultion - coarseY;
	double zd = z/(double)resoultion - coarseZ;

	//printf("%d %d %d   \t%f %f %f\n", x, y, z, xd, yd, zd);

	int ind000 = coarseX + coarseY*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind001 = coarseX + coarseY*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind010 = coarseX + (coarseY + 1)*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind011 = coarseX + (coarseY + 1)*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind100 = (coarseX + 1) + coarseY*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind101 = (coarseX + 1) + coarseY*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);
	int ind110 = (coarseX + 1) + (coarseY + 1)*numVoxelsPerDim
		+ coarseZ*pow(numVoxelsPerDim,2);
	int ind111 = (coarseX + 1) + (coarseY + 1)*numVoxelsPerDim
		+ (coarseZ + 1)*pow(numVoxelsPerDim,2);

	double c00 = map[ind000]*(1-xd) + map[ind100]*xd;
	double c01 = map[ind001]*(1-xd) + map[ind101]*xd;
	double c10 = map[ind010]*(1-xd) + map[ind110]*xd;
	double c11 = map[ind011]*(1-xd) + map[ind111]*xd;

	//printf("%f %f %f %f\n", map[ind000], map[ind100], map[ind001], map[ind101]);

	double c0 = c00*(1-yd) + c10*yd;
	double c1 = c01*(1-yd) + c11*yd;

	double c = c0*(1-zd) + c1*zd;

	return c;
}