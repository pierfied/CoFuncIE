#include <stdio.h>
#include "data_handler.h"
#include "galaxy_structs.h"
#include "density_map.h"
#include "map_likelihood.h"
#include "hamiltonian.h"
#include "omp.h"
#include "radial_likelihood.h"

int main(){

	omp_set_num_threads(16);

	/*double x = 1, p = 0, M = 1, epsilon = 0.1;
	int N = 1, n = 1000;
	double (*force)(double,int,void*);
	force = &springForce;
	hoParams params;
	params.k = 1;

	leapfrogIntegrator(&x, &p, &M, N, n, epsilon, force, &params);*/

	/*
	double x[10], p[10], M[10], k[10];
	double epsilon = 0.1;
	int N = 10, n = 1000;
	hoParams params;
	double (*force)(double*,int,void*);
	force = &springForce;

	int i;
	for(i = 0; i < N; i++){
		x[i] = 1;
		p[i] = 0;
		M[i] = 1;
		k[i] = i;
	}

	params.k = k;

	leapfrogIntegrator(x, p, M, N, n, epsilon, force, &params);

	exit(0);*/



	int numGals;
	galaxy *gals = readData(&numGals);

	printf("%d\n", numGals);

	cartesianGalaxy *cartGals = convertToCartesian(gals, numGals);
	

	int i;
	/*for(i = 0; i < 10; i++){
		printf("x: %e\ty: %e\tz: %e\n", cartGals->x, cartGals->y, cartGals->z);
		cartGals++;
	}

	double maxZ = (gals++)->z_red;
	double minZ = maxZ;
	double avgZ = maxZ;
	for(i = 1; i < numGals; i++){
		if(gals->z_red < minZ) minZ = gals->z_red;
		if(gals->z_red > maxZ) maxZ = gals->z_red;
		avgZ += gals->z_red;
		gals++;
	}

	printf("Max: %e\nMin: %e\nAvg: %e\n", maxZ, minZ, avgZ/numGals);*/

	//galaxy *newGals = convertFromCartesian(cartGals,gals,numGals);
	
	double *r, *z;
	int numInterpPts;
	setupRedshiftInterp(&r, &z, &numInterpPts);
	printf("Query = %f\n", interpRedshift(500, r, z, numInterpPts));
	printf("Query = %f\n", interpDist(0.033, r, z, numInterpPts));

	int xStart = 600;
	int yStart = 0;
	int zStart = 0;
	int numVoxelsPerDim = 10;
	int boxLength = 400;

	galaxy *trimmedGals = trimGalaxyList(gals, &numGals, xStart, yStart,
		zStart, boxLength);

	int wob;
	double maxZ = 0;
	double minZ = 10;
	for(wob = 0; wob < numGals; wob++){
		if(trimmedGals[wob].z_red > maxZ){
			maxZ = trimmedGals[wob].z_red;
		}

		if(trimmedGals[wob].z_red < minZ){
			minZ = trimmedGals[wob].z_red;
		}
	}

	printf("Min z: %f\nMax z: %f\n", minZ, maxZ);

	printf("Here!\n");

	galaxy *gal2 = trimmedGals+10000;
	cartesianGalaxy *gal = convertToCartesian(gal2, 1);
	printf("x: %f\ty: %f\tz: %f\n", gal->x, gal->y, gal->z);

	int *voxels;
	double *map = generateMap(trimmedGals, numGals, numVoxelsPerDim,
		xStart, yStart, zStart, boxLength, &voxels);

	mapData md;
	md.map = map;
	md.xStart = xStart;
	md.yStart = yStart;
	md.zStart = zStart;
	md.numVoxelsX = numVoxelsPerDim;
	md.numVoxelsY = numVoxelsPerDim;
	md.numVoxelsZ = numVoxelsPerDim;
	md.boxLengthX = boxLength;
	md.boxLengthY = boxLength;
	md.boxLengthZ = boxLength;

	drawRSample(gal2, md);

	double dist = sqrt(gal->x*gal->x + gal->y*gal->y + gal->z*gal->z);

	printf("r: %f\n", dist);

	double xUnit = gal->x/dist;
	double yUnit = gal->y/dist;
	double zUnit = gal->z/dist;

	printf("x: %f\ty: %f\tz: %f\n", xUnit, yUnit, zUnit);

	double Mo = cosmo_M;
	double k = cosmo_k;
	double lambda = cosmo_lambda;
	double h = cosmo_h;
	double Dh = 3000/h; // [Mpc]

	double Ez = 1./sqrt(Mo*pow(1+gal2->z_red,3)
		       + k*pow(1+gal2->z_red,2) + lambda);

	double sigR = Dh*Ez*gal2->z_err;

	printf("sz: %f\tsr: %f\n", gal2->z_err, sigR);

	FILE *fp1;
	fp1 = fopen("radialLikelihood","w");

	int resoultion = 10;

	double *fMap = fineMap(map, numVoxelsPerDim, resoultion);

	double R;
	double sum = 0;
	for(R = 1; R < 1500; R++){
		double x = R * xUnit;
		double y = R * yUnit;
		double z = R * zUnit;

		if(x < xStart || x > xStart + boxLength || y < yStart ||
			y > yStart + boxLength || z < zStart || 
			z > zStart + boxLength){
			continue;
		}

		double voxelLength = boxLength/(double)(numVoxelsPerDim*resoultion);

		int a = (x - xStart)/voxelLength;
		int b = (y - yStart)/voxelLength;
		int c = (z - zStart)/voxelLength;
		int index = a + numVoxelsPerDim*resoultion*b + pow(numVoxelsPerDim*resoultion, 2)*c;

		double gaussR = 1/(sigR * sqrt(2*3.14159)) * exp(-0.5*pow((R-dist)/sigR,2));

		double mapyStuff = (1+fMap[index]);

		sum += gaussR;

		int coarsex = a/resoultion;
		int coarsey = b/resoultion;
		int coarsez = c/resoultion;

		//printf("%d %d %d %f\n", coarsex, coarsey, coarsez, mapyStuff);

		printf("gauss(R = %f) = %f\tsum = %f\t1+delta = %f\n", R, gaussR, sum, mapyStuff);
		fprintf(fp1,"%f,%f,%f\n",R,gaussR, mapyStuff);

		//printf("x: %f\ty: %f\tz: %f\tR: %f\n", x, y, z, R);
	}
	fclose(fp1);

	FILE *fp2;
	fp2 = fopen("fine","w");
	int woop;
	for(i = 0; i < numVoxelsPerDim*resoultion; i++){
		for(woop = 0; woop < numVoxelsPerDim*resoultion; woop++){
			int ind = i + woop*numVoxelsPerDim*resoultion
				+ 10*pow(numVoxelsPerDim*resoultion,2);
			fprintf(fp2,"%f,",fMap[ind]);
		}
		fprintf(fp2, "\n");
	}

	FILE *fp3;
	fp3 = fopen("coarse","w");
	for(i = 0; i < numVoxelsPerDim; i++){
		for(woop = 0; woop < numVoxelsPerDim; woop++){
			int ind = i + woop*numVoxelsPerDim + pow(numVoxelsPerDim,2);
			fprintf(fp3,"%f,",map[ind]);
		}
		fprintf(fp3, "\n");
	}

	FILE *fp4;
	FILE *fp5;
	fp4 = fopen("fMap","w");
	fp5 = fopen("cMap","w");
	int a;
	for( a = 0; a < numVoxelsPerDim*resoultion; a++){
		int b;
		for(b = 0; b < numVoxelsPerDim*resoultion; b++){
			int c;
			for(c = 0; c < numVoxelsPerDim*resoultion; c++){
				int ind = a + b*numVoxelsPerDim*resoultion
					+ c*pow(numVoxelsPerDim*resoultion,2);

				int coarsex = a/resoultion;
				int coarsey = b/resoultion;
				int coarsez = c/resoultion;

				int coarseind = coarsex + coarsey*numVoxelsPerDim
					+ coarsez*pow(numVoxelsPerDim,2);

				fprintf(fp4,"%d,%d,%d,%f\n",a,b,c,fMap[ind]);
				fprintf(fp5,"%d,%d,%d,%f\n",a,b,c,map[coarseind]);
			}
		}
	}

	exit(0);

	// Declare and initialize the xi samples.
	int numSamps = 5;
	double rSamp[] = {0, 250, 500, 750, 1000};
	double xiSamp[] = {0.1, 0.25, 0.5, 0.75, 1};

	gsl_spline *spline = initCorrSpline(numSamps, rSamp, xiSamp);

	double *cov = generateCov(numVoxelsPerDim, boxLength, spline);

	for(i = 0; i < 10; i++){
		printf("%f\n", exp(*(cov+i))-1);
	}
	//exit(0);

	FILE *fp;
	fp = fopen("Catalogues/cov.txt", "w");
	int n = pow(numVoxelsPerDim, 3);
	for(i = 0; i < n; i++){
		int j;
		for(j = 0; j < n; j++){
			int index = i + n*j;
			fprintf(fp, "%e ", *(cov + index));
		}
		fprintf(fp, "\n");
	}

	double *invCov = invertCov(cov, numVoxelsPerDim);

	fp = fopen("Catalogues/inv_cov.txt", "w");
	n = pow(numVoxelsPerDim, 3);
	for(i = 0; i < n; i++){
		int j;
		for(j = 0; j < n; j++){
			int index = i + n*j;
			fprintf(fp, "%e ", *(cov + index));
		}
		fprintf(fp, "\n");
	}

	double lnLikeMap = mapLnLikelihood(map, voxels, numVoxelsPerDim,
		boxLength, spline);

	printf("Likelihood: %f\n", lnLikeMap);

	map = modifyMap(cov, invCov, map, voxels, numVoxelsPerDim);

	lnLikeMap = mapLnLikelihood(map, voxels, numVoxelsPerDim,
		boxLength, spline);

	printf("Likelihood: %f\n", lnLikeMap);

	exit(0);

	// Create a uniform map to test the likelihood.
	for(i = 0; i < n; i++){
		*(map+i) = 0;
		*(voxels+i) = 4;
	}

	// Calculate the likelihood of the flat map.
	lnLikeMap = mapLnLikelihood(map, voxels, numVoxelsPerDim, boxLength,
		spline);

	printf("Flat Map Ln Likelihood: %f\n", lnLikeMap);

	// Calculate the mass matrix.
	double *M;
	M = generateMassMatDiag(invCov, voxels, numVoxelsPerDim);
	printf("M = %f\n", M[0]);
	printf("M = %f\n", M[12]);
	printf("M = %f\n", M[122]);
	printf("M = %f\n", M[555]);
	printf("M = %f\n\n\n", M[n-1]);

	double *invM;
	invM = invertDiagMat(M, numVoxelsPerDim);
	M = invM;
	printf("M = %f\n", M[0]);
	printf("M = %f\n", M[12]);
	printf("M = %f\n", M[122]);
	printf("M = %f\n", M[555]);
	printf("M = %f\n\n\n", M[n-1]);

	double *p;
	p = generateMomenta(numVoxelsPerDim);

	for(i = 0; i < 10; i++){
		printf("p = %f\n", p[i]);
	}

	/*
	FILE *fp;
	fp = fopen("Catalogues/trimmedGals.csv", "w");

	for(i = 0; i < numGals; i++){
		fprintf(fp, "%e, %e, %e\n", trimmedGals->ra, trimmedGals->dec,
			trimmedGals->z_red);
		trimmedGals++;
	}
	fclose(fp);
	*/
}
