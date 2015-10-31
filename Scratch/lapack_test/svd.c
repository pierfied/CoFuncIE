#include <stdio.h>
#include <mkl_lapacke.h>
#include <math.h>
#include <time.h>

double *importCov(char *filename, int *n);
void print_matrix(FILE *f, double *m, int n);

void main(){
	char filename[100] = "cov.txt";

	// Load the covariance matrix.
	int n;
	double *cov = importCov(filename, &n);


	printf("Starting\n");
	time_t start = time(NULL);

	// Compute the tridagonal form of the matrix.
	double *D = malloc(n * sizeof(double));
	double *E = malloc((n-1) * sizeof(double));
	double *TAU = malloc((n-1) * sizeof(double));
	int info = LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',n,cov,n,D,E,TAU);

	// Compute the eigenvalues of the matrix.
	int LDZ = n;
	double *Z = malloc(LDZ * n * sizeof(double));
	info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR,'n',n,D,E,Z,LDZ);

	time_t end = time(NULL);
	printf("Finished in %ld second(s).\n", end-start);
	
	printf("Computation completed with code: %d\n", info);


	// Copy the upper values to the lower diagonal.
	int i,j,x,y;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			int a = n*i + j;
			int b = i + n*j;
			*(cov + b) = *(cov + a);
		}
	}

	// Output the covariance matrix.
	FILE *fp;
	fp = fopen("eigenvals.txt", "w");
	print_matrix(fp, D, n);
}

double *importCov(char *filename, int *n){
	// Open the file.
	FILE *fp;
	fp = fopen(filename, "r");

	// Get the number of lines in the file.
	char c;
	int numLines = 0;
	while((c = getc(fp)) != EOF){
		if(c == '\n'){
			numLines++;
		}
	}
	rewind(fp);

	*n = numLines;

	// Allocate the catalog of galaxies.
	double *cov = malloc(pow(numLines,2)*sizeof(double));

	// Load each galaxy.
	int i;
	double *curCovVal = cov;
	for(i = 0; i < pow(numLines,2); i++){
		fscanf(fp, "%le", curCovVal);

		curCovVal++;
	}

	return cov;
}

void print_matrix(FILE *f, double *m, int n){
	int i;
	double min = *m;
	double max = *m;
	for (i = 0; i < n; i++) {
		if(*m > max) max = *m;
		if(*m < min) min = *m;
		fprintf(f, "%e \n", *m++);
	}
	fprintf(f, "Min: %e\n", min);
	fprintf(f, "Max: %e\n", max);
}
