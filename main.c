#include <stdio.h>
#include "data_handler.c"

int main(){
	int numGals = getNumGals();
	galaxy *gals = readData();

	printf("%d\n", numGals);

	cartesianGalaxy *cartGals = convertToCartesian(gals, numGals);

	int i;/*
	for(i = 0; i < 10; i++){
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

	galaxy *newGals = convertFromCartesian(cartGals,gals,numGals);

	FILE *fp;
	fp = fopen("newGals.csv", "w");

	for(i = 0; i < numGals; i++){
		fprintf(fp, "%e, %e, %e\n", newGals->ra, newGals->dec,
			newGals->z_red);
		newGals++;
	}
	fclose(fp);
}
