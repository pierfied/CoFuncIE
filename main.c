#include <stdio.h>
#include "data_handler.c"

int main(){
	int numGals = getNumGals();
	galaxy *gals = readData();

	printf("%d\n", numGals);

	cartesianGalaxy *cartGals = convertCartesian(gals, numGals);

	int i;
	for(i = 0; i < 10; i++){
		printf("x: %e\ty: %e\tz: %e\n", cartGals->x, cartGals->y, cartGals->z);
		cartGals++;
	}
}
