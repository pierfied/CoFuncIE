#include <stdio.h>
#include <stdlib.h>

void main(){
	int *x = malloc(sizeof(int) * 100);

	int i;
	for(i = 0; i < 10; i++){
		*(x+i) = i;
	}

	for(i = 0; i < 10; i++){
		printf("i = %d\n", (x+i+1));
	}
}