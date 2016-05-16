double *fineMap(double *map, int numVoxelsPerDim, int resoultion){
	int newNumVPD = numVoxelsPerDim * resoultion;

	double *dx = malloc(sizeof(double) * pow(numVoxelsPerDim,3));
	double *dy = malloc(sizeof(double) * pow(numVoxelsPerDim,3));
	double *dz = malloc(sizeof(double) * pow(numVoxelsPerDim,3));

	int i;
	#pragma omp parallel for
	for(i = 0; i < numVoxelsPerDim; i++){
		int j;
		for(j = 0; j < numVoxelsPerDim; j++){
			int k;
			for(k = 0; k < numVoxelsPerDim; k++){
				int ind = i + numVoxelsPerDim*j
					+ pow(numVoxelsPerDim,2)*k;


				int ind2;

				if(i < numVoxelsPerDim - 1){
					ind2 = (i+1) + numVoxelsPerDim*j
						+ pow(numVoxelsPerDim,2)*k;

					dx[ind] = (map[ind2] - map[ind])/resoultion;
				}else{
					ind2 = (i-1) + numVoxelsPerDim*j
						+ pow(numVoxelsPerDim,2)*k;

					dx[ind] = dx[ind2];
				}

				if(j < numVoxelsPerDim - 1){
					ind2 = i + numVoxelsPerDim*(j+1)
						+ pow(numVoxelsPerDim,2)*k;

					dy[ind] = (map[ind2] - map[ind])/resoultion;
				}else{
					ind2 = i + numVoxelsPerDim*(j-1)
						+ pow(numVoxelsPerDim,2)*k;

					dy[ind] = dy[ind2];
				}

				if(j < numVoxelsPerDim - 1){
					ind2 = i + numVoxelsPerDim*j
						+ pow(numVoxelsPerDim,2)*(k+1);

					dz[ind] = (map[ind2] - map[ind])/resoultion;
				}else{
					ind2 = i + numVoxelsPerDim*j
						+ pow(numVoxelsPerDim,2)*(k-1);

					dz[ind] = dz[ind2];
				}
			}
		}
	}

	double *fMap = malloc(sizeof(double) * pow(newNumVPD,3));

	#pragma omp parallel for
	for(i = 0; i < newNumVPD; i++){
		int j;
		for(j = 0; j < newNumVPD; j++){
			int k;
			for(k = 0; k < newNumVPD; k++){
				int fineInd = i + newNumVPD*j + pow(newNumVPD,2)*k;

				int x = i/resoultion;
				int y = j/resoultion;
				int z = k/resoultion;

				int coarseInd = x + numVoxelsPerDim*y 
					+ pow(numVoxelsPerDim,2)*z;

				double deltaX = i - x*resoultion;
				double deltaY = j - y*resoultion;
				double deltaZ = k - z*resoultion;

				fMap[fineInd] = map[coarseInd] + deltaX*dx[coarseInd] 
					+ deltaY*dy[coarseInd] + deltaZ*dz[coarseInd];
			}
		}
	}

	free(dx);
	free(dy);
	free(dz);



	return fMap;
}