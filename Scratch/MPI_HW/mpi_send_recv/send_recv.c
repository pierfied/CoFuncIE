#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void main(){
	MPI_Init(NULL, NULL);

	//Get process rank, and world size info.
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int number;
	if(world_rank == 0){
		number = -1;
		for(number = -1; -number < world_size; number--){
			MPI_Send(&number, 1, MPI_INT, -number, 0,
				MPI_COMM_WORLD);
		}
	}else{
		MPI_Status status;
		MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
			&status);
		printf("Process %d received number %d from process %d.",
			world_rank, number, status.MPI_SOURCE);
		printf("\n");
	}

	MPI_Finalize();
}