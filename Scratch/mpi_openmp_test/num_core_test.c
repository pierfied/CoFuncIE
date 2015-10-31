#include <stdio.h>
#include <mpi.h>

void main(){
	MPI_Init(NULL, NULL);

	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	int num_cores = omp_get_num_procs();

	printf("Number of available cores on node %s: %d\n",
		processor_name, num_cores);

	MPI_Finalize();
}
