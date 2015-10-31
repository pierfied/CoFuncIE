all: mpi_convert_cart

mpi_convert_cart: mpi_convert_cart.c
	mpicc mpi_convert_cart.c -o mpi_convert_cart -lm

run: mpi_convert_cart.c
	mpicc mpi_convert_cart.c -o mpi_convert_cart -lm
	mpirun -np 9 ./mpi_convert_cart partial_data.txt

clean:
	rm mpi_convert_cart