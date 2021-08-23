CC = mpicc
cflags = -O3 -Wall
libraries = -lm

objects = jacobi_mpi.o helpers.o
executable = jacobi_mpi.x

all: $(objects)
	$(CC) $(cflags) $(objects) -o $(executable) $(libraries)

helpers.o: helpers.c
	$(CC) $(cflags) -c helpers.c -o helpers.o

jacobi_mpi.o: jacobi_mpi.c
	$(CC) $(cflags) -c jacobi_mpi.c -o jacobi_mpi.o

clean:
	rm $(objects) $(executable)

wash:
	rm JseqJacobi.o* JseqJacobi.e*
