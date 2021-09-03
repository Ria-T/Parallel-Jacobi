CC = mpicc
cflags = -O3 -Wall
libraries = -lm

objects = jacobi_mpi.o helpers.o partitioning.o
executable = jacobi_mpi.x

all: $(objects)
	$(CC) $(cflags) $(objects) -o $(executable) $(libraries)

partitioning.o: partitioning.c partitioning.h
	$(CC) $(cflags) -c partitioning.c -o partitioning.o

helpers.o: helpers.c helpers.h
	$(CC) $(cflags) -c helpers.c -o helpers.o

jacobi_mpi.o: jacobi_mpi.c
	$(CC) $(cflags) -c jacobi_mpi.c -o jacobi_mpi.o

clean:
	rm $(objects) $(executable)

wash:
	rm Jacobi.o* Jacobi.e*
	rm tabl*
