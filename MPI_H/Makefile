CC = mpicc
cflags = -fopenmp -O3 -g -Wall
libraries = -lm

objects = jacobi_mpi.o partitioning.o
executable = jacobi_mpi.x

all: $(objects)
	$(CC) $(cflags) $(objects) -o $(executable) $(libraries)

partitioning.o: partitioning.c partitioning.h
	$(CC) $(cflags) -c partitioning.c -o partitioning.o

jacobi_mpi.o: jacobi_mpi.c
	$(CC) $(cflags) -c jacobi_mpi.c -o jacobi_mpi.o

sub: all
	qsub PBSseqJ.sh

clean:
	rm $(objects) $(executable)

