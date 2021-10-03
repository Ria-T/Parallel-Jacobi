CC = mpicc
cflags = -O3 -g -Wall
libraries = -lm

objects = jacobi_mpi.o helpers.o partitioning.o table_commuting.o
executable = jacobi_mpi.x

all: $(objects)
	$(CC) $(cflags) $(objects) -o $(executable) $(libraries)

partitioning.o: partitioning.c partitioning.h
	$(CC) $(cflags) -c partitioning.c -o partitioning.o

helpers.o: helpers.c helpers.h
	$(CC) $(cflags) -c helpers.c -o helpers.o

jacobi_mpi.o: jacobi_mpi.c
	$(CC) $(cflags) -c jacobi_mpi.c -o jacobi_mpi.o

table_commuting.o: table_commuting.c
	$(CC) $(cflags) -c table_commuting.c -o table_commuting.o

sub: all wash
	qsub PBSseqJ.sh

run: all
	mpiexec -n 1 gdb $(executable) : -n 3 $(executable)

clean:
	rm $(objects) $(executable)

wash:
	rm Jacobi.o* Jacobi.e* tabl* valtabl*
