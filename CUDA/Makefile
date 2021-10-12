CC = nvcc
cflags = -O3
libraries = -lm

objects = jacobi_cuda.o
executable = jacobi_cuda.x

all: $(objects)
	$(CC) $(cflags) $(objects) -o $(executable) $(libraries)

jacobi_cuda.o: jacobi_cuda.cu
	$(CC) $(cflags) -c jacobi_cuda.cu -o jacobi_cuda.o

sub: all
	qsub PBScudaJ.sh

clean:
	rm $(objects) $(executable)
