CC = mpicc
cflags = -O3 -Wall
libraries = -lm

sourceFiles = jacobi_serial.c
executable = jacobi_serial.x

all:
	$(CC) $(cflags) $(sourceFiles) -o $(executable) $(libraries)

clean:
	rm $(executable)
