CC = mpicc
cflags = -O3 -Wall
libraries = -lm


all:
	$(CC) $(cflags) jacobi_serial.c -o jacobi_serial.x $(libraries)
clean:
	rm jacobi_serial.x
