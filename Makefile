all:
	mpicc -O3 -lm jacobi_serial.c -o jacobi_serial.x