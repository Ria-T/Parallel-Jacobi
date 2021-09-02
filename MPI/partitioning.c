#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "helpers.h"

void calculate_topology(int size[2], int processes, int dims[2]){
    double t;
    if( processes > 2 ){
        t = sqrt(processes);
        if ( t == (int)t ){
            dims[0] = t;
            dims[1] = t;
        }else{
            dims[0] = 1;
            dims[1] = processes;
        }
    }else{
        dims[0] = processes;
        dims[1] = dims[0];
    }
}

#define UP 0
#define LEFT 2
#define RIGHT 3
#define DOWN 1

int *get_neighbors(MPI_Comm *cart_comm){
    int *neighbors = malloc(sizeof(int)*4);
    MPI_Cart_shift(*cart_comm, 0, 1, neighbors + UP, neighbors + DOWN);
    MPI_Cart_shift(*cart_comm, 1, 1, neighbors + LEFT , neighbors + RIGHT);
    return neighbors;
}

void get_local_table(int size[2], int coords[2], int *rank, int world_size){
    int ndims = 2, dims[2], periods[2] = {0, 0};

    if(*rank == 0){
        calculate_topology(size, world_size, dims);
        printf("There will be %d tables of %dx%d with a %dx%d topology\n",world_size,size[0]/dims[0],size[1]/dims[1],dims[0],dims[1]);
    }
    MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Comm cart_comm;
    MPI_Cart_create( MPI_COMM_WORLD , 2, dims, periods, 1, &cart_comm);

    int *neighbors = get_neighbors(&cart_comm);

    //MPI_Cart_create(comm,ndim,dimensions,periodic,1,&comm2d);
    int cart_rank;
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm,cart_rank,2,coords);

    printf("I am %d: (%d,%d); originally %d\n   %2d\n%2d %2d %2d\n   %2d\n_\n"
    ,cart_rank, coords[0], coords[1], *rank
    ,neighbors[UP], neighbors[LEFT], cart_rank, neighbors[RIGHT], neighbors[DOWN]);

    free(neighbors);

    int coordinates[0][0];

}