#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "helpers.h"


int calculate_range(double *a, double *b, int n, int processes){
    n++;
    if(*b <= *a) return -1;
    double d = (*b - *a)/((double)processes);
    if(n==1){
        //*a = *a;
        *b = *a + d;
    }else if(n==processes){
        *a = *a + (n-1)*d;
        //*b = *b;
    }else{
        //[a+δ,a+2δ]
        *a = *a+d*(n-1);
        *b = *a+d*(n-1);
    }
    return 0;
}

void calculate_topology(int processes, int (*dims)[2]){
    double t;
    if( processes > 2 ){
        t = sqrt(processes);
        if ( t == (int)t ){
            (*dims)[0] = t;
            (*dims)[1] = t;
        }else{
            (*dims)[0] = 1;
            (*dims)[1] = processes;
        }
    }else{
        (*dims)[0] = processes;
        (*dims)[1] = (*dims)[0];
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

void get_local_table(int *n, int *m, int ***coordinates, int *rank, int world_size, MPI_Comm *cart_comm){
    int dims[2], periods[2] = {0, 0}, size[2], coords[2], bold[2][2];

    if(*rank == 0){
        calculate_topology(world_size, &dims);
        printf("There will be %d tables of %dx%d with a %dx%d topology\n",world_size,*n/dims[0],*m/dims[1],dims[0],dims[1]);
    }
    MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Cart_create( MPI_COMM_WORLD , 2, dims, periods, 1, cart_comm);

    int cart_rank;
    MPI_Comm_rank(*cart_comm, &cart_rank);
    MPI_Cart_coords(*cart_comm,cart_rank,2,coords);

    *coordinates = malloc(sizeof(int*)*2);
    (*coordinates)[0] = malloc(sizeof(int)*2); (*coordinates)[1] = malloc(sizeof(int)*2);

    (*coordinates)[0][0] = (*n/dims[0]) * coords[0] - 1;
    if((*coordinates)[0][0]!=0) (*coordinates)[0][0]++;
    (*coordinates)[0][1] = (*m/dims[1]) * coords[1] - 1;
    if((*coordinates)[0][1]!=0) (*coordinates)[0][1]++;


    (*coordinates)[1][0] = (*n/dims[0]) * (coords[0]+1) - 1;
    (*coordinates)[1][1] = (*m/dims[1]) * (coords[1]+1) - 1;
    //if(dims[1] == coords[1]+1) {printf("happend for %d %d\n",dims[0],cart_rank); coordinates[1][1] = size[1]}
    //if(dims[0] == coords[0]+1 && coordinates[1][1]!=size[0]) {printf("STUPID A (%d)\n",cart_rank);}
    if(dims[1] == coords[1]+1 && (*coordinates)[1][1]!=*m-1) {printf("STUPID B (%d)\n",cart_rank);}

    /*printf("I am %d: (%d,%d); originally %d\nWill have table [%d,%d] [%d,%d]\n   %2d\n%2d %2d %2d\n   %2d\n_\n"
    ,cart_rank, coords[0], coords[1], *rank
    ,coordinates[0][0],coordinates[0][1],coordinates[1][0],coordinates[1][1]
    ,neighbors[UP], neighbors[LEFT], cart_rank, neighbors[RIGHT], neighbors[DOWN]);*/

    printf("(%d,%d) %2d -> [%3d,%3d]-[%3d,%3d]\n",coords[0],coords[1],cart_rank
    ,(*coordinates)[0][0],(*coordinates)[0][1],(*coordinates)[1][0],(*coordinates)[1][1]);

    size[0] = *n;
    size[1] = *m;
    bold[0][0]=(*coordinates)[0][0];
    bold[0][1]=(*coordinates)[0][1];
    bold[1][0]=(*coordinates)[1][0];
    bold[1][1]=(*coordinates)[1][1];
    draw_table(bold, size, cart_rank);

    *rank = cart_rank;
}