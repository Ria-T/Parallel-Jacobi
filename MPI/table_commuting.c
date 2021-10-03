#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

int recieve_big_table(double **u, int small_size[2], int actual_size[2], int world_size, MPI_Comm *communicator, MPI_Datatype *table_row){
    int proc,one_table_size,i,j,y,istep,jstep;
    double *big_u = (double*)calloc((actual_size[0]+2)*(actual_size[1]+2), sizeof(double));
    if(big_u == NULL) return -1;

    istep = 0;
    jstep = 1;
    i = 1;
    j = 1;

    for(y=1; y<=small_size[1]; y++){
        memcpy(
            big_u + j*(actual_size[0]+2) + i,
            *u + j*(small_size[0]+2) + i,
            small_size[0]*sizeof(double)
        );
        i += istep;
        j += jstep;
    }
    //memcpy( big_u, *u, (small_size[0]+2)*(small_size[1]+1)*sizeof(double));
    free(*u);

    istep = 0;
    jstep = j;

    for(proc=1; proc<world_size; proc++){
        MPI_Recv(
                    big_u +  j*(actual_size[0]+2) + i,
                    small_size[1],
                    *table_row,
                    proc,
                    0,
                    *communicator,
                    MPI_STATUS_IGNORE
        );
        printf("%d recivied\n",proc);
        i += istep;
        j += jstep;
    }

    *u = big_u;

    return 0;
}

void send_table(double *u, int small_size[2], MPI_Comm *communicator, MPI_Datatype *table_row){
    MPI_Send(
        u,
        (small_size[0]+2)*(small_size[1]+2),
        *table_row,
        0,
        0,
        *communicator
    );
}