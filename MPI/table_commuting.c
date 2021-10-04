#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

void get_neighbors_in_order(MPI_Comm *communicator, int ***topology_ranks, int *topology_dims){
    int cords[2];

    *topology_ranks = malloc(sizeof(int) * topology_dims[0]);

    for(int i=0; i<topology_dims[0]; i++)
        (*topology_ranks)[i] = malloc(sizeof(int) * topology_dims[1]);


    for(int i=0; i<topology_dims[0]; i++){
        for(int j=0; j<topology_dims[1]; j++){
            cords[0]=j;
            cords[1]=i;
            MPI_Cart_rank(*communicator, cords, &(*topology_ranks)[i][j]);
            //MPI_Cart_shift(*communicator, 1, 1, &neighbours[i][j], &dummy_neighbour);
            printf("%d, ",(*topology_ranks)[i][j]);
        }
        printf("\n");
    }

}

void cwrite_table(double *table, int size[2]){
    //#ifdef DEBUG
int rank = 0;
    char filename[50];
    sprintf(filename,"valtabl%d",rank);
    FILE * file = fopen (filename, "a+");
    fprintf(file,"(%d) table = %dx%d\n",rank,size[0],size[1]);
    //fprintf(file," | up=%d down=%d left=%d right=%d\n",neighbors[0],neighbors[1],neighbors[2],neighbors[3]);
    if(file<0) perror("open failed!\n");
    for(int j=0; j<size[1]+2; j++){
        for(int i=0; i<size[0]+2; i++){
            fprintf(file,"%.2f ",table[j*(size[0]+2)+i]);
        }
        fprintf(file,"\n");
    }
    fprintf(file,"-------------------\n");
    fclose(file);
}

int recieve_big_table(double **u, int small_size[2], int actual_size[2], int world_size, MPI_Comm *communicator, int *topology_dims){
    int proc,one_table_size,i,j,y,istep,jstep,**topology_ranks=NULL;
    get_neighbors_in_order(communicator,&topology_ranks,topology_dims);

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

    istep = actual_size[0]/topology_dims[0];
    jstep = actual_size[1]/topology_dims[1];//j;

    i = 1;
    j = 1;

    for(int t_i=0; t_i < topology_dims[0]; t_i++){
        i = 1;
        for(int t_j=0; t_j < topology_dims[1]; t_j++){

            if(topology_ranks[t_i][t_j] != 0){
                printf("%d reciving with %d,%d\n",topology_ranks[t_i][t_j],i,j);
                for(int row=0; row<small_size[1]; row++){
                    MPI_Recv(
                        big_u +  (row+j)*(actual_size[0]+2) + i,
                        small_size[0],
                        MPI_DOUBLE,
                        topology_ranks[t_i][t_j],
                        0,
                        *communicator,
                        MPI_STATUS_IGNORE
                    );
                }
                //cwrite_table(big_u,actual_size);
            }
            i += istep;
        }
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