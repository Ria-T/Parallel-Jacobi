/* When called the function calculates:
   The size of the local table given the rank 
   The corespoding coordinates of the local/small table in the big/whole table
   The topology of a procces in the mpi_cart
   Creates the cartesian topology */
void get_local_table(int *n, int *m, int** topology_dims, int ***coordinates, int *rank, int world_size, MPI_Comm *cart_comm);

/* Given the deltax,y and the xLeft, xRight of the "big" table,
   the dimensions of the local/"small" table and the topology
   it calculates the coresponding values for the local/"small" table*/
void calculate_xy_range(double start, double end, double *a, double *b, double delta, int start_coordinates, int end_coordinates);


#define UP 0
#define LEFT 2
#define RIGHT 3
#define DOWN 1
/* Given a cartesian communicator it allocates and returns a table with the ranks of the neighbors
   The correct neighbors indexes (up, left  etc.) can be read using the above definitions*/
int *get_neighbors(MPI_Comm *cart_comm);