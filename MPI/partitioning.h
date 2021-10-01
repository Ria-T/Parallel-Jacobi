//For neighbors table

#define UP 0
#define LEFT 2
#define RIGHT 3
#define DOWN 1

int calculate_range(double *a, double *b, int n, int processes);
void get_local_table(int *n, int *m, int** topology_dims, int ***coordinates, int *rank, int world_size, MPI_Comm *cart_comm);
int *get_neighbors(MPI_Comm *cart_comm);