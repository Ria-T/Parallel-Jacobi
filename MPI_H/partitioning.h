//For neighbors table

#define UP 0
#define LEFT 2
#define RIGHT 3
#define DOWN 1

void calculate_range(double start, double end, double *a, double *b, double delta, int start_coordinates, int end_coordinates);
void get_local_table(int *n, int *m, int** topology_dims, int ***coordinates, int *rank, int world_size, MPI_Comm *cart_comm);
int *get_neighbors(MPI_Comm *cart_comm);