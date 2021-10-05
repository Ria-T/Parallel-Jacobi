void init_debug();

void draw_table(int bold[2][2], int size[2], int rank);
void write_table(double *table, int size[2], int rank, int *neighbors);
int printr(int rank, const char* fmt, ...);