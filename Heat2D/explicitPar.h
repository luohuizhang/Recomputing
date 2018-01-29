void readParam(int* iconf, double* conf);

void processToMap(int* xs, int *xe, int *ys, int *ye, int xcell, int ycell, int x_domains, int y_domains);

void initValues( double** tab, int a, int b, double temp1, double temp2);

void updateBound( double** x, int* NeighBor, MPI_Comm comm, MPI_Datatype datatype, int current, int* x1, int* y1, int* x2, int* y2, int sizex);

void Explicit( double** x0, double** x, double dt, double* res, double hx, double hy, int me, int* x1, int* y1, int* x2, int* y2, double k);

void initValues2(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init, int nproc, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains);
void initValues3(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init, int nproc, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains);
