void readParam(int* iconf, double* conf);

void processToMap(int me, int* xs, int *xe, int *ys, int *ye, int *zs, int *ze, 
                  int xcell, int ycell, int zcell, int x_domains, int y_domains, int z_domains);

void initValues(int nb, double*** tab, int a, int b, int c, double temp1, double temp2);

void updateBound(double*** x, int sizex, int sizey, int sizez, int* NeighBor, MPI_Comm comm, 
                 MPI_Datatype type1, MPI_Datatype type2, MPI_Datatype type3, int current, 
                 int* x1, int* y1, int* z1, int* x2, int* y2, int* z2);

void computeNext(double*** x0, double*** x, int sizex, int sizey, int sizez, double dt, 
                 double hx, double hy, double hz, double* res, int me, int* x1, int* y1, int* z1, 
                 int* x2, int* y2, int* z2, double k);
