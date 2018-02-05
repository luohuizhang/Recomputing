void readParam(int* iconf, double* conf);

void processToMap(int* xs, int *xe, int *ys, int *ye, int xcell, int ycell, int x_domains, int y_domains,int *x_index, int* y_index);

void initValues( double** tab, int a, int b, double temp1, double temp2);

void updateBound( double** x, int* NeighBor, MPI_Comm comm, MPI_Datatype datatype, int current, int* x1, int* y1, int* x2, int* y2, int sizex);

void Explicit( double** x0, double** x, double dt, double* res, double hx, double hy, int me, int* x1, int* y1, int* x2, int* y2, double k);
