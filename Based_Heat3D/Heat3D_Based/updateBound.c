#include "mpi.h"

void updateBound(double*** x, int size_tot_x, int size_tot_y, int size_tot_z, 
                 int* neighbor, MPI_Comm comm3d, MPI_Datatype matrix_type_oxz, MPI_Datatype 
                 matrix_type_oxy, MPI_Datatype matrix_type_oyz , int me, int* xs, int* ys, int* zs, 
                 int* xe, int* ye, int* ze) {

   /* Local variables */
   int S=0, E=1, N=2, W=3, Zd=4, Zu=5;
   int flag;
   MPI_Status status;

   /********* North/South communication ************************************/
   flag = 1;
   /* Send my boundary to North and receive from South */
   MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oxz ,neighbor[N], flag, &x[xs[me]][ye[me]+1][zs[me]], 1, 
                matrix_type_oxz, neighbor[S], flag, comm3d, &status);

   /* Send my boundary to South and receive from North */
   MPI_Sendrecv(&x[xs[me]][ye[me]][zs[me]], 1, matrix_type_oxz, neighbor[S], flag, &x[xs[me]][ys[me]-1][zs[me]], 1, 
                matrix_type_oxz, neighbor[N], flag, comm3d, &status);

   /********* Est/West communication ***************************************/
   flag = 2;
   /* Send my boundary to Est and receive from West */
   MPI_Sendrecv(&x[xe[me]][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[E], flag, &x[xs[me]-1][ys[me]][zs[me]], 1, 
                matrix_type_oyz, neighbor[W], flag, comm3d, &status);

   /* Send my boundary to West and receive from Est */
   MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[W], flag, &x[xe[me]+1][ys[me]][zs[me]], 1, 
                matrix_type_oyz, neighbor[E], flag, comm3d, &status);

   /********* Zdown/Zup communication **************************************/
   flag = 3;
   /* Send my boundary to Zup and receive from Zdown */
   MPI_Sendrecv(&x[xs[me]][ys[me]][ze[me]], 1, matrix_type_oxy, neighbor[Zu], flag, &x[xs[me]][ys[me]][zs[me]-1], 1, 
                matrix_type_oxy, neighbor[Zd], flag, comm3d, &status);

   /* Send my boundary to Zdown and receive from Zup */
   MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oxy, neighbor[Zd], flag, &x[xs[me]][ys[me]][ze[me]+1], 1, 
                matrix_type_oxy, neighbor[Zu], flag, comm3d, &status);
}
