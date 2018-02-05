#include "mpi.h"
#include <stdio.h>

/*******************************************************************/
/*             Update Bounds of subdomain with me process          */
/*******************************************************************/

void updateBound(double** x, int NeighBor[], MPI_Comm comm2d, MPI_Datatype column_type, int me, int* xs, int* ys,
                 int* xe, int* ye, int xcell)
{
   int S = 0, E = 1, N = 2, W = 3;
   int flag;
   MPI_Status status;

   /********* North/South communication **********************************/

   flag = 1;
   /* Send my boundary to North and receive from South */
   MPI_Sendrecv(&x[ys[me]][xs[me]], xcell, MPI_DOUBLE, NeighBor[N], flag, &x[ye[me]+1][xs[me]], xcell,
                MPI_DOUBLE, NeighBor[S], flag, comm2d, &status);

   /* Send my boundary to South and receive from North */
   MPI_Sendrecv(&x[ye[me]][xs[me]], xcell, MPI_DOUBLE, NeighBor[S], flag, &x[ys[me]-1][xs[me]], xcell,
                MPI_DOUBLE, NeighBor[N], flag, comm2d, &status);

   /********* Est/West communication *************************************/

   flag = 2;
   /* Send my boundary to Est and receive from West */
   MPI_Sendrecv(&x[ys[me]][xe[me]], 1, column_type, NeighBor[E], flag, &x[ys[me]][xs[me]-1], 1, column_type,
                NeighBor[W], flag, comm2d, &status);

   /* Send my boundary to West and receive from Est */
   MPI_Sendrecv(&x[ys[me]][xs[me]], 1, column_type, NeighBor[W], flag, &x[ys[me]][xe[me]+1], 1, column_type,
                NeighBor[E], flag, comm2d, &status);
}
