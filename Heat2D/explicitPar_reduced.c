#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "explicitPar.h"

#define min(a,b) a <= b ? a : b

int main(int argc, char *argv[])
{
   /* Arrays */
   double **x0;
   double **x;
   double *xfinal;
   double *xtemp;
int nproc_reduced = atoi(argv[1]);
int nproc_full=8;
   /* Sizes of the discretization */
   double dt, dt1, dt2, hx, hy, epsilon;
   double resLoc, result, t, tstart, tend, elapsed_time;
int size_x_full, size_y_full,  x_domains_full, y_domains_full;
   /* Index variables */
   int i, j, k, l;

   /* Convergence */
   int convergence = 0;

   int step, maxStep;
   int size_x, size_y, me, x_domains,y_domains;
   int iconf[5];
   int size_x_glo, size_y_glo;
   double conf[2];
 int iconf_full[5];
        double conf_full[2];


   /* MPI variables */
   int nproc, ndims;
   MPI_Comm comm, comm2d;
   int dims[2];
   int periods[2];
   int reorganisation = 0;
   MPI_Datatype column_type;
   int S = 0, E = 1, N = 2, W = 3;
   int NeighBor[4];
   int xcell, ycell, size_tot_x, size_tot_y;
   int *xs, *ys, *xe, *ye;
   double temp1_init, temp2_init, k0;
   FILE* file;

   /* temp1_init: temperature init on edges */
   temp1_init = 10.0;

   /* temp2_init: temperature init inside */
   temp2_init = -10.0;

   /* Diffusivity coefficient */
   k0 = 1;

   MPI_Init(&argc, &argv);
int world_rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);

        int color;
        if(world_rank>=nproc_reduced)
                color=1;
        else
                color=0;
        MPI_Comm_split(MPI_COMM_WORLD,  color,  world_rank, &comm);

        MPI_Comm_size(comm,&nproc);
        //printf("Full: size is %d, world id is %d \n", nproc, world_rank);
        MPI_Comm_rank(comm,&me);

   MPI_Comm_size(comm,&nproc);
   MPI_Comm_rank(comm,&me);

   /* Getting input parameters */
   if(me==0)
{
readParam(iconf_full, conf_full); 
   readParam(iconf,conf);
}
   MPI_Bcast(iconf,5,MPI_INT,0,comm);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);
 MPI_Bcast(iconf_full,5,MPI_INT,0,comm);
   size_x    = iconf[0];
   size_y    = iconf[1];
   x_domains = iconf[2];
   y_domains = iconf[3];
   maxStep   = iconf[4];
   dt1       = conf[0];
   epsilon   = conf[1];
size_x_full    = iconf_full[0];
        size_y_full    = iconf_full[1];
        x_domains_full = iconf_full[2];
        y_domains_full = iconf_full[3];

   if((me==0) && (nproc!=(x_domains*y_domains)))
    printf("Number of processes not equal to Number of subdomains\n");

   size_x_glo = size_x + 2;
   size_y_glo = size_y + 2;
   hx = 1.0/(double)(size_x_glo);
   hy = 1.0/(double)(size_y_glo);
   dt2 = 0.25*((min(hx,hy))*(min(hx,hy)))/k0;

   /* Taking a good step for convergence */
   if(dt1>=dt2)
   {
    if(me==0)
    {
     printf("\n");
     printf("  Time step too large in 'param' file - Taking convergence criterion\n");
    }
    dt=dt2;
    }
   else dt=dt1;

   size_tot_x = size_x + 2*x_domains + 2;
   size_tot_y = size_y + 2*y_domains + 2;

   /* Allocation of Contiguous 2D arrays */
   xfinal = malloc(size_x*size_y*sizeof(*xfinal));

   x0 = malloc(size_tot_y*sizeof(*x0));
   x0[0] = malloc(size_tot_x*size_tot_y*sizeof(**x0));
   x = malloc(size_tot_y*sizeof(*x));
   x[0] = malloc(size_tot_x*size_tot_y*sizeof(**x));

   for(j=1;j<=size_tot_y-1;j++)
   {
    x0[j] = x0[0] + j*size_tot_x;
    x[j] = x[0] + j*size_tot_x;
   }

   xs = malloc(nproc*sizeof(int));
   xe = malloc(nproc*sizeof(int));
   ys = malloc(nproc*sizeof(int));
   ye = malloc(nproc*sizeof(int));
   int *x_index,*y_index;
        x_index= malloc(nproc*sizeof(int));
        y_index= malloc(nproc*sizeof(int));

   /* Create 2D cartesian grid */
   periods[0] = 0;
   periods[1] = 0;

   ndims = 2;
   dims[0] = x_domains;
   dims[1] = y_domains;

   MPI_Cart_create(comm, ndims, dims, periods, reorganisation, &comm2d);

   /* Identify neighbors */
   NeighBor[0] = MPI_PROC_NULL;
   NeighBor[1] = MPI_PROC_NULL;
   NeighBor[2] = MPI_PROC_NULL;
   NeighBor[3] = MPI_PROC_NULL;

   /* Left/West and right/Est neigbors */
   MPI_Cart_shift(comm2d, 0, 1, &NeighBor[W], &NeighBor[E]);

   /* Bottom/South and Upper/North neigbors */
   MPI_Cart_shift(comm2d, 1, 1, &NeighBor[S], &NeighBor[N]);

   /* Size of each cell */
   xcell = (size_x/x_domains);
   ycell = (size_y/y_domains);

   xtemp = malloc(xcell*ycell*sizeof(*xtemp));

   /* Calculate Map for xs and ys from "me" process  */
   processToMap(xs, xe, ys, ye, xcell, ycell, x_domains, y_domains,,x_index,y_index);

   /* Create column data type to communicate with East and West neighbors */
   MPI_Type_vector( ycell, 1, size_tot_x, MPI_DOUBLE, &column_type);
   MPI_Type_commit(&column_type);

   /* Initialization */
   initValues(x0, size_tot_x, size_tot_y, temp1_init, temp2_init);

   /* Update the boundaries */
   updateBound(x0, NeighBor, comm2d, column_type, me, xs, ys, xe, ye, xcell);

   step = 0;
   t = 0.0;

   tstart = MPI_Wtime();

   /* Main loop */
   while(!convergence)
   {
    step = step + 1;
    t = t + dt ;

    /* Perform one step of the explicit scheme */
    Explicit(x0, x, dt, &resLoc, hx, hy, me, xs, ys, xe, ye, k0);

    /* Update the partial solution along the interface */
    updateBound(x0, NeighBor, comm2d, column_type, me, xs, ys, xe, ye, xcell);

    /* Reduce all cartesian subgrids for convergence */
    MPI_Allreduce(&resLoc, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

    result= sqrt(result);

    if ((result<epsilon) || (step>maxStep)) break;
   }

   j=1;
   for(i=ys[me];i<=ye[me];i++)
   {
    for(k=0;k<=xcell-1;k++) 
      xtemp[(j-1)*xcell+k] = x0[i][xs[me]+k];
    j=j+1;
   }

   MPI_Gather(xtemp, xcell*ycell, MPI_DOUBLE, xfinal, xcell*ycell, MPI_DOUBLE, 0, comm);

   tend = MPI_Wtime();
   elapsed_time = tend - tstart;

   /* Output results */
   if(me == 0)
   {
    printf("\n");
    printf("  Time step = %3.18f\n", dt);
    printf("\n");
    printf("  Convergence = %11.9f after %d steps\n", epsilon, step);
    printf("\n");
    printf("  Problem size = %d\n", size_x*size_y);
    printf("\n");
    printf("  Wall Clock = %15.6f\n", elapsed_time);

    /* Print the solution at each point of the grid */
    printf("\n");
    printf("  Computed solution in outputPar.dat\n");
    printf("\n");

    file=fopen("outputPar.dat", "w");

    for(i=1;i<=size_x+1;i++)
      fprintf(file,"%15.11f ", temp1_init);
    fprintf(file,"%15.11f\n", temp1_init);

    for(i=1;i<=y_domains;i++)
    {
     for(j=1;j<=ycell;j++)
     {
      fprintf(file,"%15.11f ", temp1_init);
      for(k=0;k<=x_domains-1;k++)
        for(l=0;l<=xcell-1;l++)
          fprintf(file,"%15.11f ", xfinal[(j-1)*xcell+l+((y_domains-i)+k*y_domains)*xcell*ycell]);
      fprintf(file,"%15.11f\n", temp1_init);
     }
    }

    for(i=1;i<=size_x+1;i++)
      fprintf(file,"%15.11f ", temp1_init);
    fprintf(file,"%15.11f\n", temp1_init);

    fclose(file);
   }

   /* Free all arrays */
   free(x);
   free(x0);
   free(xfinal);
   free(xs);
   free(xe);
   free(ys);
   free(ye); 

   MPI_Type_free(&column_type);

   MPI_Finalize();

   return 0;
}
