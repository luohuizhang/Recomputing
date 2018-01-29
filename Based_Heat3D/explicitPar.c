#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "explicitPar.h"

#define min(a,b) a <= b ? a : b

int main(int argc, char *argv[])
{
   /* Sizes for discretization */
   int size_x, size_y, size_z, me, x_domains, y_domains, z_domains;
   int size_x_glo, size_y_glo, size_z_glo;

   /* Arrays */
   double ***x;
   double ***x0;
   double *x_all;
   double *x0_all;
   double *x_alloc;
   double *x0_alloc;
   double *xfinal;
   double *xtemp;

   /* For reading parameters */
   int iconf[7];
   double conf[2];

   /* Spacing and time steps */
   double dt, dt1, dt2, hx, hy, hz, min1, min2;

   /* Current error and limit convergence */
   double resLoc, result, epsilon;

   /* Output file descriptor */
   FILE* file;

   /* Convergence pseudo-boolean */
   int convergence = 0;

   /* Index variables */
   int i, j, k, l, p, v, m;

   /* Time and step variables */
   double t;
   int step;

   /* Max step */
   int maxStep;

   /* Variables for clock */
   double time_init, time_final;
   double elapsed_time;

   /* Number of initial borders layers to */
   /* avoid artifacts problem on corners  */
   int nb_layers = 1;

   /* temp1_init: temperature init on edges */
   double temp1_init = 10.0;

   /* temp2_init: temperature init inside */
   double temp2_init = -10.0;

   /* Diffusivity coefficient */
   double k0 = 1;

   /* MPI variables */
   int sizes[3], subsizes1[3], subsizes2[3], subsizes3[3], starts[3];
   int nproc, ndims;
   MPI_Comm comm, comm3d;
   int dims[3];
   int periods[3];
   int reorganisation = 0;
   MPI_Datatype matrix_type_oxz, matrix_type_oxy, matrix_type_oyz;
   int S=0, E=1, N=2, W=3, Zd=4, Zu=5;
   int NeighBor[6];
   int xcell, ycell, zcell, size_tot_x, size_tot_y, size_tot_z;
   int *xs, *ys, *zs, *xe, *ye, *ze;

   /* MPI initialization */
   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm,&nproc);
   MPI_Comm_rank(comm,&me);

   /* Get input parameters */
   if(me==0)
    readParam(iconf, conf);

   /* Broadcast input parameters */
   MPI_Bcast(iconf,7,MPI_INT,0,comm);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);

   /* Assign input parameters to variables */
   size_x    = iconf[0];
   size_y    = iconf[1];
   size_z    = iconf[2];
   x_domains = iconf[3];
   y_domains = iconf[4];
   z_domains = iconf[5];
   maxStep   = iconf[6];
   dt1       = conf[0];
   epsilon   = conf[1];

   /* Warning message if dimensions and number of processes don't match */
   if((me==0) && (nproc!=(x_domains*y_domains*z_domains)))
    printf("Number of processes not equal to Number of subdomains\n");

   /* Various other variables */
   size_x_glo = size_x+2;
   size_y_glo = size_y+2;
   size_z_glo = size_z+2;
   hx = 1.0/(double)(size_x_glo);
   hy = 1.0/(double)(size_y_glo);
   hz = 1.0/(double)(size_z_glo);
   min1 = min(hx,hy);
   min2 = min(min1,hz);
   dt2  = 0.125*min2*min2*min2/k0;
   size_tot_x = size_x+2*x_domains+2;
   size_tot_y = size_y+2*y_domains+2;
   size_tot_z = size_z+2*z_domains+2;

   /* Take a right time step for convergence */
   if(dt1>=dt2)
   {
    if(me==0)
    {
     printf("\n");
     printf("  Time step too large in 'param' file - Taking convergence criterion\n");
    }
    dt = dt2;
   }
   else dt = dt1;

   /* Allocate 3D Contiguous arrays */
   xfinal = malloc(size_x*size_y*size_z*sizeof(*xfinal));
   x_all = malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x_all)); 
   x0_all = malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x0_all));
   x_alloc = x_all;
   x0_alloc = x0_all;

   /* Allocate size_y rows */
   x = malloc(size_tot_x*sizeof(*x));
   x0 = malloc(size_tot_x*sizeof(*x0));

   /* Loop on rows */
   for (i=0;i<size_tot_x;i++)
   {
    /* Allocate size_tot_y columns for each row */
    x[i] = malloc(size_tot_y*sizeof(**x));
    x0[i] = malloc(size_tot_y*sizeof(**x0));
    /* Loop on columns */
    for(j=0;j<size_tot_y;j++)
    {
     /* Increment size_z block on x0[i][j] address */
     x[i][j] = x_alloc;
     x0[i][j] = x0_alloc;
     x_alloc += size_tot_z;
     x0_alloc += size_tot_z;
    }
   }

   /* Allocate coordinates of processes */
   xs = malloc(nproc*sizeof(int));
   xe = malloc(nproc*sizeof(int));
   ys = malloc(nproc*sizeof(int));
   ye = malloc(nproc*sizeof(int));
   zs = malloc(nproc*sizeof(int));
   ze = malloc(nproc*sizeof(int));

   /* Create 3D cartesian grid */
   periods[0] = 0;
   periods[1] = 0;
   periods[2] = 0;

   ndims = 3;
   dims[0] = x_domains;
   dims[1] = y_domains;
   dims[2] = z_domains;

   MPI_Cart_create(comm, ndims, dims, periods, reorganisation, &comm3d);

   /* Identify neighbors */
   NeighBor[0] = MPI_PROC_NULL;
   NeighBor[1] = MPI_PROC_NULL;
   NeighBor[2] = MPI_PROC_NULL;
   NeighBor[3] = MPI_PROC_NULL;
   NeighBor[4] = MPI_PROC_NULL;
   NeighBor[5] = MPI_PROC_NULL;

   /* Left/West and right/Est neigbors */
   MPI_Cart_shift(comm3d, 0, 1, &NeighBor[W], &NeighBor[E]);

   /* Bottom/South and Upper/North neigbors */
   MPI_Cart_shift(comm3d, 1, 1, &NeighBor[S], &NeighBor[N]);

   /* Zdown/South and Zup/North neigbors */
   MPI_Cart_shift(comm3d, 2, 1, &NeighBor[Zd], &NeighBor[Zu]);

   /* Size of each cell */
   xcell = (size_x/x_domains);
   ycell = (size_y/y_domains);
   zcell = (size_z/z_domains);

   /* Allocate subdomain */
   xtemp = malloc(xcell*ycell*zcell*sizeof(*xtemp));

   /* Compute xs, xe, ys, ye, zs, ze for each cell on the grid */
   processToMap(me, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains);

   /* Create matrix data types to communicate */
   sizes[0] = size_tot_x;
   sizes[1] = size_tot_y;
   sizes[2] = size_tot_z;

   starts[0] = 0;
   starts[1] = 0;
   starts[2] = 0;

   /* Create matrix data type to communicate on vertical Oxz plane */
   subsizes1[0] = xcell;
   subsizes1[1] = 1;
   subsizes1[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes1, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxz);
   MPI_Type_commit(&matrix_type_oxz);

   /* Create matrix data type to communicate on vertical Oyz plane */
   subsizes2[0] = 1;
   subsizes2[1] = ycell;
   subsizes2[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes2, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oyz);
   MPI_Type_commit(&matrix_type_oyz);

   /* Create matrix data type to communicate on vertical Oxy plane */
   subsizes3[0] = xcell;
   subsizes3[1] = ycell;
   subsizes3[2] = 1;

   MPI_Type_create_subarray(3, sizes, subsizes3, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxy);
   MPI_Type_commit(&matrix_type_oxy);

   /* Initialize values */
   initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init);

   /* Update the boundaries */
   updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d,
               matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

   /* Initialize step and time */
   step = 0;
   t = 0.0;

   /* Starting time */
   time_init = MPI_Wtime();

   /* Main loop */
   while(!convergence)
   {  
      /* Increment step and time */
      step = step + 1;
      t = t + dt ;

      /* Perform one step of the explicit scheme */
      computeNext(x0, x, size_tot_x, size_tot_y, size_tot_z, dt, hx, hy, hz, &resLoc, me, xs, ys, zs, xe, ye, ze, k0);

      /* Update the partial solution along the interface */
      updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d, 
                  matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

      /* Sum reduction to get error */
      MPI_Allreduce(&resLoc, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

      /* Current error */
      result = sqrt(result);

      /* Break conditions of main loop */
      if ((result<epsilon) || (step>maxStep)) break;
   }

   /* Gather all subdomains */
   i = 1;
   for(k=zs[me];k<=ze[me];k++)
   {
    l = 1;
    for(j=ys[me];j<=ye[me];j++)
    {
     for(m=0;m<=xcell-1;m++)
       xtemp[(l-1)*xcell+(i-1)*xcell*ycell+m] = x0[xs[me]+m][j][k];
     l = l+1;
    }
    i = i+1;
   }
   /* Perform gathering */
   MPI_Gather(xtemp, xcell*ycell*zcell, MPI_DOUBLE, xfinal, xcell*ycell*zcell, MPI_DOUBLE, 0, comm);

   /* Ending time */
   time_final = MPI_Wtime();
   /* Elapsed time */
   elapsed_time = time_final - time_init;

   /* Print results */
   if(me == 0)
   {
    printf("\n");
    printf("  Time step = %3.18f\n",dt);
    printf("\n");
    printf("  Convergence = %11.9f after %d steps\n",epsilon,step);
    printf("\n");
    printf("  Problem size = %d\n",size_x*size_y*size_z);
    printf("\n");
    printf("  Wall Clock = %15.6f\n",elapsed_time);
    printf("\n");
    printf("  Computed solution in outputPar.dat\n");
    printf("\n");

    /* Store solution into output file */
    file=fopen("outputPar.dat","w");

    for(j=1;j<=size_y+2;j++) {
      for(i=1;i<=size_x+1;i++) {
        fprintf(file,"%15.11f ",temp1_init);
      }
      fprintf(file,"%15.11f\n",temp1_init);
    }

    fprintf(file,"\n");

    for(p=1;p<=size_z;p++) {
      for(v=1;v<=size_x+1;v++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n",temp1_init);
      for(i=1;i<=y_domains;i++) {
        for(j=1;j<=ycell;j++) {
          fprintf(file,"%15.11f ",temp1_init);
          for(k=0;k<=x_domains-1;k++) {
            for(l=0;l<=xcell-1;l++) {
              fprintf(file,"%15.11f ",xfinal[(j-1)*xcell+l+(y_domains-i)*(z_domains*xcell*ycell*zcell)+k*(y_domains*z_domains*xcell*ycell*zcell)+(p-1)*xcell*ycell]);
            }
          }
          fprintf(file,"%15.11f\n",temp1_init);
        }
      }
      for(m=1;m<=size_x+1;m++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n\n",temp1_init);
    }

    for(j=1;j<=size_y+2;j++) {
      for(i=1;i<=size_x+1;i++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n",temp1_init);
    }

    fclose(file);
   }

   /* Free arrays */
   for (i=0;i<=size_tot_x-1;i++)
   {
    free(x[i]);
    free(x0[i]);
   }
   
   free(x);
   free(x0);
   free(x_all);
   free(x0_all);
   free(xfinal);
   free(xtemp);
   free(xs);
   free(xe);
   free(ys);
   free(ye);
   free(zs);
   free(ze);

   /* Free matrices type */
   MPI_Type_free(&matrix_type_oxz);
   MPI_Type_free(&matrix_type_oxy);
   MPI_Type_free(&matrix_type_oyz);

   MPI_Finalize();

   return 0;
}
