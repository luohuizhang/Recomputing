#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "explicitSeq.h"

#define min(a,b) a <= b ? a : b

int main(void) {

   /* Sizes for the discretization */
   int size_x, size_y, size_z, size_tot_x, size_tot_y, size_tot_z;

   /* Arrays */
   double ***x;
   double ***x0;

   /* Space and time steps */
   double dt1, dt2, dt, hx, hy, hz, min1, min2;

   /* Current error and limit convergence */
   double result, epsilon;

   /* Output file descriptor */
   FILE* file;

   /* Convergence pseudo-boolean */
   int convergence = 0;

   /* Index variables */
   int i, j, k;

   /* Time and step variables */
   double t;
   int step;

   /* Max step */
   int maxStep;

   /* Variables for clock */
   clock_t time_init, time_final;
   double elapsed_time;

   /* Number of initial borders layers to */
   /* avoid artifacts problem on corners  */
   int nb_layers = 1;

   /* temp1_init: temperature init on borders */
   double temp1_init = 10.0;

   /* temp2_init: temperature init inside */
   double temp2_init = -10.0;

   /* Diffusivity coefficient */
   double k0 = 1;

   /* Get input parameters */
   printf("Size x of the square \n");
   scanf("%d",&size_x);
   printf("Size y of the square \n");
   scanf("%d",&size_y);
   printf("Size z of the square \n");
   scanf("%d",&size_z);
   printf("Max. number of steps \n");
   scanf("%d",&maxStep);
   printf("Time step\n");
   scanf("%lf",&dt1);
   printf("Convergence \n");
   scanf("%lf",&epsilon);

   /* Compute space and time steps */
   hx = 1.0/(double)(size_x+2);
   hy = 1.0/(double)(size_y+2);
   hz = 1.0/(double)(size_z+2);
   min1 = min(hx,hy);
   min2 = min(min1,hz);
   dt2 = 0.125*min2*min2*min2/k0;

   /* Take a right time step for convergence */
   if(dt1>=dt2)
   {
    printf("\n");
    printf("  Time step too large in param file - Taking convergence criterion\n");
    dt=dt2;
   }
   else dt=dt1;

   /* Allocation of 3D arrays */
   size_tot_x = size_x+2;
   size_tot_y = size_y+2;
   size_tot_z = size_z+2;

   x = malloc(size_tot_x*sizeof(*x));
   x0 = malloc(size_tot_x*sizeof(*x0));

   for(i=0;i<=size_tot_x-1;i++)
   {
    x[i] = malloc(size_tot_y*sizeof(**x));
    x0[i] = malloc((size_tot_y)*sizeof(**x0));
    for(j=0;j<=size_tot_y-1;j++)
    {
     x[i][j] = malloc(size_tot_z*sizeof(***x));
     x0[i][j] = malloc(size_tot_z*sizeof(***x0));
    }
   }

   /* Initialize values */
   initValues(nb_layers, x0, size_x, size_y, size_z, temp1_init, temp2_init);

   /* Initialize step and time */
   step = 0;
   t = 0.0;

   /* Starting time */
   time_init = clock();

   /* Main loop */
   while(!convergence)
   {
      /* Increment step and time */
      step = step + 1;
      t = t + dt;

      /* Perform one step of the explicit scheme */
      computeNext(x0, x, size_x, size_y, size_z, dt, hx, hy, hz, &result, k0);

      /* Current error */
      result = sqrt(result);

      /* Break conditions of main loop */
      if ((result<epsilon) || (step>maxStep)) break;
   }

   /* Ending time */
   time_final = clock();
   /* Elapsed time */
   elapsed_time = (time_final - time_init)*1e-6;

   /* Print results */
   printf("\n");
   printf("  Time step = %3.18f\n",dt);
   printf("\n");
   printf("  Convergence = %11.9f after %d steps\n",epsilon,step);
   printf("\n");
   printf("  Problem size = %d\n",size_x*size_y*size_z);
   printf("\n");
   printf("  Wall Clock = %15.6f\n",elapsed_time);
   printf("\n");
   printf("  Computed solution in outputSeq.dat\n");
   printf("\n");

   /* Store solution into output file */
   file=fopen("outputSeq.dat","w");

   for(k=0;k<=size_z;k++)
   {
    for(j=size_y+1;j>=0;j--)
    {
     for(i=0;i<=size_x;i++)
       fprintf(file,"%15.11f ",x0[i][j][k]);
     fprintf(file,"%15.11f\n",x0[size_x+1][j][k]);
    }
    fprintf(file,"\n");
   }

   for(j=size_y+1;j>=0;j--)
   {
    for(i=0;i<=size_x;i++)
      fprintf(file,"%15.11f ",x0[i][j][size_z+1]);
    fprintf(file,"%15.11f\n",x0[size_x+1][j][size_z+1]);
   }

   fclose(file);

   /* Free arrays */
   for(i=0;i<=size_tot_x-1;i++)
   {
    free(x[i]);
    free(x0[i]);
   }
   free(x);
   free(x0);

   return 0;
}
