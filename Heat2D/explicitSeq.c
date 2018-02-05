#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "explicitSeq.h"

#define min(a,b) a <= b ? a : b

int main(void) {

   /* Size for the discretization */
   int size_x, size_y, size_tot_x, size_tot_y;

   /* Arrays */
   double **x;
   double **x0;

   /* Spacing and time steps */
   double dt1, dt2, dt, hx, hy;

   /* Error and convergence */
   double result, epsilon;

   /* Output file descriptor */
   FILE* file;

   /* Convergence value for main loop */
   int convergence = 0;

   /* Index variables */
   int i, j;

   /* Time variable */
   double t;
   int step;

   /* Max step */
   int maxStep;

   /* Variables for clock */
   clock_t time_init, time_final;
   double elapsed_time;

   /* temp1_init: temperature init on edges */
   double temp1_init = 10.0;

   /* temp2_init: temperature init inside */
   double temp2_init = -10.0;

   /* Diffusivity coefficient */
   double k0 = 1;

   /* Getting input parameters */
   printf("Size x of the square \n");
   scanf("%d", &size_x);
   printf("Size y of the square \n");
   scanf("%d", &size_y);
   printf("Max. number of steps \n");
   scanf("%d", &maxStep);
   printf("Time step\n");
   scanf("%lf", &dt1);
   printf("Convergence \n");
   scanf("%lf", &epsilon);

   hx  = 1.0/(double)(size_x+2);
   hy  = 1.0/(double)(size_y+2);
   dt2 = 0.25*((min(hx,hy))*(min(hx,hy)))/k0;

   /* Taking a good step for convergence */
   if(dt1>=dt2)
   {
    printf("\n");
    printf("  Time step too large in param file - Taking convergence criterion\n");
    dt = dt2;
   }
   else dt = dt1;

   size_tot_x = size_x + 2;
   size_tot_y = size_y + 2;

   /* Allocation of 2D arrays */
   x = malloc(size_tot_y*sizeof(*x));
   x0 = malloc(size_tot_y*sizeof(*x0));

   for(i=0;i<=size_tot_y-1;i++)
   {
    x[i] = malloc(size_tot_x*sizeof(**x));
    x0[i] = malloc(size_tot_x*sizeof(**x0));
   }

   /* Initialization */
   initValues(x0, size_x, size_y, temp1_init, temp2_init);

   step = 0;
   t = 0.0;

   time_init=clock();

   /* Main loop */
   while(!convergence)
   {
    step = step + 1;
    t = t + dt;

    /* Perform one step of the explicit scheme */
    Explicit(x0, x, size_x, size_y, dt, hx, hy, &result, k0) ;

    result = sqrt(result);

    if ((result<epsilon) || (step>maxStep)) break;
   }

   time_final=clock();
   elapsed_time = (time_final - time_init)*1e-6;

   /* Output results */
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
   printf("  Computed solution in outputSeq.dat\n");
   printf("\n");

   file=fopen("outputSeq.dat", "w");

   for(j=size_y+1;j>=0;j--)
   {
    for(i=0;i<=size_x;i++)
      fprintf(file,"%15.11f ", x0[j][i]);
    fprintf(file,"%15.11f\n", x0[j][size_x+1]);
   }
   fclose(file);

   /* Free all arrays */
   for(i=0;i<=size_tot_y-1;i++)
   {
    free(x[i]);
    free(x0[i]);
   }
   free(x);
   free(x0);

   return 0;
}
