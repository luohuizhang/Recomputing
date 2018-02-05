void Explicit( double** x0, double** x, int size_x, int size_y, double dt, double
      hx, double hy, double* r, double k0)
{
   /* Index variables */
   int i, j;

   /* Factors for the stencil */
   double diagx, diagy, weightx, weighty, rk;

   /*
      The stencil of the explicit operator for the heat equation
      on a regular rectangular grid using a five point finite difference
      scheme in space is :

      |                                    weightx * x[i-1][j]                                    |
      |                                                                                           |
      | weighty * x[i][j-1]   (diagx * weightx + diagy * weighty) * x[i][j]   weightx * x[i][j+1] |
      |                                                                                           |
      |                                    weighty * x[i+1][j]                                    |
   */

   diagx = -2.0 + hx*hx/(2*k0*dt);
   diagy = -2.0 + hy*hy/(2*k0*dt);
   weightx = k0*dt/(hx*hx);
   weighty = k0*dt/(hy*hy);

   /* Perform an explicit update on the points within the domain */
   for(j=1;j<=size_x;j++)
     for(i=1;i<=size_y;i++)
       x[i][j] = weighty*(x0[i-1][j] + x0[i+1][j] + x0[i][j]*diagy)
               + weightx*(x0[i][j-1] + x0[i][j+1] + x0[i][j]*diagx);

   /* Copy back the computed value : x  <-- x^(n+1)             */
   /*                                x0 <-- x^n                 */
   /* and compute at the same time the 2_norm of the 'residual' */
   *r = 0.0;
   for(j=1;j<=size_x;j++)
     for(i=1;i<=size_y;i++)
     {
      rk = x0[i][j] - x[i][j];
      *r  = *r + rk*rk;
      x0[i][j] = x[i][j];
     }
}

/**************************************************************************/
/*                                                                        */
/* This subroutine setups the initial guess, i.e. the initial temperature */
/* within the domain                                                      */
/*                                                                        */
/**************************************************************************/

void initValues(double** x0, int x_dim, int y_dim, double temp1_init, double temp2_init )
{
   /* Index variables */
   int i, j;

   /* Setup temp1_init on edges */
   for(j=0;j<=x_dim+1;j++)
   {
    x0[0][j] = temp1_init;
    x0[y_dim+1][j] = temp1_init;
   }

   for(i=0;i<=y_dim+1;i++)
   {
    x0[i][0] = temp1_init;
    x0[i][x_dim+1] = temp1_init;
   }

   /* Setup temp2_init inside */
   for(j=1;j<=y_dim;j++)
     for(i=1;i<=x_dim;i++)
       x0[j][i] = temp2_init;
}
