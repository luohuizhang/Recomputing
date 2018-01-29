
/*********************************************************/
/*                                                       */
/* This subroutine computes next values on global domain */
/*                                                       */
/*********************************************************/

void computeNext(double*** x0, double*** x, int size_x, int size_y, int size_z, double dt,
                 double hx, double hy, double hz, double* r, double k0)
{
   /* Index variables */
   int i, j, k;

   /* Factors for the stencil */
   double diagx, diagy, diagz, weightx, weighty, weightz;
   
   /* Local variable for computing error */
   double rk;

   /*
   The stencil of the explicit operator for the heat equation
   on a regular rectangular grid using a seven point finite difference
   scheme in space is :

   |                                       wx * x[i-1][j][k]     wz * x[i][j][k-1]               |
   |                                                                  /                          |
   | wy * x[i][j-1][k]   (diagx * wx + diagy * wy + diagz * wz) * x[i][j][k]   wx * x[i][j+1][k] |
   |                            /                                                                |
   |                 wz * x[i][j][k+1]     wy * x[i+1][j][k]                                     |
   */

   diagx = -2.0 + hx*hx/(3*k0*dt);
   diagy = -2.0 + hy*hy/(3*k0*dt);
   diagz = -2.0 + hz*hz/(3*k0*dt);
   weightx = k0*dt/(hx*hx);
   weighty = k0*dt/(hy*hy);
   weightz = k0*dt/(hz*hz);

   /* Perform an explicit update on the points within the domain */
   for(k=1;k<=size_z;k++)
     for(j=1;j<=size_y;j++)
       for(i=1;i<=size_x;i++)
         x[i][j][k] = weightx*(x0[i-1][j][k] + x0[i+1][j][k] + x0[i][j][k]*diagx)
                    + weighty*(x0[i][j-1][k] + x0[i][j+1][k] + x0[i][j][k]*diagy)
                    + weightz*(x0[i][j][k-1] + x0[i][j][k+1] + x0[i][j][k]*diagz);

   /* Copy back the computed value : x0(n) <-- x(n) */
   /* and compute the 2_norm of the 'residual' */
   *r = 0.0;
   for(k=1;k<=size_z;k++)
     for(j=1;j<=size_y;j++)
       for(i=1;i<=size_x;i++)
       {
        rk = x0[i][j][k] - x[i][j][k];
        *r  = *r + rk*rk;
        x0[i][j][k] = x[i][j][k];
       }
}    

/**************************************************************************/
/*                                                                        */
/* This subroutine setups the initial guess, i.e. the initial temperature */
/* within the domain                                                      */
/*                                                                        */
/**************************************************************************/

void initValues(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init)
{
   /* Index variables */
   int i, j, k, l;

   /* Setup temp1_init on edges */
   for(l=1;l<=nb_layers;l++)
   {
    for(j=(l-1);j<=y_dim+1-(l-1);j++)
      for(k=(l-1);k<=z_dim+1-(l-1);k++)
      {
       x0[l-1][j][k] = temp1_init;
       x0[x_dim+1-(l-1)][j][k] = temp1_init;
      }

    for(i=(l-1);i<=x_dim+1-(l-1);i++)
      for(k=(l-1);k<=z_dim+1-(l-1);k++)
      {
       x0[i][l-1][k] = temp1_init;
       x0[i][y_dim+1-(l-1)][k] = temp1_init;
      }

    for(j=(l-1);j<=y_dim+1-(l-1);j++)
      for(i=(l-1);i<=x_dim+1-(l-1);i++)
         {
          x0[i][j][l-1] = temp1_init;
          x0[i][j][z_dim+1-(l-1)] = temp1_init;
         }
   }

   /* Setup temp2_init inside */
   for(i=nb_layers;i<=x_dim-(nb_layers-1);i++)
     for(j=nb_layers;j<=y_dim-(nb_layers-1);j++)
       for(k=nb_layers;k<=z_dim-(nb_layers-1);k++)
         x0[i][j][k] = temp2_init;
}
