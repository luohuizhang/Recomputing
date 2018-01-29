
/**************************************************************************/
/*                                                                        */
/* This subroutin computes next values in subdomain of current process me */
/* within the domain                                                      */
/*                                                                        */
/**************************************************************************/

void computeNext(double*** x0, double*** x, int size_x, int size_y, int size_z, double dt, double hx, double hy, 
                 double hz, double* r, int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze, double k0)
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
   for(k=zs[me];k<=ze[me];k++)
     for(j=ys[me];j<=ye[me];j++)
       for(i=xs[me];i<=xe[me];i++)
         x[i][j][k] = weightx *(x0[i-1][j][k] + x0[i+1][j][k] + x0[i][j][k]*diagx)
                    + weighty *(x0[i][j-1][k] + x0[i][j+1][k] + x0[i][j][k]*diagy)
                    + weightz *(x0[i][j][k-1] + x0[i][j][k+1] + x0[i][j][k]*diagz);

   /* Copy back the computed value : x0(n) <-- x(n) */
   /* and compute the 2_norm of the 'residual' */
   *r = 0.0;
   for(k=zs[me];k<=ze[me];k++)
     for(j=ys[me];j<=ye[me];j++)
       for(i=xs[me];i<=xe[me];i++)
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
   for(l=1;l<=nb_layers+1;l++)
   {
    for(j=(l-1);j<=y_dim-l;j++)
      for(k=(l-1);k<=z_dim-l;k++)
      {
       x0[l-1][j][k] = temp1_init;
       x0[x_dim-l][j][k] = temp1_init;
      }

    for(i=(l-1);i<=x_dim-l;i++)
      for(k=(l-1);k<=z_dim-l;k++)
      {
       x0[i][l-1][k] = temp1_init;
       x0[i][y_dim-l][k] = temp1_init;
      }

    for(i=(l-1);i<=x_dim-l;i++)
      for(j=(l-1);j<=y_dim-l;j++)
      {
       x0[i][j][l-1] = temp1_init;
       x0[i][j][z_dim-l] = temp1_init;
      }
   }

   /* Setup temp2_init inside */
   for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
       for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
         x0[i][j][k] = temp2_init;
}

/***************************************************************/
/*                                                             */
/* This subroutine computes the coordinates xs, xe, ys, ye, zs */
/* ze for each cell on the grid                                */
/*                                                             */
/***************************************************************/

void processToMap(int me, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{
   /* Index variables */
   int i, j, k, l, m, p, v;

   /* Computation of xs and xe with processes topology */
   for(i=0;i<=(z_domains*y_domains)-1;i++)
   {
    xs[i] = 2;
    xe[i] = xs[i]+xcell-1;
   }

   for(j=1;j<=x_domains-1;j++)
     for(k=0;k<=(z_domains*y_domains-1);k++)
     {
      xs[j*(z_domains*y_domains)+k] = xs[(j-1)*(z_domains*y_domains)]+xcell+2;
      xe[j*(z_domains*y_domains)+k] = xs[j*(z_domains*y_domains)]+xcell-1;
     }

   /* Computation of ys and ye with processes topology */
   for(i=1;i<=y_domains;i++) {
      ys[(i-1)*z_domains] = y_domains*(ycell+2)-ycell*i-2*(i-1);
      ye[(i-1)*z_domains] = ys[(i-1)*z_domains]+ycell-1;

      for(l=1;l<=z_domains-1;l++) {
         ys[(i-1)*z_domains+l] = ys[(i-1)*z_domains];
         ye[(i-1)*z_domains+l] = ys[(i-1)*z_domains+l]+ycell-1;
      }
   }

   /* Prolongation along y_domain */
   for(m=1;m<=y_domains;m++) {
      ys[(m-1)*z_domains] = y_domains*(ycell+2)-ycell*m-2*(m-1);
      ye[(m-1)*z_domains] = ys[(m-1)*z_domains]+ycell-1;

      for(i=1;i<=x_domains-1;i++) {
         ys[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[(m-1)*z_domains];
         ye[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[i*(y_domains*z_domains)+(m-1)*z_domains]+ycell-1;

         for(l=1;l<=z_domains-1;l++) {
            ys[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains];
            ye[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains+l]+ycell-1;
         }
      }
   }

   /* Computation of zs and ze with processes topology */
   for(k=0;k<=y_domains-1;k++)
   {
    v = k*z_domains;
    zs[v] = 2;
    ze[v] = 2+zcell-1;
    for(p=1;p<=x_domains-1;p++)
    {
     zs[v+p*(y_domains*z_domains)] = zs[v];
     ze[v+p*(y_domains*z_domains)] = ze[v];
    }
   }

   /* Prolongation along z_domain */
   for(m=1;m<=z_domains-1;m++)
      for(i=0;i<=y_domains-1;i++)
      {
       l = m+i*z_domains;
       zs[l] = zs[l-1]+zcell+2;
       ze[l] = zs[l]+zcell-1;
       for(v=1;v<=x_domains-1;v++)
       {
        zs[l+v*(y_domains*z_domains)] = zs[l];
        ze[l+v*(y_domains*z_domains)] = zs[l+v*(y_domains*z_domains)]+zcell-1;
       }
      }
}
