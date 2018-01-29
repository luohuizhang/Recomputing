
/**************************************************************************/
/*                                                                        */
/* This subroutin computes next values in subdomain of current process me */
/* within the domain                                                      */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include<math.h>

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
  if(me==0){
// printf("dt=%lf,diagx=%f,weightx=%f,diagy=%f,weighty=%f\n",dt,diagx,weightx,diagy,weighty);
}
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
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
        x0[i][j][k]=temp2_init;

}




/***************************************************************/
/*                                                             */
/* This subroutine computes the coordinates xs, xe, ys, ye, zs */
/* ze for each cell on the grid                                */
/*                                                             */
/***************************************************************/

void processToMap(int me, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains,int *x_index,int *y_index, int *z_index)
{
   /* Index variables */
   int i, j, k, l, m, p, v;

   /* Computation of xs and xe with processes topology */
   for(i=0;i<=(z_domains*y_domains)-1;i++)
   {
    xs[i] = 2;
    xe[i] = xs[i]+xcell-1;
    x_index[i]=0;
   }

   for(j=1;j<=x_domains-1;j++)
     for(k=0;k<=(z_domains*y_domains-1);k++)
     {
      xs[j*(z_domains*y_domains)+k] = xs[(j-1)*(z_domains*y_domains)]+xcell+2;
      xe[j*(z_domains*y_domains)+k] = xs[j*(z_domains*y_domains)]+xcell-1;
      x_index[j*(z_domains*y_domains)+k]=j;
     }

   /* Computation of ys and ye with processes topology */
   for(i=1;i<=y_domains;i++) {
      ys[(i-1)*z_domains] = y_domains*(ycell+2)-ycell*i-2*(i-1);
      ye[(i-1)*z_domains] = ys[(i-1)*z_domains]+ycell-1;
      y_index[(i-1)*z_domains]=y_domains-i;
      for(l=1;l<=z_domains-1;l++) {
         ys[(i-1)*z_domains+l] = ys[(i-1)*z_domains];
         ye[(i-1)*z_domains+l] = ys[(i-1)*z_domains+l]+ycell-1;
        y_index[(i-1)*z_domains+l]=y_domains-i; 
      }
   }

   /* Prolongation along y_domain */
   for(m=1;m<=y_domains;m++) {
      ys[(m-1)*z_domains] = y_domains*(ycell+2)-ycell*m-2*(m-1);
      ye[(m-1)*z_domains] = ys[(m-1)*z_domains]+ycell-1;
      y_index[(m-1)*z_domains]=y_domains-m;

      for(i=1;i<=x_domains-1;i++) {
         ys[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[(m-1)*z_domains];
         ye[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[i*(y_domains*z_domains)+(m-1)*z_domains]+ycell-1;
         y_index[i*(y_domains*z_domains)+(m-1)*z_domains]=y_domains-m;
         for(l=1;l<=z_domains-1;l++) {
            ys[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains];
            ye[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains+l]+ycell-1;
            y_index[i*(y_domains*z_domains)+(m-1)*z_domains+l]=y_domains-m;
         }
      }
   }

   /* Computation of zs and ze with processes topology */
   for(k=0;k<=y_domains-1;k++)
   {
    v = k*z_domains;
    zs[v] = 2;
    ze[v] = 2+zcell-1;
    z_index[v]=0;
    for(p=1;p<=x_domains-1;p++)
    {
     zs[v+p*(y_domains*z_domains)] = zs[v];
     ze[v+p*(y_domains*z_domains)] = ze[v];
     z_index[v+p*(y_domains*z_domains)]=0;
    }
   }

   /* Prolongation along z_domain */
   for(m=1;m<=z_domains-1;m++)
      for(i=0;i<=y_domains-1;i++)
      {
       l = m+i*z_domains;
       zs[l] = zs[l-1]+zcell+2;
       ze[l] = zs[l]+zcell-1;
       z_index[l]=m;
       for(v=1;v<=x_domains-1;v++)
       {
        zs[l+v*(y_domains*z_domains)] = zs[l];
        ze[l+v*(y_domains*z_domains)] = zs[l+v*(y_domains*z_domains)]+zcell-1;
        z_index[l+v*(y_domains*z_domains)]=m;
       }
      }
}





void initValues2(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init, int nproc, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{
   /* Index variables */
   int i, j, k, l,v,p;

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
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
        x0[i][j][k]=temp2_init;

 double temp;

 FILE* file;	
 file=fopen("luo.dat","r");
/*top flame*/
	


for(j=1;j<=ycell*y_domains+2;j++)
   	for(i=1;i<=xcell*x_domains+2;i++)
		fscanf(file,"%lf",&temp);
	
		
		
	
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++){
     for(v=1;v<=xcell*x_domains+2;v++)
	fscanf(file,"%lf",&temp);

     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++){
	fscanf(file,"%lf",&temp);
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++){
        for(p=0;p<=nproc;p++)
	if(xs[p]<=i&&i<=xe[p]&&ys[p]<=j&&j<=ye[p]&&zs[p]<=k&&k<=ze[p])
   	fscanf(file,"%lf",&x0[i][j][k]);
	if(p==nproc)
	 printf("can not assign x0[%d][%d][%d]!\n",i,j,k);

        }
	fscanf(file,"%lf",&temp);
	}
     for(v=1;v<=xcell*x_domains+2;v++)
	fscanf(file,"%lf",&temp);
   }
/*bottom flame*/
for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
		fscanf(file,"%lf",&temp);
	


fclose(file);


   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
	if(x0[i][j][k]==temp2_init){
	if(x0[i-1][j][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i-1][j][k])/2;
	else if (x0[i][j-1][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j-1][k])/2;
	
	else if(x0[i][j][k-1]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j][k-1])/2;
	else
		printf("temp2_init is set! x0[%d][%d][%d]!\n",i,j,k);

   /* Setup temp1_init inside
    int luo=(int)((x_dim-1)/2);
    for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
       for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
         x0[luo][j][k] = temp1_init; */
}

}

void initValues3(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init, int nproc, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{
   /* Index variables */
   int i, j, k, l,v,p;

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
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
        x0[i][j][k]=temp2_init;

 double temp;
 FILE* file;	
 file=fopen("input_mesh.dat_10000","r");
/*top flame*/
	

for(j=1;j<=ycell*y_domains/4+1;j++)
   	for(i=1;i<=xcell*x_domains/4+1;i++){
		fscanf(file,"%lf",&temp);

		}
		
	
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k+=4){
	for(v=0;v<=nproc;v++)
	if(zs[v]<=k&&k<=ze[v]){
      	for(i=1;i<=xcell*x_domains/4+1;i++){
		fscanf(file,"%lf",&temp);
	
		}
	break;}
	
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j+=4){
	
	for(v=0;v<=nproc;v++)
	if(ys[v]<=j&&j<=ye[v]&&zs[v]<=k&&k<=ze[v]){
      	
		fscanf(file,"%lf",&temp);
		
		
		break;}

   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i+=4){
        for(p=0;p<=nproc;p++)
	if(xs[p]<=i&&i<=xe[p]&&ys[p]<=j&&j<=ye[p]&&zs[p]<=k&&k<=ze[p]){
	int xindex=(i-xs[p])/4*4+xs[p]+3;
	int yindex=(j-ys[p])/4*4+ys[p]+3;
	int zindex=(k-zs[p])/4*4+zs[p]+3;
   	fscanf(file,"%lf",&x0[xindex][yindex][zindex]);

	//printf("xindex,xs[p],xe[p]:%d,%d,%d\n",xindex,xs[p],xe[p]);	
	//printf("yindex,ys[p],ye[p]:%d,%d,%d\n",yindex,ys[p],ye[p]);
	//printf("zindex,zs[p],ze[p]:%d,%d,%d\n",zindex,zs[p],ze[p]);
	}//if
	if(p==nproc)
	 printf("can not assign x0[%d][%d][%d]!\n",i,j,k);
	
        }//x for
	
	}//y for
	
   }//z for
	
//printf("before:%d\n",ftell(file));
//fseek(file,0,SEEK_END);
//printf("after:%d\n",ftell(file));
fclose(file);


   /* Setup temp1_init inside
    int luo=(int)((x_dim-1)/2);
    for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
       for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
         x0[luo][j][k] = temp1_init; */

   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
	if(x0[i][j][k]==temp2_init){	
	for(p=0;p<=nproc;p++)
	if(xs[p]<=i&&i<=xe[p]&&ys[p]<=j&&j<=ye[p]&&zs[p]<=k&&k<=ze[p]){
	int xindex=(i-xs[p])/4*4+xs[p]+3;
	int yindex=(j-ys[p])/4*4+ys[p]+3;
	int zindex=(k-zs[p])/4*4+zs[p]+3;
   	x0[i][j][k]=x0[xindex][yindex][zindex];

	//printf("xindex,xs[p],xe[p]:%d,%d,%d\n",xindex,xs[p],xe[p]);	
	//printf("yindex,ys[p],ye[p]:%d,%d,%d\n",yindex,ys[p],ye[p]);
	//printf("zindex,zs[p],ze[p]:%d,%d,%d\n",zindex,zs[p],ze[p]);
	}//if
	if(p==nproc){
	if(x0[i-1][j][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i-1][j][k])/2;
	else if (x0[i][j-1][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j-1][k])/2;
	
	else if(x0[i][j][k-1]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j][k-1])/2;
	
        }
	}		






}

void initValues4(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init, int nproc, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{
   /* Index variables */
   int i, j, k, l,v,p;

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
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
        x0[i][j][k]=temp2_init;

 double temp;
 FILE* file;	
 file=fopen("input_mesh.dat_10000","r");
/*top flame*/
	

for(j=1;j<=ycell*y_domains/4+1;j++)
   	for(i=1;i<=xcell*x_domains/4+1;i++){
		fscanf(file,"%lf",&temp);

		}
		
	
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k+=4){
	for(v=0;v<=nproc;v++)
	if(zs[v]<=k&&k<=ze[v]){
      	for(i=1;i<=xcell*x_domains/4+1;i++){
		fscanf(file,"%lf",&temp);
	
		}
	break;}
	
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j+=4){
	
	for(v=0;v<=nproc;v++)
	if(ys[v]<=j&&j<=ye[v]&&zs[v]<=k&&k<=ze[v]){
      	
		fscanf(file,"%lf",&temp);
		
		
		break;}

   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i+=4){
        for(p=0;p<=nproc;p++)
	if(xs[p]<=i&&i<=xe[p]&&ys[p]<=j&&j<=ye[p]&&zs[p]<=k&&k<=ze[p]){
	int xindex=(i-xs[p])/4*4+xs[p]+3;
	int yindex=(j-ys[p])/4*4+ys[p]+3;
	int zindex=(k-zs[p])/4*4+zs[p]+3;
   	fscanf(file,"%lf",&x0[xindex][yindex][zindex]);

	//printf("xindex,xs[p],xe[p]:%d,%d,%d\n",xindex,xs[p],xe[p]);	
	//printf("yindex,ys[p],ye[p]:%d,%d,%d\n",yindex,ys[p],ye[p]);
	//printf("zindex,zs[p],ze[p]:%d,%d,%d\n",zindex,zs[p],ze[p]);
	}//if
	if(p==nproc)
	 printf("can not assign x0[%d][%d][%d]!\n",i,j,k);
	
        }//x for
	
	}//y for
	
   }//z for
	
//printf("before:%d\n",ftell(file));
//fseek(file,0,SEEK_END);
//printf("after:%d\n",ftell(file));
fclose(file);


   /* Setup temp1_init inside
    int luo=(int)((x_dim-1)/2);
    for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
       for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
         x0[luo][j][k] = temp1_init; */




   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++){
		
	for(p=0;p<=nproc;p++)
	if(xs[p]<=i&&i<=xe[p]&&ys[p]<=j&&j<=ye[p]&&zs[p]<=k&&k<=ze[p]){
	int x_high=(i-xs[p])/4*4+xs[p]+3;
	int y_high=(j-ys[p])/4*4+ys[p]+3;
	int z_high=(k-zs[p])/4*4+zs[p]+3;
	int x_low, y_low,z_low;

	
	if(i-xs[p]>=4)
		x_low=x_high-4;
	else if(i>5)
		x_low=x_high-6;
	else
	    	x_low=1;
	if(j-ys[p]>=4)
		y_low=y_high-4;
	else if(j>5)
		y_low=y_high-6;
	else
	    	y_low=1;

	if(k-zs[p]>=4)
		z_low=z_high-4;
	else if(k>5)
		z_low=z_high-6;
	else
	    	z_low=1;

	double ii=4+i-x_high;
	double jj=4+j-y_high;
	double kk=4+k-z_high;
	
	
   	double t7=x0[x_high][y_high][z_high];
	double t8=x0[x_high][y_low][z_high];	
	double t6=x0[x_low][y_high][z_high];
	double t5=x0[x_low][y_low][z_high];
	double t4=x0[x_high][y_low][z_low];
	double t3=x0[x_high][y_high][z_low];
	double t2=x0[x_low][y_high][z_low];
	double t1=x0[x_low][y_low][z_low];

	double temp1=ii/4*(jj/4*t3+(1-jj/4)*t4)+(1-ii/4)*(jj/4*t2+(1-jj/4)*t1);
	double temp2=ii/4*(jj/4*t7+(1-jj/4)*t8)+(1-ii/4)*(jj/4*t6+(1-jj/4)*t5);
	
	x0[i][j][k]=kk/4*temp2+(1-kk/4)*temp1;

	
	

	//printf("xindex,xs[p],xe[p]:%d,%d,%d\n",xindex,xs[p],xe[p]);	
	//printf("yindex,ys[p],ye[p]:%d,%d,%d\n",yindex,ys[p],ye[p]);
	//printf("zindex,zs[p],ze[p]:%d,%d,%d\n",zindex,zs[p],ze[p]);
	}//if
	
	
        if(p==nproc){
	if(x0[i-1][j][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i-1][j][k])/2;
	else if (x0[i][j-1][k]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j-1][k])/2;
	
	else if(x0[i][j][k-1]!=temp2_init)
		x0[i][j][k]=(x0[i][j][k]+x0[i][j][k-1])/2;
	
        }
	}





}


void initValues5(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init)
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


	double sigma=1;
	double PI=3.14159265359;
	int center=x_dim/2;
	
	for(int i=0;i<x_dim;i++)
	{
		for(int j=0;j<y_dim;j++)
		{
			x0[i][j][1]=(1/(2*PI*sigma*sigma))*exp(-((i-center)*(i-center)+(j-center)*(j-center))/(2*sigma*sigma));
			x0[i][j][0]=(1/(2*PI*sigma*sigma))*exp(-((i-center)*(i-center)+(j-center)*(j-center))/(2*sigma*sigma));
		}
	}




   /* Setup temp2_init inside */
   for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
   	for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
        x0[i][j][k]=temp2_init;

}


void regen_linear(double*** x, double*** reduced_buffer, int ratio_x, int ratio_y, int ratio_z, int me, int* xs, int* xe, int* ys, int* ye, int* zs, int* ze, int xcell, int ycell, int zcell)
{
   /* Index variables */
   int i, j, k;


   for(k=0;k<zcell/ratio_z;k++)
     for(j=0;j<ycell/ratio_y;j++)
   	for(i=0;i<xcell/ratio_x;i++){
          x[xs[me]+ratio_x*i][ys[me]+j*ratio_y][zs[me]+k*ratio_z]=reduced_buffer[i][j][k];

	}
//  printf("rank is %d\n",me);	
   for(k=zs[me];k<=ze[me];k++)
     for(j=ys[me];j<=ye[me];j++)
   	for(i=xs[me];i<=xe[me];i++){
 
	
	int x_low=(i-xs[me])/ratio_x*ratio_x+xs[me];
	int y_low=(j-ys[me])/ratio_y*ratio_y+ys[me];
	int z_low=(k-zs[me])/ratio_z*ratio_z+zs[me];
        int x_high=x_low+ratio_x-1;
        int y_high=y_low+ratio_y-1;
        int z_high=z_low+ratio_z-1;
	double ii=ratio_x+i-x_high;
	double jj=ratio_y+j-y_high;
	double kk=ratio_z+k-z_high;
	
	
   	double t7=x[x_high][y_high][z_high];
	double t8=x[x_high][y_low][z_high];	
	double t6=x[x_low][y_high][z_high];
	double t5=x[x_low][y_low][z_high];
	double t4=x[x_high][y_low][z_low];
	double t3=x[x_high][y_high][z_low];
	double t2=x[x_low][y_high][z_low];
	double t1=x[x_low][y_low][z_low];


	double temp1=ii/ratio_x*(jj/ratio_y*t3+(1-jj/ratio_y)*t4)+(1-ii/ratio_x)*(jj/ratio_y*t2+(1-jj/ratio_y)*t1);
	double temp2=ii/ratio_x*(jj/ratio_y*t7+(1-jj/ratio_y)*t8)+(1-ii/ratio_x)*(jj/ratio_y*t6+(1-jj/ratio_y)*t5);
	
	x[i][j][k]=kk/ratio_z*temp2+(1-kk/ratio_z)*temp1;


}    

}
