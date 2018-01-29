#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "explicitPar.h"
#include <unistd.h>
#define min(a,b) a <= b ? a : b

int main(int argc, char *argv[])
{
	/* Sizes for discretization */
	int size_x, size_y, size_z, me, x_domains, y_domains, z_domains;
	int size_x_full, size_y_full, size_z_full, x_domains_full, y_domains_full, z_domains_full;

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
//	int ratio=1;
        int nproc_full=64;
	int nproc_reduced = atoi(argv[1]);
	/* For reading parameters */
	int iconf[7];
	double conf[2];
	int iconf_full[7];
	double conf_full[2];

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
	int *x_index, *y_index,*z_index;
	/* MPI initialization */
	MPI_Init(&argc, &argv);
	int world_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);

	int color;
	if(world_rank>=nproc_reduced)
		color=1;
	else
		color=0;
	//   int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
	MPI_Comm_split(MPI_COMM_WORLD,  color,  world_rank, &comm);
	MPI_Comm_size(comm,&nproc);
	//printf("Reduced: size is %d, world rank id is %d\n",nproc, world_rank);

	MPI_Comm_size(comm,&nproc);
	MPI_Comm_rank(comm,&me);

	/* Get input parameters */
	if(me==0){
		readParam_reduced(iconf, conf,nproc_reduced);
		readParam(iconf_full, conf_full);

	}
	/* Broadcast input parameters */
	MPI_Bcast(iconf,7,MPI_INT,0,comm);
	MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);
	MPI_Bcast(iconf_full,7,MPI_INT,0,comm);
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
size_x_full    = iconf_full[0];
        size_y_full    = iconf_full[1];
        size_z_full    = iconf_full[2];
        x_domains_full = iconf_full[3];
        y_domains_full = iconf_full[4];
        z_domains_full = iconf_full[5];

	/* Warning message if dimensions and number of processes don't match */
	if((me==0) && (nproc!=(x_domains*y_domains*z_domains)))
		printf("Reduced: Number of processes not equal to Number of subdomains\n");

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
			printf("Reduced: Time step too large in 'param' file - Taking convergence criterion\n");
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
	x_index= malloc(nproc*sizeof(int));
	y_index= malloc(nproc*sizeof(int));
	z_index=  malloc(nproc*sizeof(int));


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
	processToMap(me, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains,x_index,y_index,z_index);

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

	MPI_Datatype matrix_reduce2full;
	int subarray[3];
	subarray[0]=size_x/x_domains_full;
	subarray[1]=size_y/y_domains_full;
	subarray[2]=size_z/z_domains_full;
	MPI_Type_create_subarray(3, sizes, subarray, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_reduce2full);
	MPI_Type_commit(&matrix_reduce2full);
	/*map the reduced modle process to full model*/
	int xyz_index[4];
	xyz_index[0]=x_index[me];
	xyz_index[1]=y_index[me];
	xyz_index[2]=z_index[me];
	xyz_index[3]=me;
	for(i=0;i<nproc_full;i++)
		MPI_Send(xyz_index,4,MPI_INT,i+nproc_reduced,1,MPI_COMM_WORLD);
	int xyz_index_recv[4];
	int index=0;
	int proc_map[8][4];
	for(i=0;i<nproc_full;i++){
		MPI_Recv(xyz_index_recv,4,MPI_INT,i+nproc_reduced,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
		if(xyz_index[0]==xyz_index_recv[0]/(x_domains_full/x_domains)&&xyz_index[1]==xyz_index_recv[1]/(y_domains_full/y_domains)&&xyz_index[2]==xyz_index_recv[2]/(z_domains_full/z_domains))
		{  proc_map[index][3]=i+nproc_reduced;
			proc_map[index][0]=xyz_index_recv[0];
			proc_map[index][1]=xyz_index_recv[1];
			proc_map[index][2]=xyz_index_recv[2];

			//printf("Reduce process# %d is mapped to %d\n",me,proc_map);
			//printf("Full:%d, %d,%d,%d; Reduced: %d,%d,%d,%d\n",me,xyz_index[0],xyz_index[1],xyz_index[2],i,xyz_index_recv[0],xyz_index_recv[1],xyz_index_recv[2]);
			index++;		
		}}
	if(index<(nproc_full/nproc_reduced))
		printf("Processor %d has less map-processor than 2",me);

	/* Initialize values */
	initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init);

	//initValues4(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init,nproc, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains);
	/* Update the boundaries */
	updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d,
			matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

	/* Initialize step and time */
	step = 0;
	t = 0.0;

	/* Starting time */
	time_init = MPI_Wtime();
	FILE *luo;
	luo=fopen("reduce_cube.txt","a");
	fprintf(luo,"%d,%d,%d,%d,%d,%d,%d\n",me, xs[me], ys[me], zs[me], x_index[me], y_index[me], z_index[me]);

	fclose(luo);
	/* if(me==0){
	   x[xs[me]+1][ys[me]][zs[me]]=8.8888;
	   MPI_Send(&x[xs[me]][ys[me]][zs[me]],1,matrix_reduce2full,8,2,MPI_COMM_WORLD);
	   }*/
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
		/*if(me==0){
		  FILE* file1;
		  file1=fopen("result.dat", "w");
		  fprintf(file1,"%15.11f ", result);
		  fclose(file1);

		  }*/
		if(me==0&&step%10000==0){
			int step_full;   
			MPI_Send(&step,1,MPI_INT,8,0,MPI_COMM_WORLD);
			// while(step_full<4*step)
			MPI_Recv(&step_full,1,MPI_INT,8,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			printf("Reduced: Received data from Rank 8, step_full=%d, step=%d\n",step_full,step);
		}
		if(step==maxStep){
			for(i=0;i<(nproc_full/nproc_reduced);i++)
				MPI_Send(&x0[xs[me]+(proc_map[i][0]%(x_domains_full/x_domains))*xcell/(x_domains_full/x_domains)][ys[me]+(proc_map[i][1]%(y_domains_full/y_domains))*ycell/(y_domains_full/y_domains)][zs[me]+(proc_map[i][2]%(z_domains_full/z_domains))*zcell/(z_domains_full/z_domains)],1,matrix_reduce2full,proc_map[i][3],2,MPI_COMM_WORLD);
			break;
		}
		/* Break conditions of main loop */
		if ((result<epsilon) || (step>maxStep)) break;
	}

	printf("Reduced: Rank #%d is finished\n",me);   
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
		printf("Reduced: Time step = %3.18f\n",dt);
		printf("\n");
		printf("  Convergence = %11.9f after %d steps\n",epsilon,step);
		printf("\n");
		printf("  Problem size = %d\n",size_x*size_y*size_z);
		printf("\n");
		printf("  Wall Clock = %15.6f\n",elapsed_time);
		printf("\n");
		printf("  Computed solution in outputPar_reduced.dat\n");
		printf("\n");

		/* Store solution into output file */
		file=fopen("outputPar_reduced.dat","w");

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


		/* file=fopen("outputPar_mesh.dat","w");
		   for(k=0;k<=size_z+1;k=k+4)

		   for(j=0;j<=size_y+1;j=j+4)

		   for(i=0;i<=size_x+1;i=i+4)
		   fprintf(file,"%15.11f ",x0[i][j][k]);

		   for(p=1;p<=size_z;p=p+4) 

		   for(i=1;i<=y_domains;i++) 
		   for(j=1;j<=ycell;j++) 
		   for(k=0;k<=x_domains-1;k++) 
		   for(l=0;l<=xcell-1;l++) 
		   fprintf(file,"%15.11f ",xfinal[(j-1)*xcell+l+(y_domains-i)*(z_domains*xcell*ycell*zcell)+k*(y_domains*z_domains*xcell*ycell*zcell)+(p-1)*xcell*ycell]);




		   fclose(file);*/



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
	free(x_index);
	free(y_index);
	free(z_index);
	/* Free matrices type */
	MPI_Type_free(&matrix_type_oxz);
	MPI_Type_free(&matrix_type_oxy);
	MPI_Type_free(&matrix_type_oyz);

	MPI_Finalize();

	return 0;
}
