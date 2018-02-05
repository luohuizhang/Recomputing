# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
double *update ( int id, int p, int n_global, int n_local, int nsteps, 
		double dt, MPI_Comm comm );
void collect ( int id, int p, int n_global, int n_local, int nsteps, 
		double dt, double u_local[], MPI_Comm comm );
double dudt ( double x, double t );
double exact ( double x, double t );
void timestamp ( );
void regen_linear(  int id, int p, int n_global, int n_local, int ratio, double u_local[], double u_regen[]);
/******************************************************************************/

int main ( int argc, char *argv[] )

	/******************************************************************************/
	/*
	   Purpose:

	   WAVE_MPI solves the wave equation in parallel using MPI.

	   Discussion:

	   Discretize the equation for u(x,t):
	   d^2 u/dt^2 - c^2 * d^2 u/dx^2 = 0  for 0 < x < 1, 0 < t
	   with boundary conditions:
	   u(0,t) = u0(t) = sin ( 2 * pi * ( 0 - c * t ) )
	   u(1,t) = u1(t) = sin ( 2 * pi * ( 1 - c * t ) )
	   and initial conditions:
	   u(x,0) = g(x,t=0) =                sin ( 2 * pi * ( x - c * t ) )
	   dudt(x,0) = h(x,t=0) = - 2 * pi * c * cos ( 2 * pi * ( x - c * t ) ) 

	   by:

	   alpha = c * dt / dx.

	   U(x,t+dt) = 2 U(x,t) - U(x,t-dt) 
	   + alpha^2 ( U(x-dx,t) - 2 U(x,t) + U(x+dx,t) ).

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   17 November 2013

	   Author:

	   John Burkardt

	   Reference:

	   Geoffrey Fox, Mark Johnson, Gregory Lyzenga, Steve Otto, John Salmon, 
	   David Walker,
	   Solving problems on concurrent processors, 
	   Volume 1: General Techniques and Regular Problems,
	   Prentice Hall, 1988,
	   ISBN: 0-13-8230226,
	   LC: QA76.5.F627.

	   Local parameters:

	   Local, double DT, the time step.

	   Local, int ID, the MPI process ID.

	   Local, int N_GLOBAL, the total number of points.

	   Local, int N_LOCAL, the number of points visible to this process.

	   Local, int NSTEPS, the number of time steps.

	   Local, int P, the number of MPI processes.
	 */
{
	double dt = 0.0000125;
	int i_global_hi;
	int i_global_lo;
	int id;
	int n_global = 40961;

	int n_local;
	int nsteps = 6400;
	int p;
	double *u1_local;
	double wtime;
	int n_global_reduced=atoi(argv[2]);

	int nproc_reduced=atoi(argv[1]);
	int nproc_full=64;
	/* 
	   Initialize MPI.
	 */
	MPI_Init ( &argc, &argv );
	int world_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
	MPI_Comm comm;

	int color;
	if(world_rank>=nproc_reduced)
		color=1;
	else
		color=0;
	//   int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
	MPI_Comm_split(MPI_COMM_WORLD,  color,  world_rank, &comm);


	MPI_Comm_rank ( comm, &id );

	MPI_Comm_size ( comm, &p );

	if ( id == 0 ) 
	{
		timestamp ( );
		printf ( "\n" );
		printf ( "Full: MPI_WAVE:\n" );
		printf ( "Full:  Using %d processes.\n", p );
		printf ( "Full:  Using a total of %d points.\n", n_global );
		printf ( "Full:  Using %d time steps of size %g.\n", nsteps, dt );
		printf ( "Full:  Computing final solution at time %g\n", dt * nsteps );
	}

	wtime = MPI_Wtime ( );
	/*
	   Determine N_LOCAL.
	 */
	i_global_lo = (   id       * ( n_global - 1 ) ) / p;
	i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
	if ( 0 < id )
	{
		i_global_lo = i_global_lo - 1;
	}
	n_local = i_global_hi + 1 - i_global_lo;
	/* 
	   Update N_LOCAL values.
	 */
	u1_local = update ( id, p, n_global, n_local, nsteps, dt , comm);
	/* 
	   Collect local values into global array.
	 */
	int ratio=(n_global_reduced-1)/nproc_full;
	double *u_regen=( double * ) malloc ( (ratio+1) * sizeof ( double ) );
MPI_Recv ( &u_regen[0],   ratio+1, MPI_DOUBLE, 0, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	//printf("luoluo");
	regen_linear(  id,  p,  n_global,  n_local, ratio, u1_local,  u_regen);


	collect ( id, p, n_global, n_local, nsteps, dt, u1_local , comm);
	/*
	   Report elapsed wallclock time.
	 */
	wtime = MPI_Wtime ( ) - wtime;
	if ( id == 0 )
	{
		printf ( "\n" );
		printf ( "  Elapsed wallclock time was %g seconds\n", wtime );
	}
	/*
	   Terminate MPI.
	 */
	MPI_Finalize ( );
	/*
	   Free memory.
	 */
	free ( u1_local );
	/*
	   Terminate.
	 */
	if ( id == 0 )
	{
		printf ( "Full: Normal end of execution.\n" );
		timestamp ( );
	}

	return 0;
}
/******************************************************************************/

double *update ( int id, int p, int n_global, int n_local, int nsteps, 
		double dt,  MPI_Comm comm ) 

	/******************************************************************************/
	/*
	   Purpose:

	   UPDATE advances the solution a given number of time steps.

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   17 November 2013

	   Author:

	   John Burkardt

	   Parameters:

	   Input, int ID, the identifier of this process.

	   Input, int P, the number of processes.

	   Input, int N_GLOBAL, the total number of points.

	   Input, int N_LOCAL, the number of points visible to this process.

	   Input, int NSTEPS, the number of time steps.

	   Input, double DT, the size of the time step.

	   Output, double UPDATE[N_LOCAL], the portion of the solution
	   at the last time, as evaluated by this process.
	 */
{
	double alpha;
	double c;
	double dx;
	int i;
	int i_global;
	int i_global_hi;
	int i_global_lo;
	int i_local;
	int i_local_hi;
	int i_local_lo;
	int j;
	int ltor = 20;
	int rtol = 10;
	MPI_Status status;
	double t;
	double *u0_local;
	double *u1_local;
	double *u2_local;
	double x;
	/*
	   Determine the value of ALPHA.
	 */
	c = 1.0;
	dx = 1.0 / ( double ) ( n_global - 1 );
	alpha = c * dt / dx;

	if ( 1.0 <= fabs ( alpha ) )
	{
		if ( id == 0 )
		{
			fprintf ( stderr, "\n" );
			fprintf ( stderr, "UPDATE - Warning!\n" );
			fprintf ( stderr, "  1 <= |ALPHA| = | C * dT / dX |.\n" );
			fprintf ( stderr, "  C = %g\n", c );
			fprintf ( stderr, "  dT = %g\n", dt );
			fprintf ( stderr, "  dX = %g\n", dx );
			fprintf ( stderr, "  ALPHA = %g\n", alpha );
			fprintf ( stderr, "  Computation will not be stable!\n" );
		}
		MPI_Finalize ( );
		exit ( 1 );
	}
	/*
	   The global array of N_GLOBAL points must be divided up among the processes.
	   Each process stores about 1/P of the total + 2 extra slots.
	 */
	i_global_lo = (   id       * ( n_global - 1 ) ) / p;
	i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
	if ( 0 < id )
	{
		i_global_lo = i_global_lo - 1;
	}

	i_local_lo = 0;
	i_local_hi = i_global_hi - i_global_lo;

	u0_local = ( double * ) malloc ( n_local * sizeof ( double ) );
	u1_local = ( double * ) malloc ( n_local * sizeof ( double ) );
	u2_local = ( double * ) malloc ( n_local * sizeof ( double ) );

	t = 0.0;
	for ( i_global = i_global_lo; i_global <= i_global_hi; i_global++ ) 
	{
		x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
		i_local = i_global - i_global_lo;
		u1_local[i_local] = exact ( x, t );
	}

	for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
	{
		u0_local[i_local] = u1_local[i_local];
	}
	/* 
	   Take NSTEPS time steps.
	 */
	for ( i = 1; i <= nsteps; i++ )
	{
		t = dt * ( double ) i;
		/* 
		   For the first time step, we need to use the initial derivative information.
		 */
		if ( i == 1 )
		{
			for ( i_local = i_local_lo + 1; i_local < i_local_hi; i_local++ ) 
			{
				i_global = i_global_lo + i_local;
				x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
				u2_local[i_local] = 
					+         0.5 * alpha * alpha   * u1_local[i_local-1]
					+ ( 1.0 -       alpha * alpha ) * u1_local[i_local] 
					+         0.5 * alpha * alpha   * u1_local[i_local+1]
					+                            dt * dudt ( x, t );
			}
		}
		/* 
		   After the first time step, we can use the previous two solution estimates.
		 */
		else
		{
			for ( i_local = i_local_lo + 1; i_local < i_local_hi; i_local++ ) 
			{
				u2_local[i_local] = 
					+               alpha * alpha   * u1_local[i_local-1]
					+ 2.0 * ( 1.0 - alpha * alpha ) * u1_local[i_local] 
					+               alpha * alpha   * u1_local[i_local+1]
					-                                 u0_local[i_local];
			}
		}
		/* 
		   Exchange data with "left-hand" neighbor. 
		 */
		if ( 0 < id ) 
		{
			MPI_Send ( &u2_local[i_local_lo+1], 1, MPI_DOUBLE, id - 1, rtol, comm);
			MPI_Recv ( &u2_local[i_local_lo],   1, MPI_DOUBLE, id - 1, ltor, comm, &status );
		}
		else
		{
			x = 0.0;
			u2_local[i_local_lo] = exact ( x, t );
		}
		/* 
		   Exchange data with "right-hand" neighbor.
		 */
		if ( id < p - 1 ) 
		{
			MPI_Send ( &u2_local[i_local_hi-1], 1, MPI_DOUBLE, id + 1, ltor, 
					comm );
			MPI_Recv ( &u2_local[i_local_hi],   1, MPI_DOUBLE, id + 1, rtol, 
					comm, &status );
		}
		else
		{
			x = 1.0;
			u2_local[i_local_hi] = exact ( x, t );
		}
		/*
		   Shift data for next time step.
		 */
		for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
		{
			u0_local[i_local] = u1_local[i_local];
			u1_local[i_local] = u2_local[i_local];
		}
	}
	/*
	   Free memory.
	 */
	free ( u0_local );
	free ( u2_local );

	return u1_local;
}
/******************************************************************************/

void collect ( int id, int p, int n_global, int n_local, int nsteps, 
		double dt, double u_local[],  MPI_Comm comm ) 

	/******************************************************************************/
	/*
	   Purpose:

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   16 November 2013

	   Author:

	   John Burkardt

	   Parameters:

	   Input, int ID, the identifier of this process.

	   Input, int P, the number of processes.

	   Input, int N_GLOBAL, the total number of points.

	   Input, int N_LOCAL, the number of points visible to this process.

	   Input, int NSTEPS, the number of time steps.

	   Input, double DT, the size of the time step.

	   Input, double U_LOCAL[N_LOCAL], the final solution estimate computed
	   by this process.
	 */
{
	int buffer[2];
	int collect1 = 10;
	int collect2 = 20;
	int i;
	int i_global;
	int i_global_hi;
	int i_global_lo;
	int i_local;
	int i_local_hi;
	int i_local_lo;
	int j;
	int n_local2;
	MPI_Status status;
	double t;
	double *u_global;
	double x;

	i_global_lo = (   id       * ( n_global - 1 ) ) / p;
	i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
	if ( 0 < id )
	{
		i_global_lo = i_global_lo - 1;
	}

	i_local_lo = 0;
	i_local_hi = i_global_hi - i_global_lo;
	/* 
	   Master collects worker results into the U_GLOBAL array.
	 */
	if ( id == 0 )
	{
		/*
		   Create the global array.
		 */
		u_global = ( double * ) malloc ( n_global * sizeof ( double ) );
		/*
		   Copy the master's results into the global array.
		 */
		for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
		{
			i_global = i_global_lo + i_local - i_local_lo;
			u_global[i_global] = u_local[i_local];
		}
		/*
		   Contact each worker.
		 */
		for ( i = 1; i < p; i++ ) 
		{
			/*
			   Message "collect1" contains the global index and number of values.
			 */
			MPI_Recv ( buffer, 2, MPI_INT, i, collect1, comm, &status );
			i_global_lo = buffer[0];
			n_local2 = buffer[1];

			if ( i_global_lo < 0 )
			{
				fprintf ( stderr, "  Illegal I_GLOBAL_LO = %d\n", i_global_lo );
				exit ( 1 );
			}
			else if ( n_global <= i_global_lo + n_local2 - 1 )
			{
				fprintf ( stderr, "  Illegal I_GLOBAL_LO + N_LOCAL2 = %d\n", 
						i_global_lo + n_local2 );
				exit ( 1 );
			}
			/*
			   Message "collect2" contains the values.
			 */
			MPI_Recv ( &u_global[i_global_lo], n_local2, MPI_DOUBLE, i, collect2, 
					comm, &status );
		}
		/*
		   Print the results.
		 */
		FILE * full_file, *delta;
		full_file=fopen("full_data.dat","wb");
		delta=fopen("delta_data.dat","wb");
		t = dt * ( double ) nsteps;
		printf ( "\n" );
		printf ( "    I      X     F(X)   Exact\n" );
		printf ( "\n" );
		for ( i_global = 0; i_global < n_global; i_global++) 
		{
			x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
			printf ( "  %5d\t %8.3f\t%f\t%f\n", 
					i_global, x, u_global[i_global], exact ( x, t ) );
			fprintf(full_file,"%f",u_global[i_global]);
			fprintf(delta,"%f",exact ( x, t )-u_global[i_global]);
		}
		fclose(full_file);
		fclose(delta);
		free ( u_global );
	}
	/*
	   Workers send results to process 0.
	 */
	else
	{
		/*
		   Message "collect1" contains the global index and number of values.
		 */
		buffer[0] = i_global_lo;
		buffer[1] = n_local;
		MPI_Send ( buffer, 2, MPI_INT, 0, collect1, comm );
		/*
		   Message "collect2" contains the values.
		 */
		MPI_Send ( u_local, n_local, MPI_DOUBLE, 0, collect2, comm );
	}

	return;
}
/******************************************************************************/

double exact ( double x, double t )

	/******************************************************************************/
	/*
	   Purpose:

	   EXACT evaluates the exact solution

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   17 November 2013

	   Author:

	   John Burkardt

	   Parameters:

	   Input, double X, the location.

	   Input, double T, the time.

	   Output, double EXACT, the value of the exact solution.
	 */
{
	const double c = 1.0;
	const double pi = 3.141592653589793;
	double value;

	value = sin ( 2.0 * pi * ( x - c * t ) );

	return value;
}
/******************************************************************************/

double dudt ( double x, double t )

	/******************************************************************************/
	/*
	   Purpose:

	   DUDT evaluates the partial derivative dudt.

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   17 November 2013

	   Author:

	   John Burkardt

	   Parameters:

	   Input, double X, the location.

	   Input, double T, the time.

	   Output, double DUDT, the value of the time derivative of the solution.
	 */
{
	const double c = 1.0;
	const double pi = 3.141592653589793;
	double value;

	value = - 2.0 * pi * c * cos ( 2.0 * pi * ( x - c * t ) );

	return value;
}
/******************************************************************************/

void timestamp ( )

	/******************************************************************************/
	/*
	   Purpose:

	   TIMESTAMP prints the current YMDHMS date as a time stamp.

	   Example:

	   31 May 2001 09:45:54 AM

	   Licensing:

	   This code is distributed under the GNU LGPL license. 

	   Modified:

	   24 September 2003

	   Author:

	   John Burkardt

	   Parameters:

	   None
	 */
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	time_t now;

	now = time ( NULL );
	tm = localtime ( &now );

	strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	printf ( "%s\n", time_buffer );

	return;
# undef TIME_SIZE
}

void regen_linear(  int id, int p, int n_global, int n_local, int length, double u_local[], double u_regen[])
{
	int i;
	int ratio=(n_global-1)/p/length;

	for(i=0;i<=length;i++){
		if(id>0)
			u_local[1+i*ratio]=u_regen[i];
		else
			u_local[i*ratio]=u_regen[i];

	}
	if(id>0){
		u_local[0]=u_local[1];
		for(i=1;i<n_local-1;i++){


			int x_low=(i-1)/ratio*ratio+1;

			int x_high=x_low+ratio;

			u_local[i]=u_local[x_low]+(double)(i-x_low)/ratio*(u_local[x_high]-u_local[x_low]);

		}

	}
	else{
		for(i=1;i<n_local-1;i++){


			int x_low=i/ratio*ratio;

			int x_high=x_low+ratio;

			u_local[i]=u_local[x_low]+(double)(i-x_low)/ratio*(u_local[x_high]-u_local[x_low]);

		}

	}
}
