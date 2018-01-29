GFORTRAN    = gfortran -Wall
CC          = gcc -Wall
MPI_FORTRAN = mpif90 -Wall
MPI_CC      = mpicc -Wall
LD = -lm

.SUFFIXES : .o .c

all: laplace_full laplace_reduced

		
laplace_full : laplace_full.o 
	$(MPI_CC)  -o $@ laplace_full.o  $(LD)

laplace_reduced : laplace_reduced.o 
	$(MPI_CC)  -o $@ laplace_reduced.o $(LD)
.c.o :
	$(MPI_CC) -c $(*F).c

clean : 
	/bin/rm -f *.o laplace_full laplace_reduced
