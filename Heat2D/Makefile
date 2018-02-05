GFORTRAN    = gfortran -Wall
CC          = gcc -Wall
MPI_FORTRAN = mpif90 -Wall
MPI_CC      = mpicc -Wall
LD = -lm

.SUFFIXES : .o .c

all:  explicitPar explicitPar_reduced


	
		
explicitPar : explicitPar.o explUtilPar.o updateBound.o readParam.o
	$(MPI_CC)  -o $@ explicitPar.o explUtilPar.o updateBound.o readParam.o $(LD)
explicitPar_reduced : explicitPar_reduced.o explUtilPar.o updateBound.o readParam.o
	$(MPI_CC)  -o $@ explicitPar.o explUtilPar.o updateBound.o readParam.o $(LD)

.c.o :
	$(MPI_CC) -c $(*F).c

clean : 
	/bin/rm -f *.o  explicitPar explicitPar_reduced
