GFORTRAN    = gfortran -Wall
CC          = gcc -Wall
MPI_FORTRAN = mpif90 -Wall
MPI_CC      = mpicc -Wall
LD = -lm

.SUFFIXES : .o .c

all: explicitSeq explicitPar

explicitSeq : explicitSeq.o explUtilSeq.o
	$(CC)  -o $@ explicitSeq.o explUtilSeq.o $(LD)

explicitSeq.o : explicitSeq.c 
	$(CC) -c $(*F).c
	
explUtilSeq.o : explUtilSeq.c 
	$(CC) -c $(*F).c	
		
explicitPar : explicitPar.o explUtilPar.o updateBound.o readParam.o
	$(MPI_CC)  -o $@ explicitPar.o explUtilPar.o updateBound.o readParam.o $(LD)

.c.o :
	$(MPI_CC) -c $(*F).c

clean : 
	/bin/rm -f *.o explicitSeq explicitPar
