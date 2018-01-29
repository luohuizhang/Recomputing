GFORTRAN    = gfortran -Wall
CC          = gcc -Wall
MPI_FORTRAN = mpif90 -Wall
MPI_CC      = mpicc -Wall
LD = -lm

.SUFFIXES : .o .c

all: explicitSeq explicitPar

explicitSeq : explicitSeq.o explUtilSeq.o
	$(CC) $(LD) -o $@ explicitSeq.o explUtilSeq.o

explicitSeq.o : explicitSeq.c 
	$(CC) -c $(*F).c
	
explUtilSeq.o : explUtilSeq.c
	$(CC) -c $(*F).c
		
explicitPar : explicitPar.o explUtilPar.o updateBound.o readParam.o
	$(MPI_CC) $(LD) -o $@ explicitPar.o explUtilPar.o updateBound.o readParam.o 

.c.o :
	$(MPI_CC) -c $(*F).c

clean : 
	/bin/rm -f *.o explicitSeq explicitPar
