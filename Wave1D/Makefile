GFORTRAN    = gfortran -Wall
CC          = gcc -Wall
MPI_FORTRAN = mpif90 -Wall
MPI_CC      = mpicc -Wall
LD = -lm

.SUFFIXES : .o .c

all:  wave_full wave_reduced

		
wave_full : wave_full.o 
	$(MPI_CC)  -o $@ wave_full.o  $(LD)

wave_reduced : wave_reduced.o 
	$(MPI_CC)  -o $@ wave_reduced.o $(LD)
.c.o :
	$(MPI_CC) -c $(*F).c

clean : 
	/bin/rm -f *.o wave_full wave_reduced
