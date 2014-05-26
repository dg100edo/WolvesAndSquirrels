##########
# Makefile
##########

CC = gcc
DEBUGFLAGS = -g
MPICC = mpicc
MPIFLAGS = 
MPIHOSTS = localhost
MPINP = 1 #number of processors
MPIARGS =

all: serial omp mpi 

serial:
	$(CC) $(DEBUGFLAGS) -fopenmp -o wolves-squirrels-serial wolves-squirrels-serial.c $(LIBS)

omp:
	$(CC) $(DEBUGFLAGS) -fopenmp -o wolves-squirrels-omp wolves-squirrels-omp.c $(LIBS)

mpi:
	$(MPICC) $(DEBUGFLAGS) -fopenmp -o wolves-squirrels-mpi wolves-squirrels-mpi.c $(LIBS)

mpirun: mpi
	mpirun --host $(MPIHOSTS) -np $(MPINP) wolves-squirrels-mpi $(MPIARGS)
   
clean:
	rm -f *.o  wolves-squirrels-serial wolves-squirrels-omp wolves-squirrels-mpi
