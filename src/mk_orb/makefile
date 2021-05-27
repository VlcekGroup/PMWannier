FC	    = gfortran
MPIFLG  = -O3

# compile
all: auxfl andft cleanx

andft: 
	$(FC) -o mkorb.x variable.o read.o sub.o main.f90
auxfl:
	$(FC) -c variable.f90 read.f90 sub.f90

clean:
	rm -rf *o *mod mkorb.x

cleanx:
	rm -rf *o *mod  
