FC	    = gfortran
MPIFLG  = -O3

# compile
all: auxfl andft cleanx

andft: 
	$(FC) -o findWF.x variable.o read.o mat_alloc.o mat_prep.o find_wf.o write_wf.o main.f90 -fopenmp
auxfl:
	$(FC) -c variable.f90 read.f90 mat_alloc.f90 mat_prep.f90 find_wf.f90 write_wf.f90 -fopenmp

clean:
	rm -rf *o *mod findWF.x

cleanx:
	rm -rf *o *mod  
