FC	    = gfortran
MPIFLG  = -O3 -fcheck=all

# compile
all: auxfl andft cleanx

andft: 
	$(FC) -o wannier.x variable.o read.o mat_alloc.o mat_prep.o sdp_full.o sdp_lw.o sdp_lw_sub.o sdp_lw_sub_re.o write_wf.o r8lib.o matrix_exponential.o header.o main.f90 -fopenmp
auxfl:
	$(FC) -c variable.f90 read.f90 mat_alloc.f90 mat_prep.f90 sdp_full.f90 sdp_lw.f90 sdp_lw_sub.f90 sdp_lw_sub_re.f90 write_wf.f90 r8lib.f90 matrix_exponential.f90 header.f90 -fopenmp

clean:
	rm -rf *o *mod wannier.x

cleanx:
	rm -rf *o *mod  
