FC=mpif90  #fortran compiler
FFLAGS=-O 
SRC=src/main.f90 src/pearsonrsub.f90 src/mct.f90 src/cdfnor.f90 src/cumnor.f90 src/devlpl.f90 src/dinvnr.f90 src/gamma.f90 src/hypser.f90 src/spmpar.f90 src/stvaln.f90 src/ipmpar.f90
OBJ=${SRC:.f90=.o} #substitute .f90 with .o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $<

sim: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: #cleans all the old compilation files
	@rm -f *.mod src/*.o sim

