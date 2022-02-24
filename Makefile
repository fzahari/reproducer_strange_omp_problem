EXE=test

FC0=xlf_r
FC1=xlf90_r
FC2=xlf2003_r

FFs=-c  -qsmp=omp -g 
FFs0=-c -g 

LF= -qsmp=omp -g 

all: mod_obj obj test
fixed: mod_obj1 mod_gs mod_obj2 obj test

mod_obj: params.f90 mx_limits.f90 prec.f90 mod_nosp_basis.f90 mod_dft_partfunc.f90 mod_grid_storage.f90 mod_dft_molgrid.f90 mod_dft_fuzzycell.f90 test.f90
	$(FC2) $(FFs) -c $^

mod_obj1: mx_limits.f90 prec.f90 
	$(FC2) $(FFs) -c $^

mod_gs: mod_grid_storage.f90
	$(FC2) $(FFs0) -c $^

mod_obj2: params.f90 mod_nosp_basis.f90 mod_dft_partfunc.f90 mod_dft_molgrid.f90 mod_dft_fuzzycell.f90 test.f90
	$(FC2) $(FFs) -c $^

obj: svpleb.f radpt.f
	$(FC0) $(FFs) -c $^

test: obj
	$(FC2) -o test $(LF) *.o

.PHONY: clean
clean:
	rm *.o *.mod $(EXE)
