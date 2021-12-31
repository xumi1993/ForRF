FC = ifort
FCFLAG = -Lsaclibs -lsac -lsacio -L/usr/local/intel/2018u3/mkl/lib/intel64
INC = /usr/local/intel/2018u3/compilers_and_libraries_2018.3.222/linux/mkl/include


all: mkl_fft interpolation_mod make_rf


mkl_fft:
	$(FC) -I$(INC) $(INC)/mkl_dfti.f90  mklfft.f90 -c

interpolation_mod:
	$(FC) interpolation_mod.f90 -c

make_rf:
	$(FC) deconit.f90 interpolation_mod.o drwsac.f90 mklfft.o utils.f90 -o make_rf $(FCFLAG) -mkl

clean:
	rm make_rf *.mod *.o