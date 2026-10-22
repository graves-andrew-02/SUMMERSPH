compile-mp:
		gfortran SUMMER_SPH-Variable-Implicit_energy.f90 -o run_sph.exe -fopenmp


compile:
		gfortran SUMMER_SPH-Variable-Implicit_energy.f90 -o run_sph.exe

comp-run-mp:
		gfortran SUMMER_SPH-Variable-Implicit_energy.f90 -o run_sph.exe -fopenmp
		./run_sph.exe
