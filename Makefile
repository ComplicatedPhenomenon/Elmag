###############################################################################
# makefile for ELMAG 2.03                                          
#
mpi: modules203.f90 user203.f90 init203.f90 elmag203.f90 aux202.f90
	mpif90 -JModules -O -w modules203.f90 user203.f90 init203.f90 elmag203.f90 aux202.f90
#
single: modules203.f90 user_sp203.f90 init203.f90 elmag203.f90 aux202.f90
	gfortran -JModules -O modules203.f90 user_sp203.f90 init203.f90 elmag203.f90 aux202.f90
#
test: modules203.f90 user203.f90 init203.f90 elmag203.f90 aux202.f90
	mpif90 -JModules  -Wall -C -g  -fbacktrace -ffpe-trap={underflow,overflow,invalid,denormal} modules203.f90 user203.f90 init203.f90 elmag203.f90 aux202.f90
#
###############################################################################
