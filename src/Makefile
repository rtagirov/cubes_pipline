FC=mpif90  
#FC = "mpif90"
PROG=test

FFLAGS = -c -traceback -heap-arrays -check bounds#-O2 

INCLUDE=-I/opt/local/cfitsio/cfitsio-3.350/include
FITSLIB=-L/opt/local/cfitsio/cfitsio-3.350/lib -lcfitsio

INCLUDENET=-I/scratch/witzke/Libraries/netcdf-3.6.1/include
NETCDFLIB=-L../netcdf -lnet -L/scratch/witzke/Libraries/netcdf-3.6.1/lib -lnetcdf


OBJECTS=arrays.o io_cubes.o  initcalc.o rotate.o

####################################################################
## --with-cfitsio=/opt/local/cfitsio/cfitsio-3.350#


all: main 

main:    $(OBJECTS)  
	 $(FC) $(LFLAGS) -o test main.F90  $(OBJECTS) $(INCLUDE) $(FITSLIB) $(INCLUDENET)  $(NETCDFLIB)   
	mv test ../bin 

arrays.o : arrays.f90
	$(FC) $(FFLAGS) arrays.f90 

io_cubes.o : io_cubes.f90
	$(FC) $(FFLAGS) io_cubes.f90 $(INCLUDE) $(FITSLIB) 

rotate.o : rotate.f90
	$(FC) $(FFLAGS) rotate.f90

initcalc.o : initcalc.f90
	$(FC) $(FFLAGS) initcalc.f90
 
# ========================================================================
clean:
	rm -f *~ *.o *.mod main~  
