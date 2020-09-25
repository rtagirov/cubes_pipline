FC=ifort 
#FC= mpif90    
PROG=testcube

#ifort preprocessor flags
FPPFLAGS = #"-DMPIF" 
FFLAGS = "-c"# -traceback  -heap-arrays -check bounds" 
#  -stand f90  -assume realloc_lhs  -check all  -traceback   -fstack-protector  -assume protect_parens"  #-O2 


# NETCDF library routines
INCLUDE="-I/scratch/witzke/Libraries/netcdf-3.6.1/include"
NETCDFLIB="-L../netcdf -lnet -L/scratch/witzke/Libraries/netcdf-3.6.1/lib -lnetcdf" 

####################################################################

SUBDIRS = netcdf  src
 
  

all: testmain 

testmain:
	 @for i in $(SUBDIRS) ; do \
                cd $$i ; \
                $(MAKE)         \
                        FC=$(FC) \
                        EXEC=$(PROG) \
                        FFLAGS=$(FFLAGS) \
                        FPPFLAGS=$(FPPFLAGS) \
                        FFTLIB=$(FFTLIB) \
                        INCLUDE=$(INCLUDE) \
                        NETCDFLIB=$(NETCDFLIB) \
                        STATIC=$(STATIC) \
                        all ; cd .. ;\
        done



clean:
	@for i in $(SUBDIRS); do \
                (cd $$i; \
                $(MAKE) clean ); \
        done
	rm -f *~

