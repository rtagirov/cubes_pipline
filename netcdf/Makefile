OBJ = IOnet.o

all: libnet.a

libnet.a: $(OBJ)
	rm -f libnet.a
	ar crv libnet.a $(OBJ)

IOnet.o : IOnet.F90
	$(FC) $(FFLAGS) $(FPPFLAG)  IOnet.F90 $(INCLUDE)

clean:
	rm -f *.o core *.a *.f *~ work.pc *.mod 
 
