
# M A K E F I L E
#
# ARCHICTECTURE = linux 
#

BIN       = ../bin
PROGRAMS  = $(BIN)/pot
AUX       = aux
FC        = gfortran 	#without this line, gcc assigns FC=f77
FCOMPLR   = gfortran -O3 
FCFLAGS   = -c  -w -lgfortran # -C=all -gline #-info #note -mismatch_all is needed to resolve MPI_BCAST, etc...
#FCFLAGS   = -c -free  -mismatch_all -w  # -C=all -gline #-info #note -mismatch_all is needed to resolve MPI_BCAST, etc...
FINCLPATH = 
LDFLAGS   = 
OBJECTS   = main.o geometry.o globals.o in_out.o energy.o utils.o $(AUX)/lbfgsb.o ls_rmsd.o 
LIB_OBJECTS = geometry.o in_out.o globals.o energy.o utils.o $(AUX)/lbfgsb.o ls_rmsd.o
AR = ar cr

default: $(PROGRAMS)

$(BIN)/pot : $(OBJECTS) makefile
	$(FCOMPLR) $(LDFLAGS) $(OBJECTS) -o $@

library: $(LIB_OBJECTS) makefile
	g++ -o libEEF1.so -g -O2 -fexceptions -shared $(LIB_OBJECTS)

staticlibrary: $(LIB_OBJECTS) makefile
	$(AR) libEEF1.a $(LIB_OBJECTS)

.c.o: 
	$(CCOMPLR) $(CCFLAGS)  $(CINCPATH) $<

.F90.o: 
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) $<

#  specific dependencies

main.o: main.F90 energy.o utils.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) main.F90

utils.o: utils.F90 energy.o ls_rmsd.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) utils.F90

in_out.o: in_out.F90 globals.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) in_out.F90

globals.o: globals.F90
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) globals.F90

geometry.o: geometry.F90 in_out.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) geometry.F90

energy.o: energy.F90 geometry.o $(AUX)/lbfgsb.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) energy.F90

ls_rmsd.o: ls_rmsd.F90 globals.o
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) ls_rmsd.F90

$(AUX)/lbfgsb.o: $(AUX)/lbfgsb.f90
	$(FCOMPLR) $(FCFLAGS) $(F90INCPATH) $(AUX)/lbfgsb.f90 -o $(AUX)/lbfgsb.o

clean: 
	mv quickfit.o quickfit.keep
	rm -f *.o *~ *.mod *.g90 $(AUX)/*.o $(AUX)/*.g90 $(PROGRAMS) libEEF1.so libEEF1.a
	mv quickfit.keep quickfit.o
