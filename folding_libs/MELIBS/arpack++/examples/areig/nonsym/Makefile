# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/areig/nonsym directory makefile.

# including other makefiles.

include ../../../Makefile.inc

CPP = g++

# define OBPRM dir (modification by B Kirkpatrick 8 Jan 2003)
OBPRM_DIR = /research/blackwidowers/dsmft2/RNAFolding/obprm/OBPRM

# defining areig directory.

AREIG_DIR = $(ARPACKPP_DIR)/examples/areig

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/nonsym

# compiling and linking all examples.

all: simple ansymreg ansymshf ansymgre ansymgsh ansymgsc

# compiling and linking each nonsymmetric problem.

simple:	simple.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -I$(OBPRM_DIR)  -o simple simple.o $(OBPRM_DIR)/Clock_Class.o $(SUPERLU_LIB) $(ALL_LIBS)

simple_original:	simple_original.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o simple_original simple_original.o $(SUPERLU_LIB) $(ALL_LIBS)
simple-vec:	simple-vec.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o simple-vec simple-vec.o $(SUPERLU_LIB) $(ALL_LIBS)

ansymreg:        ansymreg.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o ansymreg ansymreg.o $(SUPERLU_LIB) $(ALL_LIBS)

ansymshf:        ansymshf.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o ansymshf ansymshf.o $(SUPERLU_LIB) $(ALL_LIBS)

ansymgre:        ansymgre.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o ansymgre ansymgre.o $(SUPERLU_LIB) $(ALL_LIBS)

ansymgsh:        ansymgsh.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o ansymgsh ansymgsh.o $(SUPERLU_LIB) $(ALL_LIBS)

ansymgsc:        ansymgsc.o
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -o ansymgsc ansymgsc.o $(SUPERLU_LIB) $(ALL_LIBS)

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) -I$(OBPRM_DIR) -c $<


