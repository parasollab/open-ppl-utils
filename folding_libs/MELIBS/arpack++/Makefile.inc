# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# This file contains some definitions used to compile arpack++ examples
# with the g++ compiler under linux.


# Defining the machine.

PLAT         = linux

# Defining the compiler.

CPP          = g++

# Defining ARPACK++ directories.
# ARPACKPP_INC is the directory that contains all arpack++ header files.
# SUPERLU_DIR and UMFPACK_DIR must be set to ARPACKPP_INC.

ARPACKPP_DIR = /research/blackwidowers/dsmft2/RNAFolding/arpack++
ARPACKPP_INC = $(ARPACKPP_DIR)/include
SUPERLU_DIR  = $(ARPACKPP_INC)
UMFPACK_DIR  = $(ARPACKPP_INC)

# Defining ARPACK, LAPACK, UMFPACK, SUPERLU, BLAS and FORTRAN libraries.
# See the arpack++ manual or the README file for directions on how to 
# obtain arpack, umfpack and SuperLU packages. 
# UMFPACK_LIB and SUPERLU_LIB must be declared only if umfpack and superlu 
# are going to be used. Some BLAS and LAPACK fortran routines are 
# distributed along with arpack fortran code, but the user should verify 
# if optimized versions of these libraries are available before installing 
# arpack. The fortran libraries described below are those required to link
# fortran and c++ code using gnu g++ and f77 compiler under linux.
# Other libraries should be defined if the user intends to compile
# arpack++ on another environment.

ARPACK_LIB   = /research/blackwidowers/dsmft2/RNAFolding/ARPACK/libarpack_SUN4.a
LAPACK_LIB   = /usr/local/lib/lapack_$(PLAT).a
UMFPACK_LIB  = /usr/local/lib/umfpack_$(PLAT).a
SUPERLU_LIB  = /research/blackwidowers/dsmft2/RNAFolding/SuperLU/superlu$(PLAT).a
BLAS_LIB     = /research/blackwidowers/dsmft2/RNAFolding/SuperLU/blas$(PLAT).a 
FORTRAN_LIBS = /usr/lib/gcc-lib/i386-redhat-linux/2.96/libg2c.a

# Defining g++ flags and directories.

CPP_WARNINGS = -Wall
CPP_DEBUG    = -g
CPP_OPTIM    = -O
CPP_LIBS     = 
CPP_INC      = /usr/local/include/g++-3

CPP_FLAGS    = $(CPP_DEBUG) -D$(PLAT) -I$(ARPACKPP_INC) -I$(CPP_INC)

# Putting all libraries together.

#ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) $(LAPACK_LIB) \
#               $(BLAS_LIB) $(FORTRAN_LIBS) 

ALL_LIBS     = $(CPP_LIBS) $(ARPACK_LIB) $(BLAS_LIB) $(FORTRAN_LIBS) 

# defining paths.

vpath %.h  $(ARPACK_INC)

