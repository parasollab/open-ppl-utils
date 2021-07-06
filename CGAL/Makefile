platform ?= LINUX_gcc

# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Directories for our build.
SOURCE_DIR  := $(CURDIR)/source
BUILD_DIR   := $(CURDIR)/build
INSTALL_DIR := $(CURDIR)/install

# Important system directories.
ifeq ($(platform), MACOS_gcc)
  SYS_BIN := /usr/local/bin
  SYS_INC := /usr/local/include
  SYS_LIB := /usr/local/lib
else ($(platform), Ubuntu)
  SYS_BIN := /usr/bin
  SYS_INC := /usr/include
  SYS_LIB := /usr/lib/x86_64-linux-gnu
else
  SYS_BIN := /usr/bin
  SYS_INC := /usr/include
  SYS_LIB := /usr/lib64

  # Look for the newer version of eigen3 on UIUC machines
  UIUC_EIGEN_DIR     := /software/eigen-3.3.7/include/eigen3
  SIGNATURE_FILE     := $(UIUC_EIGEN_DIR)/signature_of_eigen3_matrix_library
  NEW_VERSION_EXISTS := $(wildcard $(SIGNATURE_FILE))
  ifneq ($(NEW_VERSION_EXISTS), "")
    EIGEN_DIR := $(UIUC_EIGEN_DIR)
  endif
endif

# Use the default eigen directory if we didn't set it above.
EIGEN_DIR ?= $(SYS_INC)/eigen3


# CGAL Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CGAL_VER := 4.14
GIT_REPO := https://github.com/CGAL/cgal.git
GIT_OPS  := --depth 1 --branch releases/CGAL-$(CGAL_VER)

BUILD_SHARED_LIBS := ON
ifeq ($(platform), MACOS_gcc)
  BUILD_SHARED_LIBS := OFF
endif

LIB_EXTENSION := so
ifeq ($(platform), MACOS_gcc)
  LIB_EXTENSION := dylib
endif

ifeq ($(platform), MACOS_gcc)
  CMAKE_BIN := cmake
else ($(platform), Ubuntu)
  CMAKE_BIN := cmake
else
  CMAKE_BIN := cmake3
endif

CMAKE_OPS := -G"Unix Makefiles" \
  -DBLAS_DEFINITIONS= \
  -DBLAS_INCLUDE_DIR=$(SYS_INC) \
  -DBLAS_LIBRARIES=$(SYS_LIB)/libblas.so \
  -DBLAS_LIBRARIES_DIR=$(SYS_LIB) \
  -DBUILD_DOC=OFF \
  -DBUILD_SHARED_LIBS=$(BUILD_SHARED_LIBS) \
  -DBUILD_TESTING=OFF \
  -DCGAL_CREATE_CMAKE_SCRIPT=$(SOURCE_DIR)/Scripts/scripts/cgal_create_cmake_script \
  -DCGAL_CTEST_DISPLAY_MEM_AND_TIME=OFF \
  -DCGAL_DEV_MODE=OFF \
  -DCGAL_CXX_FLAGS= -frounding-math -Wno-error=catch-value \
  -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS=TRUE \
  -DCGAL_ENABLE_CHECK_HEADERS=OFF \
  -DCGAL_ENABLE_PRECONFIG=ON \
  -DCGAL_HEADER_ONLY=OFF \
  -DCGAL_INSTALL_BIN_DIR=bin \
  -DCGAL_INSTALL_CMAKE_DIR=lib64/cmake \
  -DCGAL_INSTALL_DOC_DIR=doc \
  -DCGAL_INSTALL_INC_DIR=include \
  -DCGAL_INSTALL_LIB_DIR=lib64 \
  -DCGAL_INSTALL_MAN_DIR=man \
  -DCGAL_REPORT_DUPLICATE_FILES=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR) \
  -DEIGEN3_INCLUDE_DIR=$(EIGEN_DIR) \
  -DGIT_EXECUTABLE=$(SYS_BIN)/git \
  -DGMP_INCLUDE_DIR=$(SYS_INC) \
  -DGMP_LIBRARIES=$(SYS_LIB)/libgmp.$(LIB_EXTENSION) \
  -DGMP_LIBRARIES_DIR=$(SYS_LIB) \
  -DLAPACK_DEFINITIONS= \
  -DLAPACK_INCLUDE_DIR=$(SYS_INC) \
  -DLAPACK_LIBRARIES=$(SYS_LIB)/liblapack.so \
  -DLAPACK_LIBRARIES_DIR=$(SYS_LIB) \
  -DMPFR_INCLUDE_DIR=$(SYS_INC) \
  -DMPFR_LIBRARIES=$(SYS_LIB)/libmpfr.$(LIB_EXTENSION) \
  -DMPFR_LIBRARIES_DIR=$(SYS_LIB) \
  -DWITH_BLAS=ON \
  -DWITH_CGAL_Core=ON \
  -DWITH_CGAL_ImageIO=OFF \
  -DWITH_CGAL_Qt5=OFF \
  -DWITH_Coin3D=OFF \
  -DWITH_ESBTL=OFF \
  -DWITH_Eigen3=ON \
  -DWITH_GMP=ON \
  -DWITH_GMPXX=OFF \
  -DWITH_IPE=OFF \
  -DWITH_LAPACK=ON \
  -DWITH_LEDA=OFF \
  -DWITH_MPFR=ON \
  -DWITH_NTL=OFF \
  -DWITH_OpenGL=OFF \
  -DWITH_OpenNL=OFF \
  -DWITH_RS=OFF \
  -DWITH_RS3=OFF \
  -DWITH_ZLIB=OFF \
  -DWITH_demos=OFF \
  -DWITH_examples=OFF \
  -DWITH_tests=OFF

# Recipes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Build the CGAL library and generate a local 'install' here.
$(INSTALL_DIR): $(SOURCE_DIR)
	@echo Building CGAL library...
	@rm -rf $(BUILD_DIR)
	@mkdir $(BUILD_DIR)
	@cd $(BUILD_DIR) && $(CMAKE_BIN) $(CMAKE_OPS) $(SOURCE_DIR) && make install
	@rm -r $(BUILD_DIR)

# Fetch the CGAL git repo.
$(SOURCE_DIR):
	@echo Fetching CGAL git repo...
	@git clone $(GIT_REPO) $(GIT_OPS) $(SOURCE_DIR)


# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@echo Cleaning CGAL library...
	@rm -rf $(BUILD_DIR) $(INSTALL_DIR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
