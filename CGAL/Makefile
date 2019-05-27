# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Directories for our build.
SOURCE_DIR  := $(CURDIR)/source
BUILD_DIR   := $(CURDIR)/build
INSTALL_DIR := $(CURDIR)/install

# Important system directories.
SYS_BIN := /usr/bin
SYS_INC := /usr/include
SYS_LIB := /usr/lib64


# CGAL Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Version 4.9.1 is the highest version that works on cmake v2.x, which is
# required for CentOS. Check for ability to upgrade on next major CentOS
# version.
CGAL_VER := 4.9.1
GIT_REPO := https://github.com/CGAL/cgal.git
GIT_OPS  := --depth 1 --branch releases/CGAL-$(CGAL_VER)

CMAKE_OPS := -G"Unix Makefiles" \
  -DBLAS_DEFINITIONS= \
  -DBLAS_INCLUDE_DIR= \
  -DBLAS_LIBRARIES=$(SYS_LIB)/libblas.so \
  -DBLAS_LIBRARIES_DIR=$(SYS_LIB) \
  -DBUILD_DOC=OFF \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_TESTING=OFF \
  -DCGAL_CREATE_CMAKE_SCRIPT=$(SOURCE_DIR)/Scripts/scripts/cgal_create_cmake_script \
  -DCGAL_CXX_FLAGS= -frounding-math \
  -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS=TRUE \
  -DCGAL_ENABLE_CHECK_HEADERS=OFF \
  -DCGAL_ENABLE_PRECONFIG=ON \
  -DCGAL_INSTALL_BIN_DIR=bin \
  -DCGAL_INSTALL_CMAKE_DIR=lib64/CGAL \
  -DCGAL_INSTALL_DOC_DIR=doc \
  -DCGAL_INSTALL_INC_DIR=include \
  -DCGAL_INSTALL_LIB_DIR=lib64 \
  -DCGAL_INSTALL_MAN_DIR=man \
  -DCGAL_REPORT_DUPLICATE_FILES=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR) \
  -DEIGEN3_INCLUDE_DIR=$(SYS_INC)/eigen3 \
  -DGIT_EXECUTABLE=$(SYS_BIN)/git \
  -DGMP_INCLUDE_DIR=$(SYS_INC) \
  -DGMP_LIBRARIES=$(SYS_LIB)/libgmp.so \
  -DGMP_LIBRARIES_DIR=$(SYS_LIB) \
  -DLAPACK_DEFINITIONS= \
  -DLAPACK_INCLUDE_DIR=$(SYS_INC) \
  -DLAPACK_LIBRARIES=$(SYS_LIB)/liblapack.so \
  -DLAPACK_LIBRARIES_DIR=$(SYS_LIB) \
  -DMPFR_INCLUDE_DIR=$(SYS_INC) \
  -DMPFR_LIBRARIES=$(SYS_LIB)/libmpfr.so \
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
  -DWITH_MPFI=OFF \
  -DWITH_MPFR=ON \
  -DWITH_NTL=OFF \
  -DWITH_OpenGL=OFF \
  -DWITH_OpenNL=OFF \
  -DWITH_QGLViewer=OFF \
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
	@cd $(BUILD_DIR) && cmake $(CMAKE_OPS) $(SOURCE_DIR) && make install
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
