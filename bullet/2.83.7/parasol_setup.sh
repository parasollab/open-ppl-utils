#!/usr/bin/sh
#
# Setup the bullet library for use with PMPL.

# Configure directory layout.
base_dir=$(pwd)
build_dir=${base_dir}/build
install_dir=${build_dir}/install

# Remove existing build and install directories if they exist.
rm -rf ${build_dir} ${install_dir}

# Create the build and install directories, and generate the makefile.
mkdir ${build_dir} ${install_dir}
cd ${build_dir}
cmake -G"Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${install_dir} \
  -DINSTALL_EXTRA_LIBS=TRUE \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O2 -g3 -DNDEBUG" \
  ../

# Bullet takes a while to build, so use parallel make with max number of
# processors.
numprocs=$(grep -c ^processor /proc/cpuinfo)
make install -j${numprocs}
