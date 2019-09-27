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
cmake -G"Unix Makefiles" \
  -DBUILD_BULLET2_DEMOS=OFF \
  -DBUILD_BULLET3=ON \
  -DBUILD_CLSOCKET=OFF \
  -DBUILD_CPU_DEMOS=OFF \
  -DBUILD_ENET=OFF \
  -DBUILD_EXTRAS=ON \
  -DBUILD_OPENGL3_DEMOS=OFF \
  -DBUILD_PYBULLET=OFF \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_UNIT_TESTS=OFF \
  -DBULLET2_USE_THREAD_LOCKS=OFF \
  -DCLSOCKET_DEP_ONLY=OFF \
  -DCLSOCKET_SHARED=OFF \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O3 -g3 -DNDEBUG" \
  -DCMAKE_C_FLAGS_RELWITHDEBINFO="-O3 -g3 -DNDEBUG" \
  -DCMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO="" \
  -DCMAKE_INSTALL_PREFIX=${install_dir} \
  -DINSTALL_CMAKE_FILES=OFF \
  -DINSTALL_EXTRA_LIBS=ON \
  -DINSTALL_LIBS=ON \
  -DUSE_CUSTOM_VECTOR_MATH=OFF \
  -DUSE_DOUBLE_PRECISION=ON \
  -DUSE_GLUT=OFF \
  -DUSE_GRAPHICAL_BENCHMARK=OFF \
  ../


# Bullet takes a while to build, so use parallel make with max number of
# processors.
platform=${1:-LINUX_gcc}
if [ "${platform}" = 'MACOS_gcc' ]; then
  numprocs=1
else
  numprocs=$(grep -c ^processor /proc/cpuinfo)
fi;
make install -j${numprocs}
