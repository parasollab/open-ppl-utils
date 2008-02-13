include ../../../Makefile.defaults

QHULL_DIR  = ../qhull3.0/src

DEFS = -DQHULL

INCL = -I include -I ${QHULL_DIR}

LIBS = -L${QHULL_DIR} -lqhull

CLEAN = *.o lib/libvclip.a ii_files; cd examples; ${MAKE} clean; cd ../iv_examples; ${MAKE} clean

default_target: lib/libvclip.a

all: lib/libvclip.a

%.o: src/%.C
	@${MAKE} -s check_platform
	${CXX} -c ${CXXFLAGS} ${OPTS} ${DEPS} ${INCL} $< -o $@

lib/libvclip.a: mv.o PolyTree.o vclip.o
	${AR} $@ mv.o PolyTree.o vclip.o ${QHULL_DIR}/libqhull.a

examples:
	cd $@; ${MAKE}

iv_examples:
	cd $@; ${MAKE}

