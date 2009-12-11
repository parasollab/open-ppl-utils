include ../makefile_includes/Makefile.defaults

#INCLUDE = -I/usr/local/include
#--------------------------------------------------------------------
INCL = -I. -I/usr/local/include

OBJECTS	= Basic.o EulerAngle.o  Gauss.o  Matrix.o  Quaternion.o

CLEAN = example ${OBJECTS} *.a *.o *~ *.so

#default target
default_target: library

library: libmathtool.a

install: library example

libmathtool.a: ${OBJECTS}
	/bin/rm -f libmathtool.a
	${AR} libmathtool.a ${OBJECTS}

.SUFFIXES: .o .cpp
#.cpp.o:
	#@${MAKE} -s check_platform
	#${CXX} -c ${CXXFLAGS} ${OPTS} ${DEFS} ${INCL} $< -o $@

example: example.cpp
	${CXX} ${CXXFLAGS} ${OPTS} ${DEFS} ${INCL} $< -o $@ -L. -lmathtool -lm

