include ../makefile_includes/Makefile.defaults

DEFS = -DTETLIBRARY
INCL =

OBJECTS	= predicates.o tetgen.o

CLEAN = ${OBJECTS} *.a *.o *~ *.so

#default target
default_target: library

library: libtet.a

libtet.a: ${OBJECTS}
	/bin/rm -f libtet.a
	${AR} libtet.a ${OBJECTS}

tetgen: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OPTS} ${LIBS} ${OBJECTS} -o $@ -lm

.SUFFIXES: .o .cxx
.cxx.o:
	@${MAKE} -s check_platform
	${CXX} -c ${CXXFLAGS} ${OPTS} ${DEFS} ${INCL} $< -o $@
