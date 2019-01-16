include ../../utils/makefile_includes/Makefile.defaults

INCL = -I../gsl-1.14/include

OBJECTS = kabsch.o

CLEAN = kabsch.o libkabsch.a

default_target: library

kabsch.o: kabsch.h kabsch.cpp
	@${MAKE} -s check_platform
	${CXX} -c ${CXXFLAGS} ${OPTS} ${DEFS} ${INCL} -o kabsch.o kabsch.cpp

library: ${OBJECTS}
	/bin/rm -f libkabsch.a
	${AR} libkabsch.a ${OBJECTS}

