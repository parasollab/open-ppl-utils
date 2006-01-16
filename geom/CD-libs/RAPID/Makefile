
CC = g++

#CFLAGS		= -I. -I..

# IRR
CFLAGS		= -I. -I.. -O2

LDFLAGS		= -L. -L..

.SUFFIXES: .C

SRCS		= RAPID.C build.C collide.C overlap.C 

OBJECTS		= RAPID.o build.o collide.o overlap.o 

CLEAN		= sample_client $(OBJECTS) *.a *.o *~ *.so

default: library test

install: libRAPID.a libRAPID.so
library: libRAPID.a libRAPID.so

libRAPID.a: $(OBJECTS)
	/bin/rm -f libRAPID.a
	ar ruv libRAPID.a $(OBJECTS)

libRAPID.so: $(OBJECTS)
	${CC} -shared ${CFLAGS} $^  -o libRAPID.so

.C.o:
	$(CC) ${CFLAGS} -c $<

.c.o:
	$(CC) ${CFLAGS} -c $<

test:
	$(CC) -o sample_client sample_client.C -L. -lRAPID -lm

clean:
	touch clean.a
	/bin/rm -f $(CLEAN)
