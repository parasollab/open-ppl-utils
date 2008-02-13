include ../../../../makefile_includes/Makefile.defaults

CLEAN = *.o ../core qconvex qdelaunay qhalf qvoronoi qhull *~ libqhull.a ../rbox ../qhull ../qhalf ../qconvex ../qdelaunay ../qvoronoi ../user_eg ../user_eg2 ../*.exe

OBJS = user.o global.o stat.o io.o geom2.o poly2.o merge.o qhull.o geom.o poly.o qset.o mem.o

HFILES = user.h qhull.h qhull_a.h geom.h io.h mem.h merge.h poly.h qset.h stat.h

LIBS = -L. -lqhull -lm

default_target: libqhull.a

all: rbox qconvex qdelaunay qhalf qvoronoi qhull

new: clean all

libqhull.a: ${OBJS}
	${AR} libqhull.a ${OBJS}

qconvex: qconvex.o libqhull.a
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ $@.o ${LIBS}
	cp $@ ..

qdelaunay: qdelaun.o libqhull.a
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ qdelaun.o ${LIBS}
	cp $@ ..

qhalf: qhalf.o libqhull.a
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ $@.o ${LIBS}
	cp $@ ..

qvoronoi: qvoronoi.o libqhull.a
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ $@.o ${LIBS}
	cp $@ ..

qhull: unix.o libqhull.a
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ unix.o ${LIBS}
	cp $@ ..
	-chmod +x ../eg/q_test ../eg/q_eg ../eg/q_egtest
	-cd ..; ./rbox D4 | ./qhull

# compile qhull without using libqhull.a
qhullx: qconvex.o qdelaun.o qhalf.o qvoronoi.o unix.o ${OBJS}
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o qconvex qconvex.o ${OBJS} -lm 
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o qdelaunay qdelaun.o ${OBJS} -lm 
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o qhalf qhalf.o ${OBJS} -lm 
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o qvoronoi qvoronoi.o ${OBJS} -lm 
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o qhull unix.o ${OBJS} -lm 
	cp qconvex qdelaunay qhalf qvoronoi qhull ..
	-chmod +x ../eq/q_test ../eq/q_eg ../eq/q_egtest
	-cd ..; ./rbox D4 | ./qhull

rbox: rbox.o
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ $< -lm
	cp $@ ..

user_e%: user_e%.o libqhull.a 
	${CC} ${CFLAGS} ${OPTS} ${DEFS} ${INCL} -o $@ $< ${LIBS}
	cp $@ ..

unix.o: qhull.h user.h mem.h
qconvex.o: qhull.h user.h mem.h
qdelaun.o: qhull.h user.h mem.h
qhalf.o: qhull.h user.h mem.h
qvoronoi.o: qhull.h user.h mem.h
qhull.o: $(HFILES)
geom.o: $(HFILES)
geom2.o: $(HFILES)
global.o: $(HFILES)
io.o: $(HFILES)
mem.o: mem.h 
merge.o: $(HFILES)
poly.o: $(HFILES)
poly2.o: $(HFILES)
qset.o: qset.h mem.h 
stat.o: $(HFILES)
user.o: $(HFILES)

