# Must set these gl and glut locations to build 'spinning'

GLUT_CFLAGS  = -I/afs/cs.unc.edu/project/walk/release/include -I../Xw
GLUT_LDFLAGS =  -L/usr/lib/X11R5 \
                -L/afs/cs.unc.edu/project/walk/release/lib/`hw_os` \
                -L/usr/X11R6/lib/
GLUT_LDLIBS  =  -lglut -lGL -lGLU -lXmu -lXext -lX11

.SUFFIXES: .cpp

CC = g++
CFLAGS  = -O2 -I. -I../../include $(GLUT_CFLAGS)
LDFLAGS	= -L. -L../../lib $(GLUT_LDFLAGS)
LDLIBS  = -lPQP -lm       $(GLUT_LDLIBS) 

SRCS    = main.cpp model.cpp 

OBJECTS	= main.o model.o

CLEAN   = $(OBJECTS) spinning 

.cpp.o:
	$(CC) ${CFLAGS} -c $<

default: $(OBJECTS)
	$(CC) $(CFLAGS) -o spinning $(OBJECTS) -L. $(LDFLAGS) $(LDLIBS)

clean:
	/bin/rm -f $(CLEAN)

