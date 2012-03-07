#makefile for modelloader library
include ../makefile_includes/Makefile.defaults

#--------------------------------------------------------------------
INCL = -I. -I/usr/local/include -I../mathtool/

OBJS = IModel.o ModelFactory.o  ModelTool.o  MovieBYULoader.o  ObjLoader.o BVHDataLoader.o 
TARGET = libmodloader.a
all: $(TARGET)
objs: ${OBJS}
#--------------------------------------------------------------------
default_target: library

library: $(TARGET)

$(TARGET): objs
	ar ruv $@  ${OBJS}

#clean:
#	rm -f $(OBJS) $(TARGET)
CLEAN = ${TARGET} ${OBJECTS} *.a *.o *~ *.so



