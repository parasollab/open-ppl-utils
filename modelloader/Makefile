#makefile for modelloader library
include ../makefile_includes/Makefile.defaults

#--------------------------------------------------------------------
INCL = -I. -I/usr/local/include -I../mathtool/

OBJS = IModel.o ModelFactory.o MovieBYULoader.o ObjLoader.o BVHDataLoader.o
TARGET = libmodloader.a

all: $(TARGET)
#--------------------------------------------------------------------
default_target: library

library: $(TARGET)

$(TARGET): $(OBJS)
	ar ruv $@  $(OBJS)

CLEAN = $(TARGET) $(OBJECTS) *.a *.o *~ *.so
