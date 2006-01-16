CC = g++

CFLAGS		= -O2 -I.

OBJECTS		= lib/PQP.o lib/BV.o lib/Build.o lib/TriDist.o

CLEAN		= $(OBJECTS) lib/libPQP.a include/*.h

library: lib/libPQP.a

lib/libPQP.a: $(OBJECTS)
	ar ruv lib/libPQP.a $(OBJECTS)
	cp src/*.h include/.

lib/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	/bin/rm -f $(CLEAN)
