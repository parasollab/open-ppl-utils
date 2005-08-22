# This Makefile builds two example programs that do not use Inventor.

#CFLAGS=  -O2 -I ../include -I
CFLAGS=   -I ../include -I/pub/dsmft/STL 
#lkdlkd
LDFLAGS= -o32 -L ../lib -rpath ../lib -lvclip -lm

CC = g++

cubes: cubes.C
	${CC} $(CFLAGS) cubes.C $(LDFLAGS) -o cubes

file: fromfile.C
	${CC} $(CFLAGS) fromfile.C $(LDFLAGS) -o file

verts: fromverts.C
	${CC} $(CFLAGS) fromverts.C $(LDFLAGS) -o verts


print:
	${CC} $(CFLAGS) print.C $(LDFLAGS) -o print
	-rm -r ii_files

allpairs:
	${CC} $(CFLAGS) allpairs.C $(LDFLAGS) -o allpairs
	-rm -r ii_files

clean:
	-rm -f cubes
	-rm -f verts
	-rm -f file
	-rm -f print allpairs
	-rm -f -r ii_files
