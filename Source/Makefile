CCFLAGS = -lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -g
TARGET = bu2s

#
# all is the default target that gets built when you just type "make"
#
all: ${TARGET}

CC = gcc

#
# Rule for linking together the object files needed for our program
#
${TARGET}: ${TARGET}.o MT/dSFMT.o PoissonLookup.o
	${CC} ${CCFLAGS} -o $@ MT/dSFMT.o PoissonLookup.o ${TARGET}.o -lm

#
# Rule for how the MT code gets build.  Just go to that directory and
# run make from there.
#
MT/dSFMT.o: 
	cd MT; make dSFMT.o

PoissonLookup.o:
	${CC} -DDSFMT_MEXP=19937 -I/usr/local/include/gsl -o PoissonLookup.o PoissonLookup.c

#
# Default rule for how a .o (object) file gets built from a .c file
#
.c.o:
	${CC} ${CCFLAGS} -DDSFMT_MEXP=19937 -DHAVE_SSE2 -c $<

#
# Rule for cleaning up everything associated with a build: The
# executable program, as well as any object files
#
clean:
	rm -f ${TARGET}
	rm -f ${TARGET}.o
	cd MT; make clean

