CCFLAGS = -O3
GSLFLAGS = -lm -I/usr/include/gsl -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas
TARGET = bu2s

#
# all is the default target that gets built when you just type "make"
#
all: ${TARGET}

CC = gcc

#
# Rule for linking together the object files needed for our program
#
${TARGET}: ${TARGET}.o MT/dSFMT.o PoissonLookup.o LDcalculations.o
	${CC} ${CCFLAGS} -o $@ MT/dSFMT.o PoissonLookup.o LDcalculations.o ${TARGET}.o -lgsl -lgslcblas -lm

#
# Rule for how the MT code gets build.  Just go to that directory and
# run make from there.
#
MT/dSFMT.o: 
	cd MT; make dSFMT.o; cd ..

PoissonLookup.o: PoissonLookup.c bu2s.h
	${CC} ${CCFLAGS} -c ${GSLFLAGS} -DDSFMT_MEXP=19937 -o PoissonLookup.o PoissonLookup.c

LDcalculations.o: LDcalculations.c bu2s.h
	${CC} ${CCFLAGS} -c -DDSFMT_MEXP=19937 -o LDcalculations.o LDcalculations.c

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
	rm -f PoissonLookup.o
	rm -f LDcalculations.o
	cd MT; make clean

