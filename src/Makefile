CC=gcc
CFLAGS=-lm

all : imanm exanm cgmembrane block2full lapdecomp

exanm : exanm.o membranmutil.o
	$(CC) $^ $(CFLAGS) -o $@

imanm : imanm.o membranmutil.o 
	$(CC) $^ $(CFLAGS) -o $@

cgmembrane : cgmembrane.o membranmutil.o
	$(CC) $^ $(CFLAGS) -o $@

block2full : block2full_free.o membranmutil.o
	$(CC) $^ $(CFLAGS) -o $@

lapdecomp : lapdecomp_free.o membranmutil.o lapsolve.o
	$(CC) $^ $(CFLAGS) -llapack -lblas -o $@

lapsolve.o : lapsolve.c membranmutil.h lapsolve.h
	$(CC) $< -c -o $@

%.o : %.c
	$(CC) $< -c -o $@

.PHONY : clean
clean :
	-rm *.o
