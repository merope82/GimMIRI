UNAME := $(shell uname)

LIBS = -L ../cfitsio/
INCLUDE = -I ../cfitsio/

NVCC = nvcc
CC   = gcc
CFLAGS = -O3 -D_GNU_SOURCE
LINKS = -lcfitsio -lm
NVFLAGS = --compiler-options -Wall -gencode arch=compute_35,code=\"sm_35,compute_35\" -rdc=true

NVOBJECTS = GimMIRI.o setup.o free.o sort.o airy.o evolve.o kernel.o read_dev.o
CCOBJECTS = fits.o error.o print.o read_in.o getMemorySize.o 

all: GimMIRI
	@if [ -d "../" ] ; then mv GimMIRI ../GimMIRI ; fi
	@echo "Compilation successful! Executable is in the top directory."

%.o : %.c imMIRI.h externC.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

%.o : %.cu imMIRI.h externC.h GimMIRI.cuh
	$(NVCC) $(CFLAGS) $(NVFLAGS) $(INCLUDE) -c $<

GimMIRI: $(CCOBJECTS) $(NVOBJECTS) imMIRI.h GimMIRI.cuh externC.h
	$(NVCC) $(CFLAGS) $(NVFLAGS) $(INCLUDE) $(LIBS) -o GimMIRI $(CCOBJECTS) $(NVOBJECTS) $(LINKS)

clean:
	rm -f *.o

nuke:
	rm -f *.o ../GimMIRI
                        