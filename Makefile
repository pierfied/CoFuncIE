EXEC = main.cfie

OBJS = main.o data_handler.o cosmology.o density_map.o map_likelihood.o

ifdef COMPILER
ifeq ($(COMPILER),intel)
$(info Using intel compiler.)
CC = icc
CFLAGS = -openmp
OFLAGS = -openmp -lpthread -lm
endif
else
CC = gcc
CFLAGS = -fopenmp
OFLAGS = -fopenmp -lm
endif

INCL = -I./gsl/include/

LIBS = -lgsl -lgslcblas -L./gsl/lib -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64\
 -lmkl_intel_thread -lmkl_core

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(INCL) -c $< -o $@ $(CFLAGS)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC) $(OFLAGS)
	mv $(OBJS) Objects/

.PHONY: clean

clean:
	-rm -f $(OBJS) $(EXEC) *~
	-rm -f Objects/*.o
