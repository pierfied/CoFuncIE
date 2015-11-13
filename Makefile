EXEC = main.cfie

OBJS = main.o data_handler.o cosmology.o

ifdef COMPILER
ifeq ($(COMPILER),intel)
$(info Using intel compiler.)
CC = icc
CFLAGS = -openmp -lm
OFLAGS = -openmp -lm
endif
else
CC = gcc
CFLAGS = -fopenmp
OFLAGS = -fopenmp -lm
endif

INCL = 

LIBS = 

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
