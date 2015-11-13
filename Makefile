EXEC = main.cfie

OBJS = main.o data_handler.o cosmology.o

CC = gcc

CFLAGS = 

OFLAGS = -lm

INCL = 

LIBS = 

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(INCL) -c $< -o $@ $(CFLAGS)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC) $(OFLAGS)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~
