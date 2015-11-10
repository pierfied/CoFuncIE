EXEC = main.cfie

OBJS = main.o data_handler.o cosmology.o

CC = gcc

CFLAGS = 

OFLAGS = -lm

INCL = 

LIBS = 

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~