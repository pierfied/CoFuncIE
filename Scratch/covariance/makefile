EXEC = cov 

OBJS = main.o 

CC = gcc 

INCL = -I$HOME/code/gsl/include 

LIBS = -lgsl -lgslcblas -L$HOME/code/gsl/lib

.SUFFIXES : .c .o 

%.o: %.c 
	$(CC) $(INCL) -c $< -o $@ 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC) 

#$(OBJS): $(INCL) 

.PHONY : clean 

run: $(EXEC)
	./$(EXEC)

time: $(EXEC)
	time ./$(EXEC)

clean: 
	rm -f $(OBJS) $(EXEC)