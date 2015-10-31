LINKS= -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread -lm

all: test svd

test: test.c
	icc test.c -o test $(LINKS)
	
svd: svd.c
	icc svd.c -o svd $(LINKS)

run: test
	./test

clean:
	rm test
	rm svd
