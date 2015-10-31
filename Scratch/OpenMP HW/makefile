all: par_convert_cart convert_cart

par_convert_cart: par_convert_cart.c
	gcc par_convert_cart.c -o par_convert_cart -lm -fopenmp

convert_cart: convert_cart.c
	gcc convert_cart.c -o convert_cart -lm

clean:
	rm par_convert_cart
	rm convert_cart
