CC = gcc 

CPP = g++

GMP_INC_DIR=.
GMP_LIB_DIR=.

LIB_PATHS = -L/usr/local/lib/ -L/opt/OpenBLAS/lib
INC_PATHS = -I/usr/local/include -I/opt/OpenBLAS/include

CFLAGS = -std=c99 -O3 -msse -mavx -Wall -pedantic -funroll-loops --param max-unroll-times=2 -fomit-frame-pointer -fexpensive-optimizations  

LDFLAGS = $(DEBUG) $(LIB_PATHS) -static-libgcc -fopenmp  -lgmp  -lmpfr -lflint  -lm -lopenblas -lpthread
LDFLAGS0 = -Wall $(DEBUG) $(LIB_PATHS) -lgmp  -lmpfr -lflint  -lm 

flint:
	$(CC) -O3 -static main.c -o bresultant $(LDFLAGS) /usr/lib/x86_64-linux-gnu/libmpfr.a /usr/lib/x86_64-linux-gnu/libgmp.a /usr/lib/x86_64-linux-gnu/libflint.a 

