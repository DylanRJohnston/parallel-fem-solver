#CC=icc
#OPT_FLAGS=-march=native -O3 -qopt-report=5
#Uncomment to use GCC
CC=gcc-5
OPT_FLAGS=-march=native -Wa,-q -O3 #-fopt-info-all-optall
CC_FLAGS=-std=c99 -g -Wall -Werror -Wextra -pedantic -fopenmp

all: main

obj/%.o: %.c
	$(CC) $(CC_FLAGS) $(OPT_FLAGS) $^ -c -o $@ #-qopt-report-file=$@.optrpt

main: obj/main.o
	$(CC) $(CC_FLAGS) $(OPT_FLAGS) $^ -o $@

clean:
	rm -f *.optrpt obj/* main
