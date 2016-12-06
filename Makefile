CC=icc
OPT_FLAGS=-march=native -O3 -qopt-report=5
#Uncomment to use GCC
#CC=mpicc
#OPT_FLAGS=-march=native -Wa,-q -O3 #-fopt-info-all-optall
CC_FLAGS=-std=c99 -g -Wall -Werror -Wextra -pedantic -fopenmp -lmpi

all: fem1d

obj/%.o: %.c
	$(CC) $(CC_FLAGS) $(OPT_FLAGS) $^ -c -o $@ #-qopt-report-file=$@.optrpt

fem1d: obj/fem1d_new.o
	$(CC) $(CC_FLAGS) $(OPT_FLAGS) $^ -o $@

clean:
	rm -f *.optrpt obj/* fem1d fem1d.dSYM
