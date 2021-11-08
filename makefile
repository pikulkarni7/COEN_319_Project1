CC = gcc

default:  matmult_omp

matmult_omp: matmult_omp.c
	${CC} -O3 -Wall -Wextra -fopenmp -o $@ matmult_omp.c


clean:
	
	-rm -f matmult_omp
