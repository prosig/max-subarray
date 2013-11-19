#!/bin/bash

# Makefile make me fast and accurate

CC = gcc
EXE_ARGS = ./input/test_input_01.in

CFLAGS = -g -O -Wno-error -std=c99 -pthread -fopenmp
LFLAGS = 

all: max-sub-array-sequential max-sub-array-pthread max-sub-array-openmp

max-sub-array-sequential: max-sub-array-sequential.c
	$(CC) $(CFLAGS) max-sub-array-sequential.c -o max-sub-array-sequential

max-sub-array-pthread: max-sub-array-pthread.c
	$(CC) $(CFLAGS) max-sub-array-pthread.c -o max-sub-array-pthread

max-sub-array-openmp: max-sub-array-openmp.c
	$(CC) $(CFLAGS) max-sub-array-openmp.c -o max-sub-array-openmp

clean:
	rm -rf max-sub-array-sequential
	rm -rf max-sub-array-pthread
	rm -rf max-sub-array-openmp

# End of Makefile



