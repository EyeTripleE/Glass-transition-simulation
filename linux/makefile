SHELL = /bin/bash

CC=mpicxx

OPT_LEVEL = -O4

MPICXX=$(CC) -g -O4 -std=c++11 -fopenmp

all: n_body_simulation

n_body_simulation:  simulation_mpi.cpp
		$(MPICXX) -o n_body_simulation simulation_mpi.cpp

clean:
	rm n_body_simulation 
