SHELL=/bin/bash
CC=mpicxx
OPT_LEVEL=-O3
OMP=

ifeq ($(USE_OMP),1)
  OMP=-fopenmp
endif

MPICXX=$(CC) -g $(OPT_LEVEL) -std=c++11 $(OMP)

bh:  simulation_bh.cpp
		$(MPICXX) -o n_body_simulation simulation_bh.cpp

strips:  simulation_strips.cpp
		$(MPICXX) -o n_body_simulation simulation_strips.cpp

strips_cutoff:  simulation_strips_cutoff.cpp
		$(MPICXX) -o n_body_simulation simulation_strips_cutoff.cpp

tiles:  simulation_tiles.cpp
		$(MPICXX) -o n_body_simulation simulation_tiles.cpp

tiles_cutoff:  simulation_tiles_cutoff.cpp
		$(MPICXX) -o n_body_simulation simulation_tiles_cutoff.cpp

clean:
	rm n_body_simulation
