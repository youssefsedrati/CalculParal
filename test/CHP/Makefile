DEBUG = gfortran -g -O2 -Wall -Wextra -fcheck=all -fimplicit-none -fbounds-check -fbacktrace -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file

PROD = gfortran -O3 -march=native -fimplicit-none -Wall -Wline-truncation -fwhole-file

all:run

run:fonct.f90 grad.f90 main.f90
	$(PROD) fonct.f90 grad.f90 main.f90 -o run
debug:fonct.f90 grad.f90 main.f90
	$(DEBUG) fonct.f90 grad.f90 main.f90 -o run
clear:
	rm -R sol
	mkdir sol
