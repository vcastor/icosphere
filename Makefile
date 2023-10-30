gg = gfortran
FLAG = -fcheck=all -Wall -O3
TARGETS = icosphere.exe brute_force.exe

all: $(TARGETS)

SRCF90 = $(wildcard *.f90)

%.exe: %.f90
	$(gg) $(FLAG) $< -o $@

clean:
	rm -f *.exe *.o
