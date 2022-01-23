#all: i.x e.x ix.x
#all: i.x e.x
all: i.x
i.x: calc.o ising.o rnd.o
	gfortran calc.o ising.o rnd.o -o i.x
ising.o: ising.f90
	gfortran -fno-underscoring -c ising.f90
rnd.o: rnd.f90
	gfortran -fno-underscoring -c rnd.f90
e.x: exact.o calc.o
	gfortran exact.o calc.o -o e.x
exact.o: exact.f90
	gfortran -fno-underscoring -c exact.f90
calc.o: calc.f90
	gfortran -fno-underscoring -c calc.f90
clean:
	rm *.[ox]
