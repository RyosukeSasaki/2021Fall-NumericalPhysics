all: i.x e.x ix.x
i.x: ising.o rnd.o
	gfortran ising.o rnd.o -o i.x
ising.o: ising.f90
	gfortran -fno-underscoring -c ising.f90
rnd.o: rnd.f90
	gfortran -fno-underscoring -c rnd.f90
e.x: exact.o
	gfortran exact.o -o e.x
exact.o: exact.f90
	gfortran -c exact.f90
ix.x: isingx.o rnd.o xsub.o
	gfortran isingx.o rnd.o xsub.o -lX11 -o ix.x
isingx.o: isingx.f90
	gfortran -fno-underscoring -c isingx.f90
xsub.o: xsub.c
	gcc -I/usr/include/X11 -c xsub.c
clean:
	rm *.[ox]
