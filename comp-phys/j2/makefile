FILE2-3 = differential.f90 main2-3.f90
FILE3-A = differential.f90 main3-a.f90
FILE3-B = differential.f90 main3-b.f90

2-3: $(FILE2-3)
	gfortran $(FILE2-3) -o main2-3.o

3-a: $(FILE3-A)
	gfortran $(FILE3-A) -o main3-a.o

3-b: $(FILE3-B)
	gfortran $(FILE3-B) -o main3-b.o

clean:
	rm *.o *.mod