FILE2-3 = differential.f90 main2-3.f90
FILE3-A = differential.f90 main3-a.f90
FILE3-B-1 = differential.f90 main3-b-11.f90
FILE3-B-2 = differential.f90 main3-b-12.f90
FILE3-B-3 = differential.f90 main3-b-13.f90

2-3: $(FILE2-3)
	gfortran $(FILE2-3) -o main2-3.o

3-a: $(FILE3-A)
	gfortran $(FILE3-A) -o main3-a.o

3-b: $(FILE3-B-1) $(FILE3-B-2) $(FILE3-B-3)  
	gfortran $(FILE3-B-1) -o main3-b-11.o
	gfortran $(FILE3-B-2) -o main3-b-12.o
	gfortran $(FILE3-B-3) -o main3-b-13.o

run:
	@make 3-a
	./main3-a.o > 3-a-2.dat

clean:
	rm *.o *.mod