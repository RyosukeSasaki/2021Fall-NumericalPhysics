FC = gfortran
e1.x: eigen1.o
	${FC}  -o e1.x eigen1.o -lm -llapack
eigen1.o: eigen1.f90
	${FC} -c eigen1.f90
h.x: hist.o
	${FC} hist.o -o h.x
hist.o: hist.f90
	${FC} -c hist.f90
clean:
	rm *.[ox]
