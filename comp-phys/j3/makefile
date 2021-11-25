all: w.x wx.x
w.x: wave.o
	gfortran wave.o -o w.x
wx.x: wavex.o xsub.o
	gfortran wavex.o xsub.o -lX11 -o wx.x
wave.o: wave.f90
	gfortran -c wave.f90
wavex.o: wavex.f90
	gfortran -fno-underscoring -c wavex.f90
xsub.o: xsub.c
	gcc -I/usr/include/X11 -c xsub.c
clean:
	rm *.[ox]
