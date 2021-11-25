reset
unset multiplot
set terminal tikz size 6cm, 16cm
#set terminal postscript eps
num = "4"
set output "3-b-".num.".tex"
#set terminal wxt
set parametric

set multiplot layout 3,1

set size square
set xlabel "$\\phi$"
set ylabel "$\\theta$"
set zlabel "$\\omega$"
set xrange [0:5.5]
set yrange [-3:3]
set zrange [-3:3]
set xtics 0,1,6
set mxtics 2
set ytics -3,1,3
set mytics 2
set ztics -3,3,3
set mztics 3
splot "3-b-".num."1.dat" notitle

set size square
set xlabel "$\\theta$"
set ylabel "$\\omega$"
set xrange [-3:3]
set yrange [-3:3]
set xtics -3,1,3
set mxtics 2
set ytics -3,1,3
set mytics 2
plot "3-b-".num."1.dat" u 2:3 notitle

set size square
set xlabel "$\\theta$"
set ylabel "$\\omega$"
set xrange [-3:3]
set yrange [-3:3]
set xtics -3,1,3
set mxtics 2
set ytics -3,1,3
set mytics 2
plot "3-b-".num."3.dat" u 2:3 notitle

unset multiplot