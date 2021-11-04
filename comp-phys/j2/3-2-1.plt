reset
set terminal tikz size 10cm, 10cm
set output "3-b-41.tex"
#set terminal wxt
set size square
set parametric

set xlabel "$\\phi$"
set ylabel "$\\theta$"
set zlabel "$\\omega$"

#set format x "$%.1t \times 10^{%T}$"

set xrange [0:5.5]
set yrange [-3:3]
set zrange [-3:3]

set xtics 0,1,6
set mxtics 2
set ytics -3,1,3
set mytics 2
set ztics -3,1,3
set mztics 2

#set key left top box
#set key height 1
#set key width -1
#set key spacing 1.5

splot "3-b-41.dat" notitle