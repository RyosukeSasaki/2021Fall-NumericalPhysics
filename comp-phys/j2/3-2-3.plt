reset
set terminal tikz size 10cm, 10cm
set output "3-b-13.tex"
#set terminal wxt
set size square
set parametric

set xlabel "$\\theta$"
set ylabel "$\\omega$"

#set format x "$%.1t \times 10^{%T}$"

set xrange [-3:3]
set yrange [-3:3]

set xtics -3,1,3
set mxtics 2
set ytics -3,1,3
set mytics 2

#set key left top box
#set key height 1
#set key width -1
#set key spacing 1.5

plot "3-b-13.dat" u 2:3 notitle