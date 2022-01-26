reset
set terminal tikz size 10cm, 7cm
set output "6A3.tex"
#set terminal wxt
set size square
pi = 4.0*atan(1.0)

set xlabel "$\\omega/\\omega_0$"
set ylabel "Density of States"

#set format x "$%.1t \times 10^{%T}$"
set xrange [0:2]
set yrange [0:2.5]

#set xtics 0,0,0
#set mxtics 2
set ytics 0,1,3
set mytics 2

set key left top box
set key height 1
set key width 1
set key spacing 1.5

f(x) = 1/(pi*sqrt(1-x**2/4))

plot "6a1f.dat" with boxes title "Numerical", f(x) title "Theoritical"
