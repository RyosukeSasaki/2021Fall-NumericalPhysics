reset
set terminal tikz size 10cm, 7cm
set output "6B6.tex"
#set terminal wxt
#set size square
pi = 4.0*atan(1.0)

set xlabel "$\\omega/\\omega_0$"
set ylabel "Density of States"

#set format x "$%.1t \times 10^{%T}$"
set xrange [0:2]
set yrange [0:4.5]

set xtics 0,0.5,2
set mxtics 5
set ytics 0,1,5
set mytics 2
set boxwidth 0.035

set key left top box
set key height 1
set key width 3
set key spacing 1.5

#f(x) = 1/(pi*sqrt(1-x**2/4))

set style fill solid 0.5

plot\
 "6b1a-h.dat" u ($1-.0157):2 with boxes title "$\\gamma=1$",\
 "6a1e.dat"   u ($1+.0157):2 with boxes title "Result of 6-A.",\
