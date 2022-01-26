reset
set terminal tikz size 10cm, 7cm
set output "6A4.tex"
#set terminal wxt
set size square

set samples 10000
set xlabel "$j$"
set ylabel "$\\omega/\\omega_0$"

#set format x "$%.1t \times 10^{%T}$"
set xrange [-15:15]
#set yrange [0:0]

set xtics -15,5,15
set mxtics 5
#set ytics 0,0,0
#set mytics 2

#set key left top box
#set key height 1
#set key width 1
#set key spacing 1.5

plot "6a1g.dat" u 1:2 notitle lt 7 lc -1, 2*abs(sin(3.1415*x/30)) notitle lc -1 dt 2
