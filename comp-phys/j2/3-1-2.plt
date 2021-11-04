reset
set terminal tikz size 10cm, 10cm
set output "3-a-2.tex"
#set terminal wxt
set size square

set xlabel "$\\tau$"
set ylabel "Error"

#set format x "$%.1t \times 10^{%T}$"

#set xrange [0:10.1]
#set yrange [-1.1:1.1]

#set xtics 0,1,10
#set mxtics 2
#set ytics -1,0.5,1
#set mytics 5

#set key font ", 8"
#set key bottom left box
#set key height 1
#set key width 1
#set key spacing 1.5

set logscale xy

n=4
a=0.001
err_x(x) = a*x**n

fit [0.0019:0.01] err_x(x) "3-a-2.dat" u 1:3 via a, n

plot "3-a-2.dat" every 1 u 1:3 title "Error $x$",\
err_x(x)
#"3-a-2.dat" every 1 u 1:4 title "Error $v$",\