reset
set terminal tikz size 10cm, 10cm
set output "3-a-2.tex"
#set terminal wxt
set size square

set xlabel "$\\tau$"
set ylabel "Error"

set format y "$10^{%T}$"

#set xrange [0:10.1]
#set yrange [-1.1:1.1]

#set xtics 0,1,10
#set mxtics 2
#set ytics -1,0.5,1
#set mytics 5

#set key font ", 8"
set key bottom right
set key nobox
#set key height 1
#set key width 1
#set key spacing 1.5

set logscale xy

n1=4
a1=0.001
n2=4
a2=0.001
err_x(x) = a1*x**n1
err_v(x) = a2*x**n2

fit [0.0019:0.01] err_x(x) "3-a-2.dat" u 1:3 via a1, n1
fit [0.0019:0.01] err_v(x) "3-a-2.dat" u 1:4 via a2, n2

plot "3-a-2.dat" every 1 u 1:3 title "Error $x$",\
"3-a-2.dat" every 1 u 1:4 title "Error $v$",\
err_x(x) title "Fitting Line of Error $x$",\
err_v(x) title "Fitting Line of Error $v$",\