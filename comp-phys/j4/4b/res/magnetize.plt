reset
set terminal tikz size 10cm, 8cm
set output "fit_magnetize.tex"
#set terminal wxt
#set size square
set samples 10000

set xlabel "temp"
#set ylabel "$\\chi$"
set ylabel "magnetize"

set xrange [0.0:4.0]
#set yrange [-40:160]
set yrange [-1.2:1]
#set logscale y

set key left top box
set key height 1
set key width 1
set key spacing 1.5

#set arrow from 2.269,0 to 2.269,120 nohead dt 2
#set label 1 left at first 2.34,-2.11 "temp=2.269"

#a=10
#c=2.3
#n=-1.8
#f(x) = (x>c) ? a*(x-c)**n : 0

a=-1
c=2.1
n=0.8
f(x) = a*(-x+c)**n

#fit [2.3:4.0] f(x) "temp-quantity.dat" u 1:8 via a,c,n
fit [1.3:2.0] f(x) "temp-quantity.dat" u 1:6 via a,c,n

#plot "temp-quantity.dat" u 1:6:7 lt 2 lc 6 with errorbars title "Simulation"#, "exact.dat" u 1:3 lc 2 with lines title "Exact"
#plot "temp-quantity.dat" u 1:8:9 lt 2 lc 6 title "data point" with errorbars,\
#f(x) title "Fitting Curve"
plot "temp-quantity.dat" u 1:6:7 lt 2 lc 6 title "data point" with errorbars,\
f(x) title "Fitting Curve"
