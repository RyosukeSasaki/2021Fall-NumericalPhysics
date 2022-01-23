reset
#set terminal tikz size 10cm, 8cm
#set output "x.tex"
set terminal wxt
set size square

set xlabel "temp"
set ylabel "$\\chi$"

plot "temp-quantity.dat" u 1:8:9 lt 2 lc 6 with errorbars notitle

#set output
#set terminal wxt
#replot