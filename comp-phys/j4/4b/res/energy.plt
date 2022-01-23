reset
#set terminal tikz size 10cm, 8cm
#set output "exact_energy.tex"
set terminal wxt
#set size square

set xlabel "temp"
#set ylabel "energy"
set ylabel "heat capacity"

set key left top box
set key height 1
set key width 3
set key spacing 1.5

#set arrow from 2.269,-2.2 to 2.269,-0.4 nohead dt 2
#set label 1 left at first 2.34,-2.11 "temp=2.269"

plot "temp-quantity.dat" u 1:4:5 lt 2 lc 6 with errorbars title "Simulation", "exact.dat" u 1:3 lc 2 with lines title "Exact"

