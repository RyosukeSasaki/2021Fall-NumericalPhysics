reset
set terminal tikz size 12cm, 10cm
set output "4aADD1.tex"
#set terminal wxt

set xlabel "$x$"
set ylabel "Im"
set zlabel "Re"

#set view equal xyz
set ticslevel 0
set zeroaxis

set xrange [-4:4]
set yrange [-2:2]
#set zrange [-1.5:1.5]

set xtics -4,1,4
set mxtics 2
set ytics -3,1,3
set mytics 2
#set ztics -3,1,3
#set mztics 2

set key right top box
set key height 1
set key width 2
set key spacing 1.5

splot \
"../out/4aADD1.dat" u 1:(0):2 linewidth 0.5 title "$|\\psi(x)|^2$",\
"../out/4aADD1.dat" u 1:(0):3 linewidth 0.5 title "Re$(\\psi(x))$",\
"../out/4aADD1.dat" u 1:4:(0) linewidth 0.5 title "Im$(\\psi(x))$"