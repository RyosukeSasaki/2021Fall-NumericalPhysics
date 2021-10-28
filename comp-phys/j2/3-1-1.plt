reset
set terminal tikz size 10cm, 10cm
set output "3-a-1.tex"
#set terminal wxt
set size square

set xlabel "$t$"
set ylabel "$x$"

#set format x "$%.1t \times 10^{%T}$"

set xrange [0:10.1]
set yrange [-1.1:1.1]

set xtics 0,1,10
set mxtics 2
set ytics -1,0.5,1
set mytics 5

set key font ", 8"
set key bottom left box
set key height 1
set key width 1
set key spacing 1.5

plot "3-a-1.dat" every 1 u 1:2 title "数値解",\
"3-a-1.dat" every 1 u 1:3 title "解析解"