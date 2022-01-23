reset
#set terminal tikz size 10cm, 8cm
#set output ""
set terminal wxt
set size square

set xlabel "temp"
set ylabel "unchi"

#set format x "$%.1t \times 10^{%T}$"
#set xrange [2:4]
#set yrange [0:0]

#set xtics 0,0,0
#set mxtics 2
#set ytics 0,0,0
#set mytics 2

set key right top box
set key height 1
set key width 1
set key spacing 1.5
set nokey

a=8.0
b=16.0
c=32.0
eta = 2.05

f(u,v,l1,l2) = log(v/u)/log(l2/l1)

plot "all.dat" u 1:(f($2,$3,a,b)) with lines,\
"all.dat" u 1:(f($3,$4,b,c)) with lines,\

