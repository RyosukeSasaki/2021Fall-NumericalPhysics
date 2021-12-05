reset
set terminal tikz size 10cm, 8cm
set output "4ba3.tex"
#set terminal wxt
set size square

set xlabel "$\\Delta x^2$"
set ylabel "$P_{\\rm center}$"
set samples 10000

set xrange [-0.1:2.5]
set yrange [0.013:0.02]

#set xtics -0.5,0.5,2.5
#set mxtics 5
#set ytics 0,0,0
#set mytics 2

set key right top box
set key height 1
set key width 2
set key spacing 1.5

a=0.003
b=0.6
f(x) = a*x + b

fit [0:1] f(x) "../out/4ba.dat" u ($1)**2*10**4:3 via a,b
#set label 3 left at first 0.55206,0.316437 "$y=-1.061\\times10^{-2}x+0.3214$" font ",8"
#set label 3 left at first 0.5380,0.6657 "$y=-2.500\\times10^{-3}x+0.6647$" font ",8"
set label 3 left at first 0.52147,0.016248 "$y=7.809\\times10^{-3}x+0.01393$" font ",8"

plot \
"../out/4ba.dat" u ($1)**2*10**4:3 title "数値解",\
f(x) title "Fitting Line"\