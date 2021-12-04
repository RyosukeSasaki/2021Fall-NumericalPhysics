reset
set terminal tikz size 10cm, 8cm
set output "4af2.tex"
set size square

set xlabel "$\\Delta x^2$"
set ylabel "$\\langle (x-\\langle x\\rangle)^2\\rangle$"
#set ylabel "$\\langle x\\rangle$"
set samples 10000

#set format x "$%.1t \\times 10^{%T}$"
set xrange [-0.1:2.5]
#set yrange [0:0]

set xtics -0.5,0.5,2.5
set mxtics 5
#set ytics 0,0,0
#set mytics 2

set key right top box
set key height 1
set key width 2
set key spacing 1.5

a=-0.5
b=0.07
f(x) = a*x + b

fit f(x) "../out/4af.dat" u (1/$1)**2*10**4:3 via a,b
#set label 3 left at first 0.7386,0.4949 "$y=-7.856\\times10^{-3}x+0.5000$" font ",8"
set label 3 left at first 0.7050,0.07060 "$y=-2.541\\times10^{-3}x+0.07247$" font ",8"

plot \
"4af.dat" u 1:3 title "厳密解",\
"../out/4af.dat" u (1/$1)**2*10**4:3 title "数値解",\
f(x) title "Fitting Line"\