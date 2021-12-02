#set terminal tikz size 10cm, 8cm
#set output "4ae2.tex"
set terminal wxt
set size square

set xlabel "$t$"
#set ylabel "$\\langle x\\rangle$"
set ylabel "$\\langle (x-\\langle x\\rangle)^2\\rangle$"

#set format x "$%.1t \times 10^{%T}$"

#set xrange [0:0]
#set yrange [0:0]

#set xtics 0,0,0
#set mxtics 2
#set ytics 0,0,0
#set mytics 2

#set key right bottom box
set key left top box
set key height 1
set key width 2
set key spacing 1.5

plot\
"../out/4a1.dat" u 1:3 title "$\\Delta x=1/64$",\
"../out/4a2.dat" u 1:3 title "$\\Delta x=1/128$",\
"../out/4a3.dat" u 1:3 title "$\\Delta x=1/256$",\
"../out/4ad.dat" u 1:3 with lines title "厳密解"