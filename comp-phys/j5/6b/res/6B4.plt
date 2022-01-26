reset
set terminal tikz size 10cm, 8cm
set output "6B7.tex"
#set terminal wxt
set size square

set samples 10000

set xlabel "$k/a$"
set ylabel "$\\omega/\\omega_0$"

#set format x "$%.1t \times 10^{%T}$"
set xrange [-pi/2:pi/2]
#set yrange [0:0]

#set xtics 0,0,0
#set mxtics 2
#set ytics 0,0,0
#set mytics 2

#set key left top box
set key bottom outside
#set key height 1
set key width 3
#set key spacing 1.5

f1(x)=sqrt( ( (1+1) + sqrt((1+1)**2-4*1*(sin(x))**2) ) )
h1(x)=sqrt( ( (1+1) - sqrt((1+1)**2-4*1*(sin(x))**2) ) )

f2(x)=sqrt( ( (1+2) + sqrt((1+2)**2-4*2*(sin(x))**2) )/2 )
h2(x)=sqrt( ( (1+2) - sqrt((1+2)**2-4*2*(sin(x))**2) )/2 )

f4(x)=sqrt( ( (1+4) + sqrt((1+4)**2-4*4*(sin(x))**2) )/4 )
h4(x)=sqrt( ( (1+4) - sqrt((1+4)**2-4*4*(sin(x))**2) )/4 )

f8(x)=sqrt( ( (1+8) + sqrt((1+8)**2-4*8*(sin(x))**2) )/8 )
h8(x)=sqrt( ( (1+8) - sqrt((1+8)**2-4*8*(sin(x))**2) )/8 )

plot \
f8(x) lc 4 notitle, h8(x) lc 1 title "$\\gamma=8$",\
f4(x) lc 3 notitle, h4(x) lc 2 title "$\\gamma=4$",\
f2(x) lc 2 notitle, h2(x) lc 3 title "$\\gamma=2$",\
f1(x) lc 1 notitle, h1(x) lc 4 title "$\\gamma=1$",\
