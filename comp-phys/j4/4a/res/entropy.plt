reset
set terminal tikz size 10cm, 8cm
set output "entropy.tex"
#set terminal wxt
#set size square

set xrange [0:601]

set xlabel "iseed"
set ylabel "sume"

a=4e-7
b=-1.9
f(x) = a*x+b
fit f(x) "entropy.dat" via a, b

plot "entropy.dat" lt 6 notitle, f(x) lc -1 dt 2 title "$y=6.1\\times10^{-8}-1.9,\\ R=0.02$"

#set output
#set terminal wxt
#replot