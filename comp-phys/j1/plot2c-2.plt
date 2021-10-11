reset
#set term wxt
set term tikz size 10cm, 10cm
set output "2c-2.tex"

set tics format "$10^{%T}$"
set xlabel '$h^4$'
set ylabel 'Error'
set logscale xy

a=1
f(x) = a*x
fit f(x) "output2c-2" u ($1**4):2 via a

plot "output2c-2" u ($1**4):2 notitle lt 1