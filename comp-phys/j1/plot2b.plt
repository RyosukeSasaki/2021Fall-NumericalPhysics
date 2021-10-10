#set term wxt
set term tikz size 10cm, 10cm
set output "2b.tex"

set yrange [0:2]

plot "output2b" u 1:2 with lines title "Planck",\
       "output2b" u 1:3 with lines title "Rayleigh Jeans",\
       "output2b" u 1:4 with lines title "Wien",\

