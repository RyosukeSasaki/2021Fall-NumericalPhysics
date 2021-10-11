reset
#set term wxt
set term tikz size 10cm, 10cm
set output "2c-1.tex"

set xlabel '$T/\Theta_D$'
set ylabel '$C_V/3Nk_B$'
set yrange [0:1.5]

plot "output2c-1" u 1:2 with lines title "Debye",\
       "output2c-1" u 1:3 with lines title "Dulong-Petit",\
       "output2c-1" u 1:4 with lines title "低温での振る舞い",\

