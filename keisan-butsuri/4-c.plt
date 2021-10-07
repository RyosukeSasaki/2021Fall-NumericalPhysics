set term postscript
set output "test1.ps"
set logscale xy

plot "test1" with linespoints,\
       "test2" with linespoints
