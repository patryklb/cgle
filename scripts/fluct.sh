#!/usr/bin/gnuplot -persist
#set terminal postscript portrait enhanced lw 1 'Helvetica' 14
#set output 'plot.ps'


#set isosamples 100
#set samples 100
#unset key
#set title "Initial conditions"
#splot 'output.txt' using 1:2:3 with pm3d


set title "Fluctuations averaged over time"
set xlabel "x"
set ylabel "|n-n0|/n0"
plot [][] 'fluct.txt' using 1:2 w lines title '' lt 1, 'fluct.txt' using 1:3 w lines title '' lt 2,  'fluct.txt' using 1:4 w lines title '' lt 3
#------------------------------------------------------------------------------------#

unset multiplot
