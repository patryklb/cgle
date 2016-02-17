#!/usr/bin/gnuplot -persist
#set terminal postscript portrait enhanced lw 1 'Helvetica' 14
#set output 'plot.ps'

set title "Correlation function"
set xlabel "d"
set ylabel "Corr(psi)"
plot [][0:1.1] 'corr.txt' using 1:2 w lines
