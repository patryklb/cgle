#!/usr/bin/gnuplot -persist

set multiplot layout 2,2 columnsfirst
#----------------------------------GRAPH 1-------------------------------------------#
set title "Norm"
set xlabel "t"
set ylabel "Norm"
plot 'stats.txt' using 1:2 w lines
#------------------------------------------------------------------------------------#

#----------------------------------GRAPH 2-------------------------------------------#
set title "Averaged noise (real part)"
set xlabel "t"
set ylabel "Average"
plot 'stats.txt' using 1:3 w lines
#------------------------------------------------------------------------------------#

#----------------------------------GRAPH 3-------------------------------------------#
set title "Noise variance"
set xlabel "t"
set ylabel "Variance"
plot 'stats.txt' using 1:5 w lines
#------------------------------------------------------------------------------------#
unset multiplot
