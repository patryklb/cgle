#!/usr/bin/gnuplot -persist
#set terminal postscript portrait enhanced lw 1 'Helvetica' 14
#set output 'plot.ps'


#set isosamples 100
#set samples 100
#unset key
#set title "Initial conditions"
#splot 'output.txt' using 1:2:3 with pm3d

set multiplot layout 2,2 columnsfirst
#----------------------------------GRAPH 1-------------------------------------------#
set title "Module"
set xlabel "x"
set ylabel "|Psi|"
plot [][] 'initial.txt' using 1:2 w lines title '' lt 1
#------------------------------------------------------------------------------------#

#----------------------------------GRAPH 2-------------------------------------------#
set title "Phase"
set xlabel "x"
set ylabel "|Phi|"
plot [][] 'initial.txt' using 1:3 w lines title '' lt 3
#------------------------------------------------------------------------------------#

#----------------------------------GRAPH 3-------------------------------------------#
set title "Pumping profile"
set xlabel "x"
set ylabel "P(x)"
plot [][] 'initial.txt' using 1:4 w lines title '' lt 4
#------------------------------------------------------------------------------------#
unset multiplot
