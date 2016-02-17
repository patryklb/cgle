#!/usr/bin/gnuplot -persist
#set terminal png size 2000,1000 enhanced font "Helvetica,20"
#set output 'cp/3.png'

set term png size 1024, 720  enhanced font "Helvetica"
set output "module.png"
set xlabel "t"
set ylabel "x"

#set palette
#set pm3d  map
#set hidden3d 
#set pm3d map
#set view 0, 0
unset surf
set pm3d at b
set view 0,0

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')

splot [][] 'module.txt' using 1:2:3 with lines
