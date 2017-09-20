reset
set encoding utf8
set term gif animate size 1600, 600 font "arial"
set output "animate.gif"
n= 450   #n frames
set xrange [0:0.6]
i=0
t=0
load 'animate_gif.plt'
set output
