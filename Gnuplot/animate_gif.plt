t=t+0.002

set multiplot layout 1, 2 title sprintf('t=%f fm/c', t)

#set bmargin 2
#set tmargin 1 
set yrange[0:15000]
set xlabel 'x (fm)'
set ylabel 'p [MeV/fm^3]'
plot sprintf('IDEAL/output_%i.txt',i) using 1:4 title "IDEAL" lc rgb 'blue'


set ylabel 'v_x [c]'
set yrange[-1:1]
plot sprintf('IDEAL/output_%i.txt',i) using 1:3 title "IDEAL" lc rgb 'blue'

unset multiplot

 i=i+1
 if (i <n+1) reread
