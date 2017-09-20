set term pdfcairo
set output 'contour_p.pdf'
#set term png
#set output "countour_p.png"
# Set various features of the plot
#set view map
#set contour surface 
set key outside
#set pm3d
set surface
set mxtics 2
set mytics 2
#set cntrparam cubicspline  # smooth out the lines
#set pm3d interpolate 15,15 # interpolate the color (hoehere Zahl, schoenerer plot)

set isosample 400,400 
# Set a nice color palette
#set palette model RGB
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0) 
# Axes
#set zrange [0:0.7]
set xlabel '$z \; \left[ \textrm{fm} \right]$'
set ylabel '$t \; \left[ \textrm{fm} \right]$'
set zlabel '$p$ $\left[ \textrm{GeV}/\textrm{fm}^3 \right]$'
set ticslevel 0.0
set xrange [-0.3:0.3]
#set ztics 0.1, 0.1, 0.7
set view 40,200
#set view 50,300
# Now plot
splot './DATA/c_output_p.txt' using ($1-0.3):2:($3/1000) notitle with lines palette 
