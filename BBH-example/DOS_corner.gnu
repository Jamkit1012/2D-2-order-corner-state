set encoding iso_8859_1
set terminal postscript enhanced color
set output 'corner.eps'
set palette cubehelix start 3 cycles 0.5 saturation 1 negative 
set size square
set xrange [1 : 20]
set yrange [1 :20]
set pm3d
set view equal xyz
set view map
set border lw 1
unset key
unset colorbox
set pm3d interpolate 10,10
splot "corner.dat" u 1:2:($3) w pm3d
