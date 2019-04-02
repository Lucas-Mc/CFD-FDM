set title "Project2 - Undergraduate"
# Forward Euler
# Two-Stage TVD Runge-Kutta
# Backward Euler
# Crank-Nicolson
set xlabel "Spatial Position (x)"
set ylabel "Function Value (y)"
set xrange [-8:8]
set yrange [0:2]

set grid
set mouse

colors = "#E50000 #98043F #4C087F #000DBF"
file_name = "Project2_Ut_test.dat"

set terminal png size 960,470
set output 'Project2_Ut_test.png'

set multiplot layout 1,2 rowsfirst

set key left top
set title "Pressure"
plot file_name using 1:2 with lines lc rgb word(colors,1) title '0.0 sec',file_name using 1:3 with lines lc rgb word(colors,2) title '0.5 sec',file_name using 1:4 with lines lc rgb word(colors,3) title '1.0 sec',file_name using 1:5 with lines lc rgb word(colors,4) title '2.0 sec'
unset key

set key left top
set title "Velocity"
plot file_name using 1:6 with lines lc rgb word(colors,1) title '0.0 sec',file_name using 1:7 with lines lc rgb word(colors,2) title '0.5 sec',file_name using 1:8 with lines lc rgb word(colors,3) title '1.0 sec',file_name using 1:9 with lines lc rgb word(colors,4) title '2.0 sec'
unset key

unset multiplot