set title "Project1 - Question 2, Crank-Nicolson (dx = 0.025,dt = 0.00025)"
# Forward Euler
# Two-Stage TVD Runge-Kutta
# Backward Euler
# Crank-Nicolson
set xlabel "Spatial Position (x)"
set ylabel "Function Value (y)"
set xrange [-1:1]
set yrange [0:8]

set key right top
set grid
set mouse

colors = "#E50000 #98043F #4C087F #000DBF"
file_name = "Project1_22.dat"

set terminal png size 960,470
set output 'Project1_22.png'

plot file_name using 1:2 with lines lc rgb word(colors,1) title '0.0 sec',file_name using 1:3 with lines lc rgb word(colors,2) title '0.5 sec',file_name using 1:4 with lines lc rgb word(colors,3) title '1.0 sec',file_name using 1:5 with lines lc rgb word(colors,4) title '1.5 sec'
