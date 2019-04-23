set title "Project2 - Undergraduate"
# Forward Euler
# Two-Stage TVD Runge-Kutta
# Backward Euler
# Crank-Nicolson
set xlabel "Spatial Position (x)"
set ylabel "Function Value (y)"
set xrange [-8:8]
set yrange [-0.5:1]

set grid
set mouse

colors = "#E50000 #98043F #4C087F #000DBF"
file_name1 = "Project2_Ut_100.dat"
file_name2 = "Project2_Ut_200.dat"
file_name3 = "Project2_Ut_400.dat"
file_name4 = "Project2_Ut_800.dat"

set terminal png size 960,470
set output 'Project2_Ut_VelocityBoth.png'

#set multiplot layout 1,2 rowsfirst

#set key left top
#set title "Pressure"
#plot file_name1 using 1:5 with lines lc rgb word(colors,1) title '100 intervals',file_name2 #using 1:5 with lines lc rgb word(colors,2) title '200 intervals',file_name3 using 1:5 with lines #lc rgb word(colors,3) title '400 intervals',file_name4 using 1:5 with lines lc rgb #word(colors,4) title '800 intervals'
#unset key

set key left top
set title "Velocity"
plot file_name1 using 1:9 with lines lc rgb word(colors,1) title '100 intervals',file_name2 using 1:9 with lines lc rgb word(colors,2) title '200 intervals',file_name3 using 1:9 with lines lc rgb word(colors,3) title '400 intervals',file_name4 using 1:9 with lines lc rgb word(colors,4) title '800 intervals'
unset key

#unset multiplot