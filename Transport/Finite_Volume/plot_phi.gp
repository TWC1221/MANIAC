set term qt
set grid
set xlabel 'x'
set ylabel 'y'
set cblabel 'phi'
set colorbox vertical
set xrange [0.0000000:1.0000000]
set yrange [0.0000000:1.0000000]
set size ratio -1
set view map
set pm3d at b
set palette defined (0 '#000004', 0.25 '#3b0f70', 0.5 '#8c2981', 0.75 '#de4968', 1 '#fcfdbf')
splot 'phi_xy.dat' using 1:2:3 with pm3d notitle
# set contour base
# set cntrparam levels 10
# replot 'phi_xy.dat' using 1:2:3 with lines lc rgb 'black' notitle
# unset view; set view 60,30
# unset pm3d
# splot 'phi_xy.dat' using 1:2:3 with lines lw 1 notitle
