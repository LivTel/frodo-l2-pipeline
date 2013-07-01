#!/usr/bin/gnuplot

set term png size 640,480
set output "plt.png"
#set terminal x11 persist

set xdata time
set timefmt "%Y-%m-%d	%H:%M%S"
set ylabel "average FWHM (px)"
set key top left
set yrange [2:3.5]
set xrange ["2013-01-01	00:00:00":"2013-07-01	00:00:00"]

set tmargin 0
set bmargin 2
set lmargin 10
set rmargin 3
unset xtics

set multiplot layout 2,1 title ""

# can use "using 2:17:18 with yerrorbars" to get errorbars (need notitle too)
plot "Red-Low" using 2:17:18 with lines lw 2 linecolor rgb "red" title "Red-Low", \
"Blue-Low" using 2:17:18 with lines lw 2 linecolor rgb "blue" title "Blue-Low"

set xtics nomirror
set tics scale 0 font ",8"
set xlabel "date"

plot "Red-High" using 2:17:18 with lines lw 2 linecolor rgb "red" title "Red-High", \
"Blue-High" using 2:17:18 with lines lw 2 linecolor rgb "blue" title "Blue-High"

unset multiplot

