#!/usr/bin/env gnuplot

set datafile separator ","
set term pdf
set output 'speedup.pdf'

set xlabel 'size of image'
set ylabel 'speed-up'
set xtics
set ytics
set grid

set title 'Hough Circle Detection Speed-up'

plot './data/eye0.csv' using 1:4 title 'eye0.jpg' with linespoints, \
     './data/eye1.csv' using 1:4 title 'eye1.jpg' with linespoints, \
     './data/eye2.csv' using 1:4 title 'eye2.jpg' with linespoints

set term x11
