#!/usr/bin/env gnuplot

set datafile separator ","
set term pdf
set output 'comparison.pdf'

set xlabel 'size of image'
set ylabel 'speed-up'
set xtics
set ytics
set grid

set title 'Hough Circle Detection Original vs. Optimized Implementation'

plot './data/eye0.csv' using 1:2 title 'original:  eye0.jpg' with linespoints, \
     './data/eye1.csv' using 1:2 title 'original:  eye1.jpg' with linespoints, \
     './data/eye2.csv' using 1:2 title 'original:  eye2.jpg' with linespoints, \
     './data/eye0.csv' using 1:3 title 'optimized: eye0.jpg' with linespoints, \
     './data/eye1.csv' using 1:3 title 'optimized: eye1.jpg' with linespoints, \
     './data/eye2.csv' using 1:3 title 'optimized: eye2.jpg' with linespoints

set term x11
