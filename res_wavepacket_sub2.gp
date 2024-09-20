#!/usr/bin/gnuplot
set term postscript enhanced color size 30cm,15cm font 'Helvetica,26' lw 4
set output "gp_outfile.eps"

set bmargin 0.5
set lmargin 5.0
set rmargin -5.0
set cbtics font ',20'

set xlabel "R (a.u.)"
set ylabel "t (s)"
set zlabel "P (a.u.)"

set xrange [5.8:6.8]
set yrange [-1.2e-15:1.e-13]

unset key

set view map
set size ratio 0.5 0.8,1

splot ARG1 u 1:2:4 w pm3d
