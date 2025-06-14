set xlabel "X"
set ylabel "V"
set title "Lid driven cavity validation for Re=1000"
set datafile separator ","
set autoscale xy
set mouse
set grid
plot 'benchmark_v.csv' using 1:2 title "Ghia 1982 data" ,"result.csv" using 2:1 with lines lt rgb "blue" title "power law scheme results"



