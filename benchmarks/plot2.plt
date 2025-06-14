set xlabel "Y"
set ylabel "U"
set title "Lid driven cavity validation for Re=1000"
set datafile separator ","
set autoscale xy
set grid
plot 'benchmark_u.csv' using 1:2 title "Ghia 1982 data" ,"result2.csv" using 2:1 with lines lt rgb "blue" title "power law scheme results"

