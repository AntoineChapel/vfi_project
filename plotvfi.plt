set term png
set output "plot.png"
set title "Value Function Iteration"
set nokey
set grid
plot "datavfi.txt" using 1:2 with lines title "k'", "datavfi.txt" using 1:1 with lines title "45°"
