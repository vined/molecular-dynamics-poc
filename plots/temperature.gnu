unset log
unset label
set encoding utf8
set title "Temperature"
set xlabel "t, ps"
set ylabel "T, K"
plot "out/temperature.dat" title 'Temperature' with lines
