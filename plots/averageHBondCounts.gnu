unset log
unset label
set encoding utf8
set title "Average count of hydrogen bonds"
set xlabel "Bonds"
set ylabel "Probability"
plot "out/temperature.dat" title 'Temperature' with lines
