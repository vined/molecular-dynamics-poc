unset log
unset label
set encoding utf8
set title "Pressure"
set xlabel "t, ps"
set ylabel "T, K"
plot "out/energies.dat" using 0:3 title 'Pressure' with lines
