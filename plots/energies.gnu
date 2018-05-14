unset log
unset label
set encoding utf8
set title "Energy"
set xlabel "t, fs"
set ylabel "Energy, E"
plot "out/energies.dat" using 0:1 title 'Kinetic energy' with lines, '' using 0:2 title 'Potential energy' with lines lt 2, '' using 0:($1+$2) title 'Total energy' with lines lt 3
