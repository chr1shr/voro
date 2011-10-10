# This gnuplot script will generate EPS figures using the code output for the
# four example files
set term postscript eps color solid
set xlabel 'x'
set ylabel 'y'
set output'input1.eps'
set key top left
plot 'input1.gnu' lc rgbcolor "#0088ff" t "Voronoi cells", 'input1.bd' lw 4 lt 3 t "Boundary", 'input1.net' lt 4 t "Connections", 'input1' u 2:3 w p pt 7 lt 1 t "Generators"
set output 'input2.eps'
set key top right
plot [*:600] 'input2.gnu' lc rgbcolor "#0088ff" t "Voronoi cells", 'input2.bd' lw 4 lt 3 t "Boundary", 'input2.net' lt 4 t "Connections", 'input2' u 2:3 w p pt 7 lt 1 t "Generators"
set output 'input3.eps'
plot 'input3.gnu' lc rgbcolor "#0088ff" t "Voronoi cells", 'input3.bd' lw 4 lt 3 t "Boundary", 'input3.net' lt 4 t "Connections", 'input3' u 2:3 w p pt 7 lt 1 t "Generators"
set output 'input4.eps'
set key bottom
plot 'input4.gnu' lc rgbcolor "#0088ff" t "Voronoi cells", 'input4.bd' lw 4 lt 3 t "Boundary", 'input4.net' lt 4 t "Connections", 'input4' u 2:3 w p pt 7 lt 1 t "Generators"
set output
