set term pdf
set size ratio -1
set pointsize 0.8
set key top left
set output 'out1.pdf'
plot [150:600] [200:750] 'out1.gnu' w l t 'Voronoi cells', 'out1.bd' u 2:3 w l lt 3 lw 3 t 'Boundary', 'out1.pts' u 2:3 w p pt 7 lt 4 t 'Generators'
set output 'out2.pdf'
set key top right
plot [125:640] 'out2.gnu' w l t 'Voronoi cells', 'out2.bd' u 2:3 w l lt 3 lw 3 t 'Boundary', 'out2.pts' u 2:3 w p pt 7 lt 4 t 'Generators'
set output 'out3.pdf'
set key top left
plot [0:475] [650:920] 'out3.gnu' w l t 'Voronoi cells', 'out3.bd' u 2:3 w l lt 3 lw 3 t 'Boundary', 'out3.pts' u 2:3 w p pt 7 lt 4 t 'Generators'
set output 'out5.pdf'
set key bottom right
plot [0:425] [600:1050] 'out5.gnu' w l t 'Voronoi cells', 'out5.bd' u 2:3 w l lt 3 lw 3 t 'Boundary', 'out5.pts' u 2:3 w p pt 7 lt 4 t 'Generators'
set output
