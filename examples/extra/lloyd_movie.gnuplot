set term pngcairo size 800,800
set pointsize 0.6
set output 'FILENAME'
unset key
set title "Frame N"
plot [0:1] [0:1] 'lloyd/lloyd_p.N' u 2:3 w p pt 7 lt 4, 'lloyd/lloyd_v.N' lt 3
set output
