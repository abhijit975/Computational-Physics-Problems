set term png size 640,480
set output "fitdata.png"
set size square
unset key
set xlabel "ln(b)"
set ylabel "ln[N(b)]"
set title "Fractal dimension D_f = 1.40"
f(x) = -1.40355*x+2.131 
plot "dimension.dat" u 1:2 w p lc rgb 'red', f(x) w l lc rgb 'blue'
unset output
