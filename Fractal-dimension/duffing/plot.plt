set terminal png size 640,480
set output "duffing.png"
unset key
#set size square
set xlabel "x"
set ylabel "p"
plot "duffing.dat" u 1:2 w p ps 0.2 pt 7 lc rgb 'blue'
unset output

