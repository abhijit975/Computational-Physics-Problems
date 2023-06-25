Q1.

To create the executable file please type cc *.c -lm -o q1.x
Executing q1.x creates a data file duffing.dat which can be plotted using plot.plt file
by writing the command in terminal gnuplot plot.plt. It creates the duffing.png file with
the image of the strange attractor.

Then copy the duffing.dat file into the dimension folder. Compiling and executing the
dimension.c file creates a data file data.dat which can be fitted in gnuplot using log(b)
and log(N(b)). Linear fit yields slope = -1.40. 

Q2.

To create the executable file please type cc *.c -lm -o dla.x
Executing the dla.x file creates two data files data.dat and gyration.dat
Plotting the data.dat file yields the dla plot.
The gyration.dat file lists Log(Rg) and Log(N) which when fitted yields a slope of 1.84.


