
# Set the output file and the terminal type
set terminal png
set output 'Gauss-siedle.png'

# Set the title, labels, and key
set title "Convergence History of Gauss-siedle Iterative Method"
set xlabel "Iteration"
set ylabel "Error"
set key outside

# Use logarithmic scale for the y-axis
set logscale y
#set logscale x
# Plot data from RESVEC.dat
plot 'RESVEC.dat' using 1 with lines title 'Convergence History'

# Export to png
set output
