
# Set the output file and the terminal type
set terminal png
set output 'richardson_alpha.png'

# Set the title, labels, and key
set title "Convergence History of relaxition Iterative Method"
set xlabel "Iteration"
set ylabel "Error"
set key outside

# Plot data from RESVEC.dat
plot 'RESVEC.dat' using 1 with lines title 'Convergence History'

# Export to png
set output
