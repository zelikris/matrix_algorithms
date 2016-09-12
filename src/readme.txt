The program can be accessed through the use of flags as command line arguments. The first argument will be the path to the data file to be processed.

Format: project.jar ./data/1_H.dat -function -parameters

To run non-parametrized functions use

Example: java -jar project.jar ./data/1_H.dat -qr_fact_house

and the result will print to the screen. The functions that operate this way are:
(Using an nxn square matrix A)
-lu_fact
-qr_fact_house
-qr_fact_givens

(Using an an augmented matrix A|b)
-solve_lu
-solve_qr_house
-solve_qr_givens

For the paramaterized functions, each additional parameter begins with - and is in a specific order

Format: project.jar ./data/1_Hb.dat -jacobi_iter -"initial vector" -tolerance -maxIterations

The initial vector is enclosed by parentheses, and is a list of comma separated values.

Example: java -jar project.jar ./data/iterations.dat -jacobi_iter -"1.0, 0.0, 0.0" -0.0000001 -25

Two functions use this format:
-jacobi_iter
-gs_iter

The power method format is similar, but it takes another vector as an argument. It's use is

Format: project.jar ./data/power.dat -power_method -"initial vector" -"auxiliary vector" -tolerance -maxIterations

Example: java -jar project.jar ./data/power.dat -power_method -"1.0, 0.0" -"0.0, 1.0" -0.000001 -50