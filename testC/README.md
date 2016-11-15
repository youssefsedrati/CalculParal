Calcul Parallèle
===================

Compilation
-------------
* Use the command
```
make
```

Utilisation
-------------------
* Set the parameters
	* Modify the file `param.dat`. This file is written as follows
		* First line :   `Nx Ny`
		* Second line :  `Lx Ly D`
		* Third line :  `method max_iter eps`
			* `method` indicates the computation method : 
				* To use Jacobi, set this value as 0. 
				*  To use Conjugate Gradient, set this value as 1
			* `max_iter` indicates the maximum number of iterations
			* `eps` indicates the error value that is used to determine the convergence

* Launch 
```
mpiexec -np <nb_processes> ./main.out
```

* Plot the solution

	* After launching the program using `<nb_processes>` processes, the following command with gather the results in one datafile and plot the computed solution
	```
	./plot.rb <nb_processes>
	``` 
	*	Change the variable `plotfile` in the file `plot.rb` to change the reference curve (that corresponds to the real solution expected)
		*	Set this variable to `gnuplot.cmd`for `sin(x) + cos(y)`
		*	Set this variable to `gnuplot2.cmd`for `x*(1-x)*y*(1-y)`  

* Plot speedups
 	* The following command will compute the speedup curves using the set parameters, with number of processes starting from 1 up to `<max_processes>`
	```
	./benchmarks.rb <max_processes>
	``` 
	*	This will also generate solution files that can be plotted using gnuplot