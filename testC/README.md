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
 	