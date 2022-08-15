# Smallest-ball-of-balls
Matlab code of a dual simplex-type algorithm for the smallest enclosing ball of balls.


## Introduction
This is the Matlab code of a dual simplex-type algorithm for computing the smallest enclosing ball of a set of balls and other closely related problems. Our algorithm employs a pivoting scheme resembling the simplex method for linear programming, in which a sequence of exact curve searches is performed until a new dual feasible solution with a strictly smaller objective function value is found. We utilize the Cholesky factorization and procedures for updating it, yielding a numerically stable implementation of the algorithm. Our algorithm can efficiently solve instances of dimension 5000 with 100000 points, often within minutes.

For further information on the algorithm, please consult the respective paper (see Reference / Citation).

## Running the code
To run the code, call function `main`:  function `[x, S, output] = main (P, options)`

### Input of main 
- P ~ matrix containing the data points - each column corresponds to a point;
- options ~ struct containing the options to apply to the algorithm.

Use `createinputstruct.m` to create the options struct. The options include:
- feasTol ~ Primal feasibility tolerance;
- feasOption ~ 1: it uses the most infeasible point; 2: it uses the first infeasible point found;
- maxTime ~ maximum running time allowed;
- iterLog ~ Whether to print a log of each iteration on the console during the algorithm's run.

### Output of main  
- x ~ optimal solution;
- S ~ basis (indexes of the points);
- output ~ a structure with output data.

The output struct includes:
- it ~ number of iterations
- bu ~ number of basis updates
- t ~ running time
- mbs ~ maximum size observed of a basis
- obs ~ size of the optimal basis
- status ~ Either MAXTIME or OPTIMAL

## Reference / Citation:
[Cavaleiro, M., Alizadeh, F. A dual simplex-type algorithm for the smallest enclosing ball of balls. Computational Optimization and Applications 79, 767-787 (2021)](https://link.springer.com/article/10.1007/s10589-021-00283-6)
