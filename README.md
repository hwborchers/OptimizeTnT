# OptimizeTnT

> This repository will compile tips and tricks for solving optimization tasks of all kinds, continuous or discrete, constrained or unconstrained, smooth or nonsmooth, also special cases and rarely talked-about problem cases, but some comparisons between solvers available in R packages will be included.

### Tips 'n' Tricks

**[Linear equality constraints](LinearEqualityConstraints.md)**  
If there are *only* linear equality constraints (such as `sum(x) == 1` or `A x = b`) there is a 'trick' that will enable this problem to be solved without constraints at all.

**[Stereographic optimization](StereographicOptimization.md)**  
We will look at the *stereographic projection* and how it can help to solve optimization problems under the 'unit-length' constraint `sum(x^2) == 1`.

**[The Remez problem](RemezProblem.md)**  
The Remez problem is the task of finding or calculating the best polynomial approximation of a continuous function on a closed interval such that the maximum absolute distance between the polynomial and the function is minimized.|

**Selected optimization solvers**  
The Optimization Task View presents a long list of optimization solvers available in R packages. Here is a selection of modern and state-of-the-art solvers, taken from the task view, for different optimization problem classes.

**The 'Historize' routine**

**Compare Nelder-Mead Solvers**\
There are several implementations of Nelder-Mead optimization solvers in R packages. We will compare them in terms of accuracy and run-time behavior.

**Solving 'Minimax' optimization problems**

**Visualizing optimization solutions**

**The Bin Packing problem**  
Solve the *Bin Packing Problem* (BPP) as a 'Mixed-Integer Linear Programming Problem' and use this model to compare some MILP solvers available in packages on CRAN.

**Solve optimization problems with 'higher-accuracy' precision**  
If floating-point precision is not enough, there is the possibility to get more accurate optimization results by using 'big' floating-point numbers.


### Optimization Test Functions
