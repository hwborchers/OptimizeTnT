
# Recommended Optimization Solvers

The Optimization Task View presents a long list of optimization solvers available in R packages. Here is a selection of modern and state-of-the-art solvers, taken from the task view, for different optimization problem classes. Each solver is accompanied by a usage line, a link to its help page, and maybe to an appropriate Wikipedia page.

**Contents**

* [Unconstrained optimization](#unconstrained-optimization)
* [Constrained optimization](#constrained-optimization)
* [Quadratic optimization](#quadratic-optimization)
* [Least-squares problems](#least-squares-problems)
* [Convex optimization](#convex-optimization)
* [Derivative-free optimization](#derivative-free-optimization)
* [Global and stochastic optimization](#global-and-stochastic-optimization)
* [Linear and mixed-integer programming](#linear-and-mixed-integer-programming)
* [Multivariate root finding](#multivariate-root-finding)
* [Appendix I: Base R solvers](#appendix-i-r-solvers)
* [Appendix II: External Solvers](#appendix-ii-external-solvers)
* [References](#references)


## Continuous Optimization

*All optimization solvers in Base R have been moved to Appendix I. This page will focus on more modern and state-of-the-art solvers available in packages on CRAN or Github.*

----

### Unconstrained optimization

> In *numerical optimization* the goal is to minimize (or maximize) an objective function that depends on real variables: $\min_x f(x)$. *Unconstrained* means there are no restrictions on these variables, except possibly box (aka bounds) constraints $l_i \le x_i \le u_i$ for all $i$.

#### ucminf::**ucminf**

```r
ucminf(par, fn, gr = NULL, ..., control = list(), hessian = 0)
```
[help page](https://cran.r-universe.dev/ucminf/doc/manual.html#ucminf)

An algorithm for general-purpose unconstrained non-linear optimization. The algorithm is of quasi-Newton type with BFGS updating and soft line searching. It is effective and accurate.

Note: The interface of `ucminf` is designed for easy interchange with ‘optim’.

----

#### lbfgs::**lbfgs**

```r
lbfgs(call_eval, call_grad, vars, ...,
      max_iterations = 0, ftol = 1e-4,
      orthantwise_c = 0                 # > 0 to start the OWL-QN algorithm
)
```
[help page](https://cran.r-universe.dev/lbfgs/doc/manual.html#lbfgs) and
[libLBFGS](http://www.chokkan.org/software/liblbfgs/index.html) page.

`lbfgs` performs function optimization using the L-BFGS and *Orthant-Wise Limited-memory Quasi-Newton* optimization (OWL-QN) algorithms. It's a wrapper to the *libLBFGS* library by Naoaki Okazaki, based on an implementation of the L-BFGS method written by Jorge Nocedal.

*Note*: There are many more options to set if needed, see the help page. It has no box constraints. The OWL-QN algorithm is particularly suited for higher-dimensional problems (but can be extremely slow).

----

#### nloptr::**lbfgs**

```r
lbfgs(x0, fn, gr = NULL, lower = NULL, upper = NULL,
      nl.info = FALSE, control = list(), ...)

varmetric(x0, fn, gr = NULL, rank2 = TRUE,  # rank-2|-1 update method
          lower = NULL, upper = NULL,
          nl.info = FALSE, control = list(), ...)

tnewton(x0, fn, gr = NULL, lower = NULL, upper = NULL,
        precond = TRUE, restart = TRUE,
        nl.info = FALSE, control = list(), ...)
```
lbfgs [help page](https://cran.r-universe.dev/nloptr/doc/manual.html#lbfgs) |
varmetric [help page](https://cran.r-universe.dev/nloptr/doc/manual.html#varmetric) |
tnewton [help page](https://cran.r-universe.dev/nloptr/doc/manual.html#tnewton)

Package 'nloptr' supports the following quasi-Newton methods:

* `lbfgs`: Low-storage version of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) method. 
* `varmetric`: Shifted limited-memory variable-metric algorithm.
* `tnewton`: Truncated Newton method, using a conjugate-gradient approach.

Relevant options are `xto_rel = 1e-6`, `maxeval = 1000`, but also
`stopval = -Inf` or `check_derivatives = FALSE`,
see `nloptr.print.options()` for all options!

*Note*: Package 'nloptr' is a wraper for [NLopt](https://nlopt.readthedocs.io/en/latest/), an open-source library for nonlinear optimization.

----

----

### Constrained optimization

> There are many different types of constraints. Here we list solvers that are capable of handling general non-linear constraints, equality and inequality constraints. The preferred method for doing this is the 'augmented Lagrangian' approach.
-- [wp:Augmented_Lagrangian_method](https://en.wikipedia.org/wiki/Augmented_Lagrangian_method)


#### alabama::**auglag**

```r
auglag(par, fn, gr, hin, hin.jac, heq, heq.jac, 
       control.outer=list(), control.optim = list(), ...
)
```
[help page](https://cran.r-universe.dev/alabama/doc/manual.html#auglag)

`auglag 'implements the Augmented Lagrangian Minimization Algorithm for optimizing smooth nonlinear objective functions with constraints. Linear or nonlinear equality and inequality constraints are allowed.

*Note*: There is also function `constrOptim.nl`, applying an Augmented Lagrangian Adaptive Barrier Minimization Algorithm; can be used to replace `constrOptim` in Base R.

----

#### NlcOptim::**solnl**

```r
solnl(X = NULL, objfun = NULL, confun = NULL, A = NULL, B = NULL,
      Aeq = NULL, Beq = NULL, lb = NULL, ub = NULL, tolX = 1e-05,
      tolFun = 1e-06, tolCon = 1e-06, maxnFun = 1e+07, maxIter = 4000
)
```
[help page](https://cran.r-universe.dev/NlcOptim/doc/manual.html#solnl)

`solnl` implements the Sequential Quadratic Programming (SQP) method to find solutions for general nonlinear optimization problems (with nonlinear objective and constraint functions). Linear or nonlinear equality and inequality constraints are allowed. 

*Note*: The SQP method is described in detail in the book "Numeric Optimization" by Nocedal and Wright. [wp:SQP](https://en.wikipedia.org/wiki/Sequential_quadratic_programming)

----

#### BB::**spg**

```r
spg(par, fn, gr=NULL, method=3, lower=-Inf, upper=Inf, 
    project=NULL, projectArgs=NULL, 
	control=list(), quiet=FALSE, alertConvergence=TRUE, ...
)
```
[help page](https://cran.r-universe.dev/BB/doc/manual.html#spg)

`spg` applies Barzlai-Borwein spectral methods for large-scale nonlinear optimization subject to simple constraints. Based on code from the [TANGO](https://www.ime.usp.br/~egbirgin/tango/) project.
[wp:BB](https://en.wikipedia.org/wiki/Barzilai-Borwein_method)

*Note*: Use the projection function `projectLinear` for equality and inequality constraints `A x - b >= 0`. 

----

#### nloptr::**slsqp**

```r
slsqp(x0, fn, gr = NULL, lower = NULL, upper = NULL,
      hin = NULL, hinjac = NULL, heq = NULL, heqjac = NULL,
      nl.info = FALSE, control = list(), ...
)

cobyla(x0, fn, lower = NULL, upper = NULL, hin = NULL,
       nl.info = FALSE, control = list(), ...
)
```
`slsqp` [help page](https://cran.r-universe.dev/nloptr/doc/manual.html#slsqp),
`cobyla` [help page](https://cran.r-universe.dev/nloptr/doc/manual.html#cobyla), and
[NLopt](https://nlopt.readthedocs.io/en/latest/) Documentation

`slsqp` implements a Sequential (least-squares) quadratic programming (SQP) algorithm for nonlinearly constrained, gradient-based optimization, supporting both equality and inequality constraints.

`cobyla` implements an algorithm for derivative-free optimization with nonlinear inequality constraints (no equality constraints). It constructs successive linear approximations  and constraints via simplices.

----

----

### Quadratic optimization

> Quadratic optimization (or: Quadratic Programming, QP) problems are optimization problems of the form $\min_x 0.5\, x'Px + q' x$ with linear or quadratic constraints. These problems are ubiquitous and appear in almost all data fitting and least-squares approximation tasks.
-- [wp:Quadratic programming](https://en.wikipedia.org/wiki/Quadratic_programming)

#### osqp::**solve_osqp**

```r
solve_osqp(P = NULL, q = NULL,
           A = NULL, l = NULL, u = NULL,
           pars = osqpSettings()
)
```
[help page](https://cran.r-universe.dev/osqp/doc/manual.html#solve_osqp) and [OSQP](https://osqp.org/) home page.

`solve_osqp` solves Quadratic Programming (QP) problems with linear inequality constraints, using the OSQP library of the University of Oxford. `P` and `q` as above, linear inequalities as `l <= Ax <= u`, `A` a (sparse) matrix.

----

#### piqp::**solve_piqp**

```r
solve_piqp(P = NULL, c = NULL,
           A = NULL, b = NULL, G = NULL, h = NULL,
           x_lb = NULL, x_ub = NULL,
           settings = list(), backend = c("auto", "sparse", "dense")
)
```
[help page](https://cran.r-universe.dev/piqp/doc/manual.html#solve_piqp) and [PIQP](https://predict-epfl.github.io/piqp/) home page.

`solve_piqp` solves Quadratic Programming (QP) problems with linear equality *and* inequality constraints, using the very fast "Proximal Interior Point Quadratic Programming" solver of the EPFL. `c` is `q` above, equality constraints `A x = b`, inequality constraints `G x <= h` and box constraints `x_lb <= x <= x_ub`.

----

----

### Least-squares problems

> Least-squares is a standard approach in regression analysis to approximate the solution of linear or nonlinear systems. Nonlinear least-squares problems are usually solved by iterative procedures, quite similar to optimization problems.
-- [wp:Least squares](https://en.wikipedia.org/wiki/Least_squares)

#### nlsr::**nlsr**

```r
nlsr(formula = NULL, data = NULL, start = NULL,
     control = NULL, trace = FALSE, subset = NULL,
     lower = -Inf, upper = Inf,  weights = NULL, ...
)
```
[help page](https://cran.r-universe.dev/nlsr/doc/manual.html#nlsr)

`nlsr` provides solutions to a nonlinear least squares problem using the Nash Marquardt tools, i.e., internally applying the Levenberg-Marquardt algorithm. Conceived as a replacement for `nls`.

----

#### minpack.lm::**nls.lm**

```r
nls.lm(par, lower=NULL, upper=NULL, fn, jac = NULL,
       control = nls.lm.control(), ...
)
```
[help page](https://cran.r-universe.dev/minpack.lm/doc/manual.html#nls.lm)

`nls.lm` minimizes the sum-of-squares of the vector returned by the function fn, by a modification of the Levenberg-Marquardt algorithm. The user may also provide a function `jac` which calculates the Jacobian.

----

----

### Convex optimization

> A convex optimization problem is one where the objective function and all constraints are convex. There are very effective interior-point methods to solve convex problems quickly and with high accuracy. Convex problems have only one local optimum, i.e., a local solution is also a global solution.
-- [wp:Convex optimization](https://en.wikipedia.org/wiki/Convex_optimization)

#### CVXR:**solve**

```r
n <- length(data)
x <- Variable(n)                            # initialize n variables
objective <- Minimize(p_norm(data - x, 2))  # define the objective function
constraint <- list(diff(x) >= 0)            # and the list of constraints
problem <- Problem(objective, constraint)   # check the problem for convexity
result <- solve(problem)                    # and send it to a convex solver
result$getValue(x)                          # retrieve the results
```
[help page](https://cran.r-universe.dev/CVXR/doc/manual.html#solve) and [CVXR](https://cvxr.rbind.io/) home page.

'CVXR' is an R package that implements the disciplined convex programming (DCP) approach by Steven Boyd and colleagues to verify the problem’s convexity. It provides an object-oriented modeling language that allows the user to formulate convex optimization problems in a natural mathematical syntax. When the model has been verified, it will be sent to a convex solver such as ECOS or SCS.

*Note*: The code above is an example of how the modeling language looks like. It solves an isotonic least-squares regression $x_1 \le x_2 \le \ldots \le x_n$. There is no single `solve` function whose usage could be described. See the [Tutorial Examples](https://cvxr.rbind.io/examples/) in the documentation.

*Remark*: There are more convex solvers in R packages (cf. 'clarabel', 'sdpt3r', 'cccp', etc.), and their usage can indeed be a quite complex task even for an experienced user.

----

----

### Derivative-free optimization

> 'Derivative-free' means numerical optimization without the use of gradients (or hessians). The most popular candidate is Nelder-Mead, a simplex-based direct search method.
-- [wp:Nelder-Mead](https://en.wikipedia.org/wiki/Nelder-Mead_method) and [wp:Derivative-free optimization](https://en.wikipedia.org/wiki/Derivative-free_optimization)

> *Note*: Gradient-free solvers without constraints can be equipped with inequality constraints by wrapping them into a function that returns `Inf` (or some large value) if constraints are violated.

> *Remark*: Saying an algorithm is derivative-free does *not* mean it will minimize nonsmooth objective functions effectively!

#### dfoptim::**nmk**

```r
nmk(par, fn, control = list(), ...)

nmkb(par, fn, lower=-Inf, upper=Inf, control = list(), ...)
```
[help page](https://cran.r-universe.dev/dfoptim/doc/manual.html#nmk)

`nmk` is an implementation of Nelder-Mead based on MATLAB code by C.T. Kelley for his book "Iterative Methods for Optimization". The function `nmkb` allows box constraints for Nelder-Mead.

----

#### nloptr::**neldermead**

```r
neldermead(x0, fn, lower = NULL, upper = NULL,
           nl.info = FALSE, control = list(), ...
)
```
[help page](https://cran.r-universe.dev/nloptr/doc/manual.html#neldermead)

An implementation of almost the original Nelder-mead simplex algorithm, but provides explicit support for box constraints. The function is quite fast and accurate.

----

#### pracma::**anms**

```r
anms(fn, x0, ..., tol = 1e-10, maxfeval = NULL)
```
[help page](https://cran.r-universe.dev/pracma/doc/manual.html#anms)

`anms` is an implementation of Nelder-Mead with adaptive parameters, based on the ideas of Gao & Han 2012. It is particularly suitable for higher-dimensional problems (20-30 variables).

----

----

## Global and Stochastic Optimization

### Differential Evolution

> Differential Evolution (DE) optimizes a problem by maintaining a population of candidate solutions and creating new candidate solutions by combining existing ones according to their simple formulae, and then keeping whichever candidate solution has the best score or fitness.
-- [wp:DE](https://en.wikipedia.org/wiki/Differential_evolution) and [wp:Global optimization](https://en.wikipedia.org/wiki/Global_optimization)

#### DEoptim::**DEoptim**

```r
DEoptim(fn, lower, upper,
        control = DEoptim.control(), ..., fnMap=NULL
)
```
[help page](https://cran.r-universe.dev/DEoptim/doc/manual.html#DEoptim)

`DEoptim` implements the Differential Evolution (DE) algorithm similar to what is described in the book "Differential Evolution – A Practical Approach to Global Optimization" by Price et al. The package relies on an interface to a C implementation of DE, which is effective and fast.

The control options `NP` (for population size) and `itermax` should probably be set 5-10 times higher than the default. The recommended strategies are `2` (default) or `3` in `DEoptim.control(strategy = ..)`.

----

#### DEoptimR::**JDEoptim**

```r
JDEoptim(lower, upper, fn,
         constr = NULL, meq = 0,
         NP = 10*length(lower),
         maxiter = 200*length(lower), ...
)
```
[help page](https://cran.r-universe.dev/DEoptimR/doc/manual.html#JDEoptim)

`JDEoptim` implements in pure R a Differential Evolution (DE) algorithm with "self-adapting control parameters", as suggested by Brest et al. 2006 (with Java code). It incorporates very general nonlinear (equality and inequality) constraints (though equality constraints seem to be problematic).

*Note*: It has many more control options than are shown on the usage line. The default values for population size `NP` and `maxiter` are too small.

----

### Simulated Annealing

Simulated Annealing is a stochastic search technique for global optimization. It is often used for functions with very many local minima. The algorithm is inspired by the annealing process in metallurgy.  
-- [wp:Simulated Annealing](https://en.wikipedia.org/wiki/Simulated_annealing)

#### GenSA::**GenSA**

```r
GenSA(par = NULL, fn, lower, upper, control = list(), ...)
```
[help page](https://cran.r-universe.dev/GenSA/doc/manual.html#GenSA)

`GenSA` searches for a global minimum of complicated non-linear objective functions with a large number of optima. It implements a generalized "Simulated Annealing" (SA) approach in C++, that is effective though not always as successful as Differential Evolution.

----

### Genetic Algorithms

----

### Particle-swarm Optimization

> Particle-swarm Optimization (PSO) is a heuristic approach that optimizes a problem by iteratively trying to improve a set of particles (candidate solutions) with regard to a given measure of quality. The original idea was to simulate social interactions in a swarm of, e.g., fish or birds when looking for food.
-- [wp:PSO](https://en.wikipedia.org/wiki/Particle_swarm_optimization)

#### pso::**psoptim**

```r
psoptim(par, fn, gr = NULL, ..., lower = -1, upper = 1,
        control = list()
)
```
[help page](https://cran.r-universe.dev/pso/doc/manual.html#psoptim)

`psoptim` is a general implementation of particle swarm optimization usable as a direct replacement for optim. Control elements include the maximum number of iterations or the swarm size.

Note: This implementation of "Particle Swarm Optimization" is *not* derivative-free!

----

### CMA-ES

> Covariance matrix adaptation evolution strategy (CMA-ES) is a special of evolution strategy for numerical optimization of non-linear or non-convex continuous optimization problems. For more detailed information see
-- [wp:CMA-ES](https://en.wikipedia.org/wiki/CMA-ES)

#### rCMA::**cmaOptimDP**

```r
cmaOptimDP(cma, fitFunc, isFeasible = function(x) {     TRUE },
  maxDimPrint = 5, iterPrint = 10, verbose = 2)
```
[help page](https://cran.r-universe.dev/rCMA/doc/manual.html#cmaOptimDP)

*Note*: `cmaOptimDP` utilizes Java code by Hansen. The package is therefore depending on 'rJava' that on some systems may lead to installation problems.

----

#### adagio::**pureCMAES**

```r
pureCMAES(par, fun, lower = NULL, upper = NULL, sigma = 0.5,
          stopfitness = -Inf, stopeval = 1000*length(par)^2, ...)
```
[help page](https://cran.r-universe.dev/pracma/doc/manual.html#pureCMAES)

The Octave/MATLAB function `pureCMAES` has been converted to R. This function will be used for higher-dimensional search spaces (up to 30-50 variables), though it can become very slow. The results are in general quite good.

*Note*: There are more implementations of CMA-ES, see the CRAN [Optimization Task View](https://CRAN.R-project.org/view=Optimization).

----

----

## Linear and Mixed-integer Programming

### Linear Programming

> linear programming is a technique for the optimization of a linear objective function, subject to linear equality and inequality constraints.
-- [wp:Linear Programming](https://en.wikipedia.org/wiki/Linear_programming)

#### highs::**highs_solve**

```r
highs_solve(Q = NULL, L, lower, upper,
            A, lhs, rhs, types,
            maximum = FALSE, offset = 0, control = highs_control()
)
```
[help page](https://cran.r-universe.dev/highs/doc/manual.html#highs_solve)
and the [HiGHS](https://highs.dev/) Web page.

'highs' wraps the HiGHS solver, currently among the best open-source mixed integer linear programming solvers. Solves linear and quadratic problems of the form $0.5 x'Qx + Lx$ with linear and bounds constraints:
$\text{lhs} \le Ax \le \text{rhs}$, $\text{lower} \le x \le \text{upper}$, and different types of variables (1 or 'C' continuous, 2 or 'I' integer, 3 or 'SC' semi-continuous, etc.).

*Note*: Quadratic problems do not allow for integer variables.

----

#### rcbc::**cbc_solve**

```r
cbc_solve(obj, mat,
  row_ub = rep.int(Inf, nrow(mat)),       # nrow(mat) no. of constraints
  row_lb = rep.int(-Inf, nrow(mat)),
  col_lb = rep.int(-Inf, ncol(mat)),      # ncol(mat) no. of variables
  col_ub = rep.int(Inf, ncol(mat)),
  is_integer = rep.int(FALSE, ncol(mat)), # type of variables
  max = FALSE,
  cbc_args = list(),
  initial_solution = NULL
)
```
[help page](https://rdrr.io/github/yuehmeir2/myFormAssembler/man/cbc.html)

CLP and CBC are open-source [COIN-OR](https://www.coin-or.org/) projects. CLP solves linear programs and CBC adds handling of integer variables. 'rcbc' provides an interface for mixed-integer linear programs (without integer variables CLP will be called).

*Remark*: 'rcbc' has to be installed from
`https://dirkschumacher.github.io/rcbc/`.

----

----

## Multivariate Root Finding

> Finding roots of multivariate functions could in principle be done by minimizing the square of the function values. However there are more specialized procedures to perform this task, inspired by Gauss-Newton or spectral approaches.

#### nleqslv::**nleqslv**

```r
nleqslv(x, fn, jac=NULL, ...,
    method = c("Broyden", "Newton"),
    global = c("dbldog", "pwldog", "cline", "qline", "gline", "hook", none"),
    xscalm = c("fixed","auto"),
    jacobian=FALSE, control = list()
)
```
[help page](https://cran.r-universe.dev/nleqslv/doc/manual.html#nleqslv)

Solves a system of nonlinear equations using a Broyden or a Newton method with a choice of global strategies such as line search and trust region. `fn` must be a function $f: \mathbb{R}^n \to \mathbb{R}^n$.

----

#### BB::**sane**

```r
sane(par, fn, method=2, control=list(),
       quiet=FALSE, alertConvergence=TRUE, ...
)

dfsane(par, fn, method=2, control=list(),
       quiet=FALSE, alertConvergence=TRUE, ...
)
```
`sane` [help page](https://cran.r-universe.dev/BB/doc/manual.html#sane) and
`dfsane` [help page](https://cran.r-universe.dev/BB/doc/manual.html#dfsane)

`sane` provides a non-monotone spectral approach for solving large-scale nonlinear systems of equations, and `dfsane` provides a derivative-free version of this spectral approach. The function value must have the same length as the input vector.

----

----

## Appendix I: R Solvers

### Base R Solvers

> R provides `optimize` (univariate optimization) and `optim` or `nlminb` (multivariate optimization) for function minimization. Besides that, there is `nlm` which in general is not recommended. 

#### **optimize**

```r
optimize(f, interval, ..., lower = min(interval), upper = max(interval),
         maximum = FALSE, tol = .Machine$double.eps^0.25
)
```
[help page](./man/optimize.html)

Univariate (i.e., onedimensional) optimization aims to find the optimum of a continuous, real-valued function in an interval of the real numbers (with respect to its first argument). `optimize` uses a combination of golden section search and successive parabolic interpolation.

*Note*: When there is more than one local minimum, it is uncertain which one will be chosen.

----

#### **optim**

```r
optim(par, fn, gr = NULL, ..., method = c("Nelder-Mead"),
      lower = -Inf, upper = Inf,
      control = list(), hessian = FALSE
)
```
[help page](./man/optim.html)

General-purpose multivariate optimization; different methods, Nelder-Mead, BFGS, and SANN are provided.,The gradient function if needed, can be user-supplied or will be calculated with finite differences. It includes options for handling box constraints.

The `optim` function supports the following methods:

* `optim(method="Nelder-Mead")` -- Nelder-Mead, *derivative-free* and 
   robust, but relatively slow. 
* `optim(method="BFGS")` -- 'BFGS' implementation, a quasi-Newton method. 
* `optim(method="L-BFGS-B")` -- 'BFGS' method, allows for box constraints
* `optim(method="CG")` conjugate gradients method, *not* recommended.
* `optim(method="SANN")` variant SA method, *not* recommended.
* `optim(method="Brent")` -- univariate optimization, calls `optimize`,
  see above.

See the [optim Cheat Sheet](OptimCheatsheet.pdf) for more information about `optim` and the other optimization solvers in Base R, `nlm` and `nlminb`.

*Note*: Default is the derivative-free Nelder-Mead algorithm. If the objective function has derivatives, applying a 'BFGS' variant is definitely to be preferred.  

----

#### **nlminb**

```r
nlminb(start, objective, gradient = NULL, hessian = NULL, ...,
       scale = 1, control = list(), lower = -Inf, upper = Inf
)
```
[help page](./man/nlminb.html)

Uses the [PORT](https://netlib.org/port/) library, a Fortran implementation of quasi-Newton BFGS, for unconstrained and box-constrained optimization.

*Note*: In general, `optim` is to be preferred, or use `ucminf`, see above.

----

#### **constrOptim**

```r
constrOptim(theta, f, grad, ui, ci, mu = 1e-04, control = list(),
            method = if(is.null(grad)) "Nelder-Mead" else "BFGS",
            outer.iterations = 100, outer.eps = 1e-05, ...,
            hessian = FALSE
)
```
[help page](./man/constrOptim.html)

Minimize a function subject to linear inequality constraints using an adaptive barrier algorithm. The feasible region is defined by `ui %*% theta - ci >= 0`. The starting value must be in the interior of the feasible region, but the minimum may be on the boundary.

----

#### **nls**

```r
nls(formula, data, start, control, algorithm,
    trace, subset, weights, na.action, model,
    lower, upper, ...
)
```
[help page](./man/nls.md)

`nls` determines the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model. Lower and upper bounds can only be used with the `port` algorithm. The default is Gauss-Newton, another possibility is `plinear` for partially linear least-squares models).

*Note*: `nls` has problems with 'singular gradients'; instead, use one of the two alternatives, `nlsr` or `minpack.lm`, listed above.

----

----

## Appendix II: External Solvers

### Ipopt through JuliaCall

### The ROI Infrastructure

### Utilizing the NEOS Server

----

----

## References

Boyd, S., and L. Vandenberghe. Convex Optimization. Cambridge University Press, 2004.  
URL: https://web.stanford.edu/~boyd/cvxbook/

Nocedal, J., and S.J. Wright. Numerical Optimization. Springer Verlag, New York 1999.

Price, K.V., R.M. Storn and J.A. Lampinen. Differential Evolution -- A Practical Approach to Global Optimization. Springer Verlag, Berlin Heidelberg 2005.

Schwendinger, F., and H.W. Borchers (2023). CRAN Task View: Optimization and Mathematical Programming. Version 2023-08-19.
URL: https://CRAN.R-project.org/view=Optimization.

Vanderbei. Linear Programming: Foundations and Examples. 5th Edition, Springer Verlag, New York 2020.

Computational Optimization Open Textbook. Cornell University, New York 2020.  
URL: https://optimization.cbe.cornell.edu/index.php?title=Main_Page
