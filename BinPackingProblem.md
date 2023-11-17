# Bin Packing Problem

Author: Hans W Borchers  
Date: November 15, 2023

## Introduction

### The problem

The **Bin Packing problem** (BPP) is a *combinatorial* optimization problem that addresses the following question: 

- Given a finite number $n$ of items with weights $w_1, \ldots, w_n$  
  and a set of $k$ 'bins', all of capacity $b$, 

- put the items into the bins such that no bin is overloaded.

- What is the minimal number of bins (of this capacity) that is needed for the task?

Obviously this is an important problem in many practical applications, for example in transportation and packaging industries.

The problem is **strongly NP-hard** and cannot be solved in polynomial time, as Garey and Johnson (1979) have shown in their path-breaking book "Computers and Intractability". Still there exist sophisticated algorithms that can solve quite large instances of the problem in reasonabble time.

### Approximate solutions

There also exist approximate and heuristic algorithms that may come quite close to the true solution. Function `bpp_approx` in package {adagio} provides three such heuristic approaches as 'first fit', best fit, and 'worst fit'. These are *offline* algorithms, ie. the set of items is known in advance and can for instance be sorted decreasingly before applying the algorithm.

Assume we are given the following sizes of items plus a capacity of our bins.


::: {.cell}

```{.r .cell-code}
# Examples1:
S1 <- c(50, 3, 48, 53, 53, 4, 3, 41, 23, 20, 52, 49)        # len=12, k=5
S2 <- c( 8, 19, 40, 51, 21, 48, 85, 68, 86, 79,
        93, 93, 21, 73, 72, 77,  7, 41, 14, 70)             # len=20, k=12
S3 <- c(100, 98, 96, 93, 91, 87, 81, 59, 58, 55, 50, 43,
        22, 21, 20, 15, 14, 10,  8,  6,  5,  4,  3,  1, 0)  # len=25, k=11

# Take S2, sort if is.unsorted(S)
S <- S2
S <-sort(S, decreasing = TRUE)

n <- length(S)
b <- 100
```
:::


The calling `bpp_approx` with different methods gives (contrary to its name, the "worstfit" method appears to return a good approximation)


::: {.cell}

```{.r .cell-code}
library(adagio)
sol <- bpp_approx(S, b, method = "bestfit")
m <- sol$nbins
cat("Method: best fit", '\tnbins =', m, '\t filled =', round(100*sol$filled, 1),
    "%", '\n')
```

::: {.cell-output .cell-output-stdout}
```
Method: best fit 	nbins = 12 	 filled = 88.8 % 
```
:::
:::



### Theoretical model

Following Martello and Toth in their book "Knapsack problems: Algorithns and Computer Implementations" resp. [Chapter 8](http://www.or.deis.unibo.it/kp/Chapter8.pdf) of this book -- the Bin Packing Problem can be modelled as a linear program.

For details see also the entry in the [Wikipedia: Bin packing problem](https://en.wikipedia.org/wiki/Bin_packing_problem).

Given a finite set $S$ of natural numbers and a capacity $b$, decision variables $x_{ij}$ determine whether item $j$ is put into the $j$-th bin, and variables $y_k$ are $1$ if the $k$-th bin is used (not empty). The constraints are expressed as
$\sum_j x_{ij} = 1$ and $\sum_i S_i x_{ij} \le b y_j$. The objective to minimize is $\sum_j y_j$ as the number of non-empty bins.

### Exact solutions

We will need two types of decision variables, $x_{ij} = 1$ if the i-th item is put into the j-th bin, and $0$ otherwise, and decision variables $y_l = 1$ if the l-th bin is used, and equal $0$ otherwise (ie., there are no items in this bin). If $n$ is the number of items and $m$ the number of bins, then the total number of decision variables is $m \times n + m = (n+1)\, m$.

The number of items is know in advance; the maximum number of bins is obviously $m = n$. This will result in a high number of binary variables. A better estimate is to use the number of bins calculated by an approximate binpacking algorithm such as `approx_bpp()` in 'adagio'.

The Linear Programming tools require the decision variables to be arranged in sequence, that is as a vektor. So we will utilize a vector `x` of length `n*m+m` representing the decision. For example, `x[1:m]` determines in which bin the first item is to be placed.

We use the values for `S`, `n=25` from above and `m=15` to show that it will be reduced. There are `n` constraints declaring that each item has to be put in exactly one bin, and `m` constraints limiting the sum of weights of items put into them. The matrix `A` shall represent these constraints, so let's start filling it.


::: {.cell}

```{.r .cell-code}
n <- length(S)                  # S the (unsorted) weights of items
m <- 15  # sol$nbins            # maximum number of needed bins
A = matrix(0, n+m, (n+1)*m)     # n+m constraint for n*m+m variables

# sum(j=1..m, x_ij) == 1 for all i
for (i in 1:n) A[i, ((i-1)*m + 1):(i*m)] = 1
```
:::


Next the constraint that the sum of weights of item in each bin has to be smaller than the capacity.


::: {.cell}

```{.r .cell-code}
# sum(i=1..n, s_i*x_ij) <= b*y_j for all j
inds = seq(1, m*(n+1), by = m)
for (i in 1:m) A[n+i, inds+(i-1)] = c(S, -b)
```
:::


And finally there are the coefficients in `L` of the objective function, the `types` and lower and upper bounds `0` and `1`, requesting that the decision 
variables are *binary* (`0`and `1` their only only allowed values), and constraints on the sum of weights for each bin.


::: {.cell}

```{.r .cell-code}
L = c(rep(0, m*n), rep(1, m))   # coefficients of the objective
types = rep('I', m*(n+1))       # all variables are integer
lb = (rep(0, m*(n+1)))          #     of values 0 or 1, ie. binary
ub = c(rep(1, m*(n+1)))

lhs = c(rep(1, n), rep(-b, m))  # 1 <= x_i1+...+x_im <= 1 for all i
rhs = c(rep(1, n), rep(0, m))   # 
```
:::


Now we solve the Bin Packing problem with the Highs solver.


::: {.cell}

```{.r .cell-code}
## system.time: 2.01 sec
library(highs)
system.time(
sol <- highs_solve(L = L, lower = lb, upper = ub,
                   A = A, lhs = lhs, rhs = rhs, types = types)
)
```

::: {.cell-output .cell-output-stdout}
```
   user  system elapsed 
  1.295   0.028   1.322 
```
:::

```{.r .cell-code}
sol$objective_value; sol$status_message
```

::: {.cell-output .cell-output-stdout}
```
[1] 12
```
:::

::: {.cell-output .cell-output-stdout}
```
[1] "Optimal"
```
:::
:::


To see how much each bin is filled, we use the 'primal solution' and reshape this into the matrix as it was originally intended. The (n+1)-th row represents the decision variable $y_i$, their sum is the number of bins filled.


::: {.cell}

```{.r .cell-code}
M <- matrix(sol$primal_solution, nrow=n+1, byrow = TRUE)

sum(M[n+1, ])
```

::: {.cell-output .cell-output-stdout}
```
[1] 12
```
:::
:::


Three bins are empty; the solution of a bin packing problem is obviously not unique. How can we generate the full solution from this, meaning which item weight goes into which bin? First find the bins that are filled.


::: {.cell}

```{.r .cell-code}
inds <- which(M[n+1, ] == 1)

for (j in inds) {
    bj <- M[1:n, j]*S; bj <- bj[bj != 0]
    cat("bin", j, '\t', sum(bj), '\t', bj, '\n')
}
```

::: {.cell-output .cell-output-stdout}
```
bin 1 	 99 	 77 14 8 
bin 2 	 87 	 68 19 
bin 3 	 94 	 73 21 
bin 6 	 93 	 93 
bin 7 	 85 	 85 
bin 8 	 72 	 72 
bin 9 	 70 	 70 
bin 10 	 93 	 86 7 
bin 11 	 93 	 93 
bin 12 	 99 	 51 48 
bin 13 	 81 	 41 40 
bin 14 	 100 	 79 21 
```
:::
:::


The other bins are empty. 


## The `binpacking` function

By summarizing all the steps we have done so far, we define a function 'binPacking()' that executes all intermediate calculations from the input and displays a result and a short report on the command line.

Input variables are

* the set `S` of integers (or positive real numbers)
* the capacity of the bins (the same for all of them)
* a logical variable indicating to generate a report (or not)

The output is a list with the following components:

* nbins -- the number of bins used
* load -- the load for each bin
* bins -- for each bin the weights put into this bin
* status -- a status message ('Optimal' or 'Infeasible')

Because an approximate solution has already been found, the status 'Infeasible' should never occur -- if it does, there is a programming or input error somewhere.



::: {.cell}

```{.r .cell-code}
# library(adagio)
# library(highs)

binpacking <- function(S, cap, printing = TRUE) {
    stopifnot(is.numeric(S), is.numeric(cap), length(cap) == 1)
    if (any(S > cap) || any(S <= 0))
        stop("Set 'S' must be integers 0 < S[i] <= cap.")
    if (any(diff(S) > 0)) S <- sort(S, decreasing = TRUE)
    b <- cap         # bin capacity
    n <- length(S)   # no. of items

    # Best approximative number of bins
    k1 <- adagio::bpp_approx(S, b, "firstfit")$nbins
    k2 <- adagio::bpp_approx(S, b, "bestfit")$nbins
    k3 <- adagio::bpp_approx(S, b, "worstfit")$nbins
    k <- min(k1, k2, k3)

    # Coefficients of the objective function
    L <- c(rep(0, n*k), rep(1, k))      # coefficients of y

    # Set decision variables as binary
    types <- rep('I', n*k + k)          # all are integer variables
    lb = c(rep(0, n*k + k))             # make them binary
    ub = c(rep(1, n*k + k))

    # Define the constraint matrix
    A = matrix(0, n+k, n*k + k)         # constraint matrix
    for (i in 1:n) A[i, ((i-1)*k + 1):(i*k)] = 1
    inds = seq(1, (n-1)*k + 1, by = k)
    # for (i in 1:n) A[n+i, c(inds+(i-1), n*k + i)] = c(S, -b)  # WRONG ?
    for (i in 1:k) A[n+i, inds + (i-1)] = S
    for (i in 1:k) A[n+i, n*k + i] = -b
    lhs = c(rep(1, n), rep(-Inf, k))    # sum(x) <= b*y 
    rhs = c(rep(1, n), rep( 0, k))

    # Apply the MILP solver from package 'highs'
    st <- system.time(
      sol <- highs_solve( L = L, lower = lb, upper = ub,
                          A = A, lhs = lhs, rhs = rhs, types = types)
    )
    elapsed <- st[3]
    if (sol[["status_message"]] != "Optimal")   # "Infeasible"
        stop("No optimal solution found. (Reason unknown.)")

    # Prepare results for output
    primal <- round(sol[["primal_solution"]])
    M = matrix(primal, nrow = n+1, byrow = TRUE)
    
    l <- 0; Lst <- list()
    for (j in 1:k) {
        if (M[n+1, j] != 0) {       # this bin is being used
            z = M[1:n, j] * S       # fill with values from S
            z = z[z != 0]           # reduce to items in this bin
            l <- l+1
            Lst[[l]] <- z
        }
    }
    loading <- sapply(Lst, sum)       # sum weights in each bin
    perc.filling <- round(100*sum(loading)/length(loading)/cap, 1)

    # Print details to the command line
    if (printing) {
        cat("= BIN PACKING PROBLEM =====", '\n')
        cat("  S = {", max(S), ", ..., ", min(S), "}",
            "of lenght", n, '\n')
        cat("  sum(S) =", sum(S), '\n')
        cat("  capacity =", cap, '\n')
        cat("  Approxi solution:", k, '\n')
        cat("  Optimal solution:", sol$objective_value, '\n')
        cat("  Bin loading: ", loading, '\n')
        cat("  Perc. filling:", perc.filling, "%", '\n')
        cat("  Status message:", sol$status_message, '\n')
        cat("  Computing time:", elapsed, "[s]", '\n')
        cat("= =========================", '\n\n')
    }

    # Return solution as a list
    return(list(nbins = sol$objective_value,
                load = loading,
                bins = Lst,
                status = sol$status_message)
    )
}
```
:::



## Examples

### Simple examples

These are some examples found in survey articles; they do not pose problems to our function `binpacking`.


::: {.cell}

```{.r .cell-code}
cap <- 100
S1 <- c(50, 3, 48, 53, 53, 4, 3, 41, 23, 20, 52, 49)    # length 12

sol <- binpacking(S1, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 53 , ...,  3 } of lenght 12 
  sum(S) = 399 
  capacity = 100 
  Approxi solution: 4 
  Optimal solution: 4 
  Bin loading:  100 99 100 100 
  Perc. filling: 99.8 % 
  Status message: Optimal 
  Computing time: 0.002 [s] 
= ========================= 
```
:::
:::

::: {.cell}

```{.r .cell-code}
cap <- 100
S2 <- c( 8, 19, 40, 51, 21, 48, 85, 68, 86, 79,     # length 20
        93, 93, 21, 73, 72, 77,  7, 41, 14, 70)

sol <- binpacking(S2, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 93 , ...,  7 } of lenght 20 
  sum(S) = 1066 
  capacity = 100 
  Approxi solution: 12 
  Optimal solution: 12 
  Bin loading:  95 85 100 86 93 91 68 91 73 92 93 99 
  Perc. filling: 88.8 % 
  Status message: Optimal 
  Computing time: 0.626 [s] 
= ========================= 
```
:::
:::

::: {.cell}

```{.r .cell-code}
cap <- 100
S3 <- c(100, 98, 96, 93, 91, 87, 81, 59, 58, 55, 50, 43,    # length 24
         22, 21, 20, 15, 14, 10,  8,  6,  5,  4,  3,  1)
sol <- binpacking(S3, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 100 , ...,  1 } of lenght 24 
  sum(S) = 1040 
  capacity = 100 
  Approxi solution: 11 
  Optimal solution: 11 
  Bin loading:  96 70 98 93 96 100 100 91 98 98 100 
  Perc. filling: 94.5 % 
  Status message: Optimal 
  Computing time: 0.005 [s] 
= ========================= 
```
:::
:::


Here is a small example, where the optimal solution and the approximate solution differ. Apparently this happens when the optimal solution almost fills the bins to their maximum.


::: {.cell}

```{.r .cell-code}
cap = 100
S4 <- c(49, 41, 34, 33, 29, 26, 26, 22, 20, 19)     # length 10
sol <- binpacking(S4, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 49 , ...,  19 } of lenght 10 
  sum(S) = 299 
  capacity = 100 
  Approxi solution: 4 
  Optimal solution: 3 
  Bin loading:  99 100 100 
  Perc. filling: 99.7 % 
  Status message: Optimal 
  Computing time: 0.052 [s] 
= ========================= 
```
:::
:::


### More  and self-generated examples

Here is a more extended example


::: {.cell}

```{.r .cell-code}
S <- c( 99,99,96,96,92,92,91,88,87,86,
        85,76,74,72,69,67,67,62,61,56,
        52,51,49,46,44,42,40,40,33,33,
        30,30,29,28,28,27,25,24,23,22,
        21,20,17,14,13,11,10, 7, 7, 3 )     # length 50
cap <- 100
sol <- binpacking(S, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 99 , ...,  3 } of lenght 50 
  sum(S) = 2434 
  capacity = 100 
  Approxi solution: 25 
  Optimal solution: 25 
  Bin loading:  100 99 97 94 100 99 88 100 99 99 99 96 99 99 98 99 100 86 100 92 98 96 99 100 98 
  Perc. filling: 97.4 % 
  Status message: Optimal 
  Computing time: 0.526 [s] 
= ========================= 
```
:::
:::

 
We can also generate our own examples by simply generating random numbers and select a capacity that appears difficult to fill.


::: {.cell}

```{.r .cell-code}
S <- c(100:1)               # Gauss' school excercise
cap <- 505
sol <- binpacking(S, cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 100 , ...,  1 } of lenght 100 
  sum(S) = 5050 
  capacity = 505 
  Approxi solution: 11 
  Optimal solution: 10 
  Bin loading:  505 505 505 505 505 505 505 505 505 505 
  Perc. filling: 100 % 
  Status message: Optimal 
  Computing time: 3.477 [s] 
= ========================= 
```
:::
:::


### Splitting into subsets

Task: Split the set of prime numbers
`S = {19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97}`
in two disjoint subsets S1 and S2, such that sum(S1) = sum(S2) !


::: {.cell}

```{.r .cell-code}
S <- c(19, 23, 29, 31, 37, 41, 43, 47, 53,
       59, 61, 67, 71, 73, 79, 83, 89, 97)
cap = sum(S)/2      # 501

sol <- binpacking(rev(S), cap)
```

::: {.cell-output .cell-output-stdout}
```
= BIN PACKING PROBLEM ===== 
  S = { 97 , ...,  19 } of lenght 18 
  sum(S) = 1002 
  capacity = 501 
  Approxi solution: 3 
  Optimal solution: 2 
  Bin loading:  501 501 
  Perc. filling: 100 % 
  Status message: Optimal 
  Computing time: 0.028 [s] 
= ========================= 
```
:::
:::


To see how these two subsets are filled with numbers from 1 to 100, look at `sol$bins`.


::: {.cell}

```{.r .cell-code}
sol$bins
```

::: {.cell-output .cell-output-stdout}
```
[[1]]
[1] 97 79 73 67 59 43 31 29 23

[[2]]
[1] 89 83 71 61 53 47 41 37 19
```
:::
:::


You might try to add the primes from 2 (or 3) to 17; there is also an exact spitting for these sets!

Such a result can also be found by applying a 'subset sum' algorithm. For more than two subsets it is more tricky.


## Using other LP solvers

We will once again utilize the following example that is given in detail again to fix the notation. The reader can skip this and go directly to the comparison of the different LP solvers.


::: {.cell}

```{.r .cell-code}
S <- c(100, 98, 96, 93, 91, 87, 81, 59, 58, 55, 50, 43,
        22, 21, 20, 15, 14, 10,  8,  6,  5,  4,  3,  1, 0)

b <- 100         # bin capacity
n <- length(S)   # no. of items

# max no. of bins needed
k <- adagio::bpp_approx(sort(S, decreasing = TRUE), b)$nbins
```
:::


There are `n*n + k` decision variables, the $n^2$ `x`s and k `y`s. The linear part of the objective includes only `y`.

Our variables are *integer*, but we want them to be *binary*, so we define lower and upper boundaries 0 and 1 for them.


::: {.cell}

```{.r .cell-code}
L <- c(rep(0, n*k), rep(1, k))  # coefficients of y
types <- rep('I', n*k + k)      # all are integer variables

lb = c(rep(0, n*k + k))
ub = c(rep(1, n*k + k))
```
:::


We define the matrix `A` that puts together all the restrictions on the binary decision variables.


::: {.cell}

```{.r .cell-code}
A = matrix(0, n+k, n*k + k)

for (i in 1:n) A[i, ((i-1)*k + 1):(i*k)] = 1

inds = seq(1, (n-1)*k + 1, by = k)
# for (i in 1:n) A[n+i, c(inds+(i-1), n*k + i)] = c(S, -b)  # WRONG !

for (i in 1:k) A[n+i, inds + (i-1)] = S
for (i in 1:k) A[n+i, n*k + i] = -b
```
:::


The conditions for matrix `A` are


::: {.cell}

```{.r .cell-code}
lhs = c(rep(1, n), rep(-Inf, k))  # sum(x) <= b*y 
rhs = c(rep(1, n), rep( 0, k))   # 
```
:::


Now we apply our solvers to this model.

### highs

HiGHS is high performance serial and parallel software for solving large-scale sparse linear programming (LP), mixed-integer programming (MIP) and quadratic programming (QP) models, developed in C++11, with interfaces to C, C#, FORTRAN, Julia and Python.

HiGHS is freely available under the MIT licence


::: {.cell}

```{.r .cell-code}
library(highs)

## system.rime: 2.01 sec
system.time(
  sol <- highs_solve(L = L, lower = lb, upper = ub,
                     A = A, lhs = lhs, rhs = rhs, types = types)
)
```

::: {.cell-output .cell-output-stdout}
```
   user  system elapsed 
  0.005   0.000   0.005 
```
:::

```{.r .cell-code}
##    user  system elapsed 
##   0.951   0.000   0.950 
```
:::


and retrieve the solution.


::: {.cell}

```{.r .cell-code}
sol$objective_value   # the no. of bins used
```

::: {.cell-output .cell-output-stdout}
```
[1] 11
```
:::

```{.r .cell-code}
# the matrix of item to bin assignments
m <- matrix(sol$primal_solution[1:(n*k)], nrow=n, byrow = TRUE)

# the filling of the bins
filling <- apply(m * S, 2, sum)
sort(filling, decreasing = TRUE)
```

::: {.cell-output .cell-output-stdout}
```
 [1] 100 100 100  98  98  98  96  96  93  91  70
```
:::

```{.r .cell-code}
##  [1] 100  99  98  93  93  93  91  91  88  79  73  68
```
:::


### lpSolve

Lp_solve is freely available software for solving linear, integer and mixed integer programs. The current version is 5.5.

```r
library(lpSolve)
k = 13

system.time(
  sol <- lp("min",
    objective.in = L,
    const.mat = A,
    const.dir = c(rep("==", n), rep("<=", n)),
    const.rhs = rhs,
    all.bin = TRUE
  )
)
```

does not terminate !

### Rglpk

The GNU Linear Programming Kit (GLPK) is a software package intended for solving large-scale linear programming (LP), mixed integer programming (MIP), and other related problems.

```r
library(Rglpk)

system.time(
  sol <- Rglpk_solve_LP(obj = L,
                        mat = A,
                        dir = c(rep("==", n), rep("<=", k)),
                        rhs = rhs,
                        # bounds = NULL,
                        types = "B",
                        max = FALSE  #, control = list()
  )
)
```

### Rsymphony

SYMPHONY is an open-source COIN-OR project and a generic MILP solver, callable library, and extensible framework for implementing customized solvers for mixed-integer linear programs (MILPs).

```r
library(Rsymphony)

system.time(
  sol <- Rsymphony_solve_LP(obj = L,
                            mat = A,
                            dir = c(rep("==", n), rep("<=", k)),
                            rhs = rhs,
                            types = "B")
)
```

     user  system elapsed 
    9.970   0.321  10.288 

solves it correctly, but is much slower.

### Cbc, Clp

CBC is an open-source MILP solver. It uses many of the COIN-OR components and is designed to be used with CLP. It is available as a library and as a standalone solver.

```r
library(rcbc)

system.time(
  sol <- cbc_solve(
    obj = L,
    mat = A,
    row_ub = c(rep(1, n), rep( 0, k)),
    row_lb = c(rep(1, n), rep(-Inf, k)),
    col_lb = rep.int(0, n*k + k),
    col_ub = rep.int(1, n*k + k),
    is_integer = rep.int(TRUE, n*k + k),
    # max = FALSE,
    # cbc_args = list(),
    initial_solution = NULL
  )
)
```
     user  system elapsed 
    0.357   0.206   0.329 

### adagio::binpacking

This is a Fortran implementation of Martello and Toth's well known algorithm for bin packing, published in their book "Knapsack Problems". It is integrated into the 'knapsack' package available from R-Forge. Please note that it is free for academic use and otherwise under ACM license -- which is the reason it's not on CRAN (the CRAN administrators don't like ACM licenses).

```r
library(knapsack)

S2 <- c( 8, 19, 40, 51, 21, 48, 85, 68, 86, 79,
        93, 93, 21, 73, 72, 77,  7, 41, 14, 70)
system.time(
  {S <- sort(S, decreasing=TRUE)
  sol <- binpacking(S, b)}         # $nbins, $xbins
)
```

     user  system elapsed 
    0.001   0.000   0.007 

The distribution of weights into bins and the filling of bins can be retrieved as follows:

```r
sol$fill <- numeric(sol$nbins)

for (i in 1:sol$nbins) {
    is <- which(sol$xbins == i)
    # cat(i, ':', is, '\n')
    sol$fill[i] <- sum(S[is])
}

sort(sol$fill, decreasing = TRUE)
```

    [1] 100 100 100  99  98  93  93  92  81  72  70  68

This is a specialized executable solver in Fortran and famous for its speed. It can be used for medium-sized problems where even solvers such as CBC or Highs will not succeed.
    
## Conclusion

Success resp. running times of the different solvers are summarized in the following table.

    | package   | solver             | timing [secs] | comment              |
    |-----------|--------------------|---------------|----------------------|
    | lpSolve   |                 lp |            NA |    did not terminate |
    | Rglpk     |     Rglpk_solve_LP |            NA |    did not terminate |
    | Rsymphony | Rsymphony_solve_LP |     9.97 secs |                      |
    | rcbc      |          cbc_solve |     0.36 secs |                      |
    | highs     |        highs_solve |     0.95 secs |                      |
    | adagio    |         binpacking |     0.01 secs | BPP solver (Fortran) |

