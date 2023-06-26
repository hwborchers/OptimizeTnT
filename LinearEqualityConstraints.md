## Linear equality constraints only

If there are only *linear equality* constraints there is another 'trick'
from linear algebra that will enable this problem to be solved without
constraints at all, even if the objective function itself is nonlinear.

Assume $f$ is a nonlinear function and the task is to minimize this
function under linear equality constraints only. The linear equalities
shall be given in matrix notation as $A x = d$ with matrix $A$ and
(column) vector $d$. Then let $N$ be the *kernel* or *null space* of the
matrix $A$, i.e. the subspace of all vectors $x$ such that $A x = 0$,
let $b_1, \ldots, b_m$ a basis of $N$ and $x_0$ a least-squares solution
of the linear equation problem $A x = d$.

The original constrained minimization task can now be formulated as an
optimization problem *without* constraints! Define a (nonlinear)
function $g$ on $\mathrm{R}^m$ as
$$g(x_1, \ldots, x_m) = f(x_0 + x_1 b_1 + \ldots + x_m b_m)$$ The point
$\bar{x} = x_0 + x_1 b_1 + \ldots + x_m b_m$ is a least-squares solution
of $A x = d$ and the optimization problem will return a minimum for the
function $g$. If $A x_0 \ne d$ there is anyway no solution to the
original constrained problem with equality constraints.

How can this be done in ? We attempt to minimize Rosenbrock's function
(in 3 dimensions) along the line $x_1 + x_2 + x_3 = 1$. For determining
the null space we use function in package .

    > require(pracma)
    > A <- matrix(c(1, 1, 1), nrow = 1)
    > N <- null(A); N
               [,1]       [,2]
    [1,] -0.5773503 -0.5773503
    [2,]  0.7886751 -0.2113249
    [3,] -0.2113249  0.7886751

    > x0 <- qr.solve(A, 1)      # (1, 0, 0) exact solution of linear equation

The columns in $N$ are a basis of the nullspace. We define a new
function in two dimensions as

    > g <- function(x)
    +     f(x0 + x[1]*N[, 1] + x[2]* N[, 2])

Optimizing this functions yields

    > sol <- optim(c(0, 0), g, method = "Nelder-Mead")
    > sol$par
    [1] 0.4832432 0.2591742

and the minimum of the original problem can be uncovered as

    > xmin <- x0 + sol$par[1]*N[, 1] + sol$par[2]* N[, 2]
    > xmin; f(xmin)
    [1]  0.5713651 0.3263519 0.1022830
    [1]  0.6393138

This approach can easily be implemented as a solver for nonlinear
optimization problems with linear equality constraints. Function in does
this for least-squares problems with equality constraints (see below
section ??). Unfortunately, this cannot be easily extended to the case
when there are additional (linear) inequality constraints, because the
inequalities would have to be expressed in the basis of the null space.
