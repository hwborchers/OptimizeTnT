
## Minimax functions and constraints

### Introduction

Minimax problems are problems of the form $\min_x \max_i f_i(x)$ where the $f_i$ are independent functions on the same domain of definition. So the task is to find the minimum of a function $F(x)$ that itself is defined as the maximum of some other functions. These kinds of problems occur often in operations research, decision theory, portfolio management, etc.

The difficulty with minimax functions is that they in general are not smooth at points $x$ where some of the $f_i(x)$ are equal and in most cases they are not smooth at the minimum itself. Optimization solvers that invest in gradients are in danger of being driven away from these non-smooth points, and even an algorithm such as *Nelder-Mead* will have difficulties.

Unfortunately, R does not have specialized solvers for non-smooth problems. But we will show here how minimax problems can be transformed into smooth problems with constraints, problems that some of the available solvers can solve quite effectively.

### Example: Hald-Madsen

As an example, we will look at the Hald-Madsen [??] function defined as the maximum of the 21 functions $f_i(x)$ with

$$f_i(x_1, \ldots, x_5) =|\frac{x_1 + x_2 t_i}{1 + x_3 t_i + x_4 t_i^2 + x_5 t_i^3} - exp(t_i)|, \qquad i = 0, \ldots, 20$$

and $t_i = i/10 - 1$. In R this function can be defined as

```{r}
fHaldMadsen <- function(x) {
    t <- seq(-1.0, 1.0, length.out = 21)
    f <- (x[1]+x[2]*t) / (1+x[3]*t+x[4]*t^2+x[5]*t^3) - exp(t)
    max(abs(f))
}
```

Typically, such functions are not only non-differential at single points, but along the edges where two of the functions $f_i(x)$ cross each other.


## Solving a minimax problem

Instead, we turn it into a smooth problem with added variables and constraints. The function will be as simple as $f(x_1, \ldots, x_5, x_6) = x_6$ with one more variable $x_6$ added. We request that $x_6$ will be minimized under the condition that it still is greater than all the constraint functions $f_i(x_1, \ldots, x_5)$ defined before. These functions have to be reformulated as constraints $x_6 - |f_i(x_1, \ldots, x_5)| \ge 0$. Because of the absolute value condition, each of the inequalities is actually two inequalities.

```{r}
f6 <- function(x) x[6]

hin <- function(x) {
    t <- seq(-1.0, 1.0, length.out = 21)
    h <- (x[1] + x[2]*t)/(1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
    c(x[6] - h, x[6] + h)
}
```

Now one of the available solvers handling constraints can be applied. Here we will use the COBYLA algorithm, implemented in package *nloptr*.

```{r}
x0 <- c(rep(0.5, 5), 0)
opt <- nloptr::cobyla(x0, f6, hin = hin)

opt$par
```

```
[1]  0.9998776288  0.2535884402 -0.7466075719  0.2452015021 -0.0374902911
[6]  0.0001223712
```

Component [6] is our $x_6$, that is the minimum we are looking for. The first five components represent the solution.

This is the true and reported minimum. It may have become clear that non-smooth problems are difficult to solve, but minimax problems can be treated as constrained problems and *can be solved quite fast and accurately* (if the constraints themselves are smooth functions).


## Classical approaches

Here we apply some 'normal' approaches to non-smooth objective functions. We will see that gradient-based and gradient-free solvers are not a good choice for minimax problems. Stochastic solvers, such as Differential Evolution algorithms, may sometimes come quite close to the true minimum.

### Gradient-based approach

Let us try to solve this with a gradient-based procedure in $[-1, 1]^5$ with a starting point $1.0, 0.5, -0.5, 0.5, 0.0$ that is not so far away from the true minimum.

```{r}
x0 <- c(1.0, 0.5, -0.5, 0.5, 0.0)
opt <- optim(rep(0.0, 5), fHaldMadsen,
             method = "L-BFGS-B")

opt$par; opt$value
```

```
[1]  1.03411891  0.63377623 -0.28766165 -0.08473760 -0.02268006
[1] 0.04188082
```

We will see that this is not the true minimum. The minimum value will not change if we use smaller tolerances. With method "BFGS" (without bounds) the minimum will only be slightly lower.

### Gradient-free approach

We try gradient-free approaches such as Nelder-Mead or Hooke-Jeeves. Nelder-Mead without bounds is available through `optim()`:


```{r}
opt <- optim(x0, fHaldMadsen,
             method = "Nelder-Mead")

opt$par; opt$value
```

```
[1]  1.0054041  0.6191019 -0.3082727 -0.1936943  0.1066206
[1] 0.03161394
```

Instead, we can make use of solvers in *dfoptim*, a package for gradient-free optimization.

```{r}
opt <- dfoptim::nmk(x0, fHaldMadsen)

opt$par; opt$value
```

```
[1]  1.0177232 -0.4001709 -1.6546850  1.5588244 -0.6838769
[1] 0.08544269
```

or of the adaptive Nelder-Mead implementation in *pracma*.

```{r}
opt <- pracma::anms(fHaldMadsen, x0)

opt$xmin; opt$fmin
```

```
[1]  1.003411798  0.344928111 -0.645439154  0.151495264 -0.009215026
[1] 0.004458874
```

This result is the best of all our attempts by now, still far away from the true minimum.

### Solve with stochastic optimizer

To make it more plausible, a stochastic global solver such as could be used. Beware that a high number of iterations will be needed and the solution is not overly accurate.

```{r}
    library(DEoptim)
```

```
DEoptim package
Differential Evolution algorithm in R
Authors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich
```

```{r}
    sol <- DEoptim(fHaldMadsen, lower = rep(-1, 5), upper = rep(1, 5),
                   DEoptim.control(itermax = 5000, trace = FALSE))
    sol$optim
```

```
$bestmem
       par1        par2        par3        par4        par5 
 0.99987763  0.25358844 -0.74660757  0.24520150 -0.03749029 

$bestval
[1] 0.0001223713
```

The reader may try other optimization routines with randomly selected starting points to convince himself that we have found the minimum.

