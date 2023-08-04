## The 'unit-length' Constraint

The unit-length restriction or $sum(x^2) = 1$ constraint will appear in many numerical, technical, and data analysis applications. It is a quadratic constraint and cannot be easily eliminated, compared to eliminating a variable given the linear constraint $sum(x) = 1$.

We will look at different approaches to solve optimization problems under the unit-length constraint.

```r
  library(adagio)
  fn <- fnRosenbrock
```


### Stereographic Projection

The **stereographic projection** is implemented in package 'pracma' as `stereographic()` with its inverse as `stereographic_inv()`. Imagine the sphere $S^n \in R^{n+1}$ defined by

$$
  S^n = \\{x \in R^{n+1} |\\; ||x|| = 1\\}
$$

The stereographic projection takes the south pole $Sp = (0,...,0, -1)$ as the center of projection. Each point on the sphere will be projected onto a point on the tangent plane at the north pole $Np = (0,...,0, +1)$.

The inverse projection is defined similarly. The procedure defines a bijection (actually a *diffeomorphism*) between $S$ minus the south pole and a plane representing $R^{n-1}$.


### Stereographic projection and optimization

As an example objective function we will take the Rosenbrock function in 'adagio' in 30 dimensions (to make it a bit more challenging). For utilizing the stereographic projection define a function `fn()` of $n-1$ variables on the tangent plane.

```r
  fn = function(x) {
      x1 <- stereographic_inv(c(x, 1))
      adagio::fnRosenbrock(x1)
}
```

We will start near the north pole (in 30 dimensions)
```r
n <- 30
set.seed(1799)
x0 <- runif(n-1)    # one freedom lost on the tangent plane
```

and apply `anms()` from the 'pracma' package, that is

```r
library(pracma)
system.time( sol <- pracma::anms(fn = fn, x0) )
##    user  system elapsed 
##   2.197   0.000   2.203 
```

The result is displayed here. To find the minimal point we have to project the intermediate solution from the tangent plane back to the sphere.

```r
sol$fmin
## [1] 26.2277

my_xmin <- c(stereographic_inv(c(sol$xmin, 1)))
##  [1]  7.474359e-01 5.641682e-01 3.265019e-01 1.160691e-01 2.336713e-02
##       ...
## [26]  1.005602e-02 1.005595e-02 1.005211e-02 9.860475e-03 9.585063e-05

sum(my_xmin^2)
## [1] 1
```

and we can verify the minimum by calling the Rosenbrock function on it.

```r
adagio::fnRosenbrock(my_xmin)
## [1] 26.2277
```


### General constraint optimization

There exist several CRAN packages that will solve constrained optimization
problems: 'alabama' (with `auglag()`), 'NlcOptim' (`solnl()`), or 'Rsolnp' 
(`solnp()`).  We will utilize one of the most reliable and accurate constraint solvers, `slspq()` in the 'nloptr' package.

```r
library(nloptr)
p0 <- c(x0, 1)
heq <- function(x) 1 - sum(x*x)

system.time(ans <- slsqp(p0, fnRosenbrock, heq=heq))
##    user  system elapsed 
##   0.097   0.000   0.097 

ans$value
## [1] 26.2277
```

We can compare these two solutions by plotting the function along the line from one point to the other; the function `flineviz` in 'adagio will do this for us. Here we compare the minimal values.

```r
print(c(sol$fmin, ans$value), digits=12)
## [1] 26.2276990802 26.2276990802
```
