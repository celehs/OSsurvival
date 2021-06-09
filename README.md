
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Optimal Surrogate Survival: Quantifying the Feasibility of Shortening Clinical Trial Duration Using Surrogate Markers

<!-- badges: start -->

<!-- badges: end -->

The goal of Optimal Surrogate Survival (OSsurvival) is to
nonparametrically estimate the proportion of treatment on the primary
outcome explained by a surrogate (PTE) of an optimally transformation of
surrogate marker measured at an earlier time. The primary outcome
measured at a later time may be subject to censoring.

## Installation

The development version of this package can be installed from GitHub
with:

    library(devtools)
    devtools::install_github("celehs/OSsurvival")

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(OSsurvival)
# load the data
data("sysdata")

# time at which the surrogate is measured
t.0 = data.example$t.0

# time at which the primary outcome is measured
t = data.example$t

# observed survival time
xob = data.example$data$xob

# surrogate information at t.0
s.ob = data.example$data$s.ob

# event indicator
deltaob = data.example$data$deltaob

# treatment indicator
aob = data.example$data$aob

# main estimation function
# varind: whether to estimate variance; re:number of replications for resampling
out = pte.survival(xob, s.ob, deltaob, aob, t, t.0, varind=0, re=100)

# estimated PTE
out$pte.est
#> [1] 0.8413843

# estimated g1
out$g1.est
#> [1] -0.0180641

# estimated g2(s) at equally spaced s point
plot(out$sgrid, out$gs.est, type="l", xlab = "Surrogate Marker", ylab = "Optimal Transformation")
```

<img src="man/figures/README-example-1.png" width="100%" /> The PTE
result indicates that this is a moderate to high surrogate marker in
this setting.

## Some background

The observable data for analysis consist of \(n\) sets of independent
and identically distributed random vectors
\(\{D_i = (X_i,\ \delta_i,\ S_i I(X_i \ge t_0),\ A_i),\ i = 1, ..., n\}\),
where \(T_i = T_i^{(1)}A_i + T_i^{(0)}(1-A_i)\),
\(C_i = C_i^{(1)}A_i + C_i^{(0)}(1-A_i)\), \(X_i=min(T_i, C_i)\), the
primary outcome \(Y_i=I(X_i>t)\), the surrogate
\(S_i = S_i^{(1)}A_i + S_i^{(0)}(1-A_i)\) is only observed for those
with \(X_i > t_0\).

The function outputs estimates (and standard error estimates if
indicated) of \(PTE =\Delta_{g_{opt}(S)} / \Delta\), the proportion of
treatment effect explained quantity based on the ratio between the
treatment effect on the optimal transformation of the potential
surrogate marker and the treatment effect on the primary outcome, where
\(\Delta\) is the treatment effect on the primary outcome \(Y\);
\(\Delta_{g_{opt}(S)}\) is the treatment effect on the optimal
transformation of the surrogate
\(g_{opt}(S)=I(T\leq t_0)\ g_{1, opt}+I(T>t_0)\ g_{2, opt}\approx I(T>t_0)\ g_{2, opt}\).

## Citation

Wang X, Cai T, Tian L, Bourgeois F, Parast L. Quantifying the
Feasibility of Shortening Clinical Trial Duration Using Surrogate
Markers. Under revision.
