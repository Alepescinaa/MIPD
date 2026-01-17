
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MIPD

<!-- badges: start -->

<!-- badges: end -->

The goal of MIPD is to …

## Installation

You can install the development version of MIPD from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Alepescinaa/MIPD")
```

## Example

This is a basic example:

``` r
library(MIPD)
```

load dataset

``` r
data(toy_data_MM)
```

To inspect help file main function:

``` r
?run_mipd
```

run imputation and illness death model

``` r

fit <- run_mipd(
  toy_data_MM,
  m = 5,
  cov_vector = c("cov1", "cov2", "cov3"),
  clock_assumption = "forward",
  distribution = "gompertz",
  inner_cores = 5,
  max_iter = 20,
  eps = 1e-3
)
#> [[1]]
#>    cov1.1    cov2.1    cov3.1 
#> 0.4660372 0.6959313 0.1671973 
#> 
#> [[2]]
#>    cov1.2    cov2.2    cov3.2 
#> 0.4054468 0.6827780 0.2468225 
#> 
#> [[3]]
#>    cov1.3    cov2.3    cov3.3 
#> 0.6122589 0.2680092 0.2573598
#> Iteration 1 — Total criteria: 3.28216 (c1: 3.15011, c2: 0.06603,  c3: 0.06603)
#> [[1]]
#>    cov1.1    cov2.1    cov3.1 
#> 0.4918652 0.7016724 0.1696413 
#> 
#> [[2]]
#>    cov1.2    cov2.2    cov3.2 
#> 0.4117957 0.6832464 0.2474701 
#> 
#> [[3]]
#>    cov1.3    cov2.3    cov3.3 
#> 0.5815368 0.2678758 0.2562380
#> Iteration 2 — Total criteria: 0.11623 (c1: 0.07829, c2: 0.01897,  c3: 0.01897)
#> [[1]]
#>    cov1.1    cov2.1    cov3.1 
#> 0.4950971 0.7023917 0.1718428 
#> 
#> [[2]]
#>    cov1.2    cov2.2    cov3.2 
#> 0.4115215 0.6829650 0.2467260 
#> 
#> [[3]]
#>    cov1.3    cov2.3    cov3.3 
#> 0.5823834 0.2692077 0.2587982
#> Iteration 3 — Total criteria: 0.02874 (c1: 0.02057, c2: 0.00408,  c3: 0.00408)
#> [[1]]
#>    cov1.1    cov2.1    cov3.1 
#> 0.4947344 0.7026579 0.1710610 
#> 
#> [[2]]
#>    cov1.2    cov2.2    cov3.2 
#> 0.4120585 0.6829276 0.2469169 
#> 
#> [[3]]
#>    cov1.3    cov2.3    cov3.3 
#> 0.5809241 0.2693909 0.2581687
#> Iteration 4 — Total criteria: 0.00993 (c1: 0.00566, c2: 0.00213,  c3: 0.00213)
#> [[1]]
#>    cov1.1    cov2.1    cov3.1 
#> 0.4946542 0.7026185 0.1710491 
#> 
#> [[2]]
#>    cov1.2    cov2.2    cov3.2 
#> 0.4120727 0.6829041 0.2469058 
#> 
#> [[3]]
#>    cov1.3    cov2.3    cov3.3 
#> 0.5812089 0.2694629 0.2581858
#> Iteration 5 — Total criteria: 0.00051 (c1: 0.00029, c2: 0.00011,  c3: 0.00011)
```

See model results:

``` r
print(fit[[1]]) # trans 1 
#> Pooled flexsurvreg (Rubin's rules)
#> m imputations: 5 
#> 
#>            est     se       df      LCL      UCL
#> shape   0.1100 0.0048 427.7631   0.1006   0.1195
#> rate  -13.1391 0.4492 275.1557 -14.0234 -12.2548
#> cov1    0.5129 0.4075  34.7012  -0.3146   1.3404
#> cov2    0.6975 0.1243  24.1002   0.4411   0.9540
#> cov3    0.1671 0.0926  50.5372  -0.0189   0.3531

print(fit[[2]]) # trans 2
#> Pooled flexsurvreg (Rubin's rules)
#> m imputations: 5 
#> 
#>            est     se         df      LCL      UCL
#> shape   0.1171 0.0026 26697.5933   0.1121   0.1221
#> rate  -12.5334 0.2370 15008.0791 -12.9979 -12.0688
#> cov1    0.4040 0.1857  2514.6780   0.0398   0.7681
#> cov2    0.6825 0.0555   502.8034   0.5735   0.7915
#> cov3    0.2474 0.0434  1901.6366   0.1623   0.3324

print(fit[[3]]) # trans 3
#> Pooled flexsurvreg (Rubin's rules)
#> m imputations: 5 
#> 
#>           est     se        df      LCL     UCL
#> shape  0.0804 0.0064 1136.0075   0.0677  0.0930
#> rate  -8.8907 0.6101 1207.3289 -10.0876 -7.6937
#> cov1   0.5891 0.3801  156.3696  -0.1617  1.3399
#> cov2   0.2395 0.1281   54.7036  -0.0172  0.4962
#> cov3   0.2527 0.0962  133.0243   0.0625  0.4429
```
