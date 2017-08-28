
<!-- README.md is generated from README.Rmd. Please edit that file -->
ri2 makes conducting randomization inference easy and (with the blessing of the original authors) is the successor package to [ri](https://cran.r-project.org/web/packages/ri/index.html).

ri2 has specific support for the following:

1.  All randomization schemes in [randomizr](http://randomizr.declaredesign.org).
2.  Difference-in-means and OLS-adjusted estimates of ATE estimates using R-native formula syntax.
3.  Multi-arm trials.
4.  ANOVA-style hypothesis tests (e.g., testing interaction term under null of constant effects),

Additionally, ri2 provides:

1.  Accommodation for arbitrary randomization schemes
2.  Accommodation for arbitrary (scalar) test statistics

If you'd like to install the most current development release, you can do so using the `devtools` package. Since there's some `Rcpp` in `ri2`, you'll need some development tools as described on the [devtools page](https://www.rstudio.com/products/rpackages/devtools/).

-   on a Mac, you have to first install Xcode (from the app store is easiest!)
-   on Windows, install `Rtools` with `install.packages("Rtools")`.

Then use the following code:

``` r
install.packages("devtools")
devtools::install_github("DeclareDesign/randomizr")
devtools::install_github("acoppock/ri2")
```

Here is the basic syntax for a two-arm trial:

``` r
library(ri2)
#> Loading required package: randomizr
#> Loading required package: estimatr
N <- 100
declaration <- declare_ra(N = N, m = 50)

Z <- conduct_ra(declaration)
X <- rnorm(N)
Y <- .9 * X + .2 * Z + rnorm(N)
df <- data.frame(Y, X, Z)

ri_out <-
  conduct_ri(
    formula = Y ~ Z,
    declaration = declaration,
    assignment = "Z",
    sharp_hypothesis = 0,
    data = df
  )

plot(ri_out)
```

![](README-unnamed-chunk-3-1.png)

``` r
summary(ri_out)
#> # A tibble: 1 x 5
#>   coefficient    estimate p_value null_ci_lower null_ci_upper
#>         <chr>       <dbl>   <dbl>         <dbl>         <dbl>
#> 1           Z -0.02348029   0.936    -0.6029264     0.5557554
```

The development of ri2 is supported by a Standards Grant from [EGAP](http://egap.org).
