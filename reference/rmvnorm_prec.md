# Simulate multivariate normal variables given a mean vector and precision matrix

Simulate multivariate normal variables given a mean vector and precision
matrix

## Usage

``` r
rmvnorm_prec(mu, prec, n.sims, seed = 123)
```

## Arguments

- mu:

  vector of parameter means

- prec:

  joint precision matrix

- n.sims:

  number of draws

- seed:

  seed number

## Value

length(mu) by n.sims matrix of parameter draws
