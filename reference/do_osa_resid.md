# SAMのOSA residualを計算する関数

SAMのOSA residualを計算する関数

## Usage

``` r
do_osa_resid(
  samres,
  method = "oneStepGaussianOffMode",
  subset = NULL,
  trace = 2
)
```

## Arguments

- samres:

  sam object

- method:

  Method to calculate OSA (see details).

- subset:

  Index vector of observations that will be added one by one during OSA.
  By default `1:length(observations)` (with `conditional` subtracted).

- trace:

  Logical; Trace progress? More options available for
  `method="oneStepGeneric"` - see details.
