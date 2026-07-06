# Plot age-aggregated biomass factors

Draws a stacked bar chart of the age-aggregated contributions returned
by \[decompose_biomass_factors()\]. Positive and negative contributions
are stacked on opposite sides of zero.

## Usage

``` r
plot_biomass_factors(x, type = c("percent", "absolute"), scale = 1)
```

## Arguments

- x:

  Output from \[decompose_biomass_factors()\].

- type:

  Output scale. \`"percent"\` (default) uses \`percent_aggregated\`;
  \`"absolute"\` uses \`age_aggregated\`.

- scale:

  Positive divisor applied to values when \`type = "absolute"\`.

## Value

A \`ggplot\` object. The \`effect\` variable in the plot data is a
factor whose levels follow the row order of the aggregated matrix.
