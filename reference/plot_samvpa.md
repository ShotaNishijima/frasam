# Plot SAM or VPA results

Plot SAM or VPA results

## Usage

``` r
plot_samvpa(
  vpa_sam_list,
  CI = 0.95,
  scenario_name = NULL,
  alpha = 0.4,
  size = 1,
  base_size = 14,
  log_scale = FALSE,
  legend_name = "Scenario",
  legend_nrow = 1,
  legend_position = "top",
  what.plot = c("biomass", "SSB", "Recruitment", "U"),
  years = NULL,
  ncol = 2,
  scale_recruitment = 1000,
  scale_biomass = 1000,
  scale_ssb = 1000,
  scale_catch = 1000
)
```

## Arguments

- vpa_sam_list:

  A SAM or VPA result object, or a list of result objects.

- CI:

  Confidence interval width. Set 0 to omit confidence intervals.

- scenario_name:

  Scenario names used in the legend.

- alpha:

  Alpha value for confidence interval ribbons.

- size:

  Line width.

- base_size:

  Base font size.

- log_scale:

  If `TRUE`, use a log scale on the y-axis.

- legend_name:

  Legend title.

- legend_nrow:

  Number of rows in the legend.

- legend_position:

  Legend position.

- what.plot:

  Statistics to plot.

- years:

  Years to include in the plot.

- ncol:

  Number of columns in the facet plot.

- scale_recruitment:

  Divisor for recruitment values in the plot.

- scale_biomass:

  Divisor for biomass values in the plot.

- scale_ssb:

  Divisor for spawning stock biomass values in the plot.

- scale_catch:

  Divisor for catch values in the plot.
