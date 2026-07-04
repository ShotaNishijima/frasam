# Update a future projection plot legend for SAM

\`frasyr::plot_futures()\` labels the historical assessment line as VPA.
This helper updates only the legend labels of the returned ggplot object
so the historical line is shown as SAM when a SAM result was supplied.

## Usage

``` r
plot_update2sam(
  plot,
  sam_label = "SAM",
  scenario_labels = NULL,
  legend_title = "",
  ncol_legend = 2
)
```

## Arguments

- plot:

  A ggplot object returned by
  [`frasyr::plot_futures()`](https://rdrr.io/pkg/frasyr/man/plot_futures.html).

- sam_label:

  Label used for the historical SAM line.

- scenario_labels:

  Optional labels for future scenarios. Supply either a named character
  vector whose names match the current scenario names, or an unnamed
  vector with the same length as the non-SAM scenarios.

- legend_title:

  Legend title.

- ncol_legend:

  Number of columns in the colour legend.

## Value

A ggplot object with updated colour and fill legend labels.
