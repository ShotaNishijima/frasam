library(frasam)

test_that("plot_cond_nll draws horizontal bars and value labels", {
  x <- tibble::tibble(
    type = factor(c("Process_N", "Process_F", "Survey")),
    nll = c(-12.34, -2.01, 8.26)
  )

  p <- plot_cond_nll(x)

  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$data$type), as.character(x$type))
  expect_equal(p$data$label, c("-12.3", "-2.0", "8.3"))
  expect_equal(p$data$label_hjust, c(1.1, 1.1, -0.1))
  expect_equal(length(p$layers), 3L)
})

test_that("plot_cond_nll can omit labels and validates input", {
  x <- data.frame(type = c("Process_N", "Survey"), nll = c(-1, 2))

  p <- plot_cond_nll(x, show_value = FALSE)
  expect_equal(length(p$layers), 2L)

  expect_error(plot_cond_nll(data.frame(type = "x")), "get_cond_nll")
  expect_error(plot_cond_nll(x, show_value = NA), "show_value")
})
