make_short_simulated_future_data <- function(n_sim = 2) {
  sam_bh_prec <- readRDS(testthat::test_path("testdata", "sam_bh_prec.rds"))
  use_sam_tmb(TmbFile = sam_bh_prec$input$cpp.file.name, overwrite = FALSE)

  mu <- sam_bh_prec$obj$env$last.par.best
  precision <- sam_bh_prec$rep$jointPrecision
  par_sim <- rmvnorm_prec(mu, precision, n.sims = n_sim, seed = 2)

  sam_sim <- lapply(seq_len(n_sim), function(i) {
    suppressWarnings(update_sam(sam_bh_prec, par_sim[, i]))
  })

  faa_current <- rowMeans(
    sam_bh_prec$faa[, as.character(2011:2013), drop = FALSE]
  )

  future_list <- lapply(seq_along(sam_sim), function(i) {
    x <- sam_sim[[i]]
    invisible(capture.output(
      result <- frasyr::make_future_data(
        x,
        nsim = 1,
        nyear = 3,
        future_initial_year_name = 2014,
        start_F_year_name = 2014,
        start_biopar_year_name = 2015,
        start_random_rec_year_name = 2015,
        waa_year = 2014,
        waa_catch_year = 2014,
        maa_year = 2014,
        M_year = 2014,
        faa_year = NULL,
        currentF = faa_current,
        futureF = faa_current,
        start_ABC_year_name = 2015,
        Pope = FALSE,
        res_SR = make_SRres(x),
        more_process_error = x$sigma.logN[-1],
        scale_ssb = 1 / x$input$scale,
        scale_R = x$input$scale_number,
        seed_number = 20260706 + i,
        bias_correction = FALSE
      )
    ))
    result
  })

  list(
    par_sim = par_sim,
    sam_sim = sam_sim,
    future_list = future_list,
    future = frasyr::unlist_future_data(future_list)
  )
}