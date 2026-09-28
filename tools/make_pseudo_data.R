load_local_frasam <- function() {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all()
    return(invisible(TRUE))
  }

  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all()
    return(invisible(TRUE))
  }

  r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)
  for (r_file in r_files) {
    source(r_file, local = globalenv())
  }
  invisible(TRUE)
}

if (!requireNamespace("frasyr", quietly = TRUE)) {
  stop("Package 'frasyr' is required to run this script.")
}

library(frasyr)

load_local_frasam()

res_sam <- get(load("tools/sam_masaba_P2025.rda"))
use_sam_tmb(TmbFile = "sam2", overwrite = FALSE)

input <- res_sam$input
input$dat$maa[1:7, ] <- c(0, 0, 0.5, 1, 1, 1, 1)
input$dat$waa[1:7, ] <- c(100, 200, 280, 350, 380, 400, 400)
input$dat$M[] <- 0.6
# input$abund[6:7] <- "B"
# input$p0.list <- NULL
# input$AR <- 0
# input$SR <- "RW"
# input$catch_prop <- NULL


use_yrs <- as.character(1975:2014)

input$dat$caa <- input$dat$caa[,use_yrs]
input$dat$waa <- input$dat$waa[,use_yrs]
input$dat$M <- input$dat$M[,use_yrs]
input$dat$index <- input$dat$index[1:5,use_yrs]
input$dat$maa <- input$dat$maa[,use_yrs]
input$abund <- c("N","N","N","SSB","B")
input$b.fix <- input$b.fix[1:5]
input$min.age <- input$min.age[1:5]
input$max.age <- input$max.age[1:5]
input$catch_prop <- NULL
# input$dat$catch.prop

dat <- input$dat

res_tmp <- sam(
  dat,
  last.catch.zero = TRUE,
  abund = c("N","N","N","SSB","B"),
  min.age=c(0,0,1,0,0),
  max.age = c(0,0,1,6,6),
  rec.age = 0,
  index.key=0:4,
  b.est=FALSE,
  SR = "RW",
  varC = c(0,0,1,1,1,2,2),
  varF = c(0,0,1,1,1,1,1),
  varN = c(0,1,1,1,1,1,1),
  rho.mode=3,
  bias.correct = FALSE,
  silent = TRUE
)

# res_tmp$opt
# res_tmp$rep
# OK
# res_tmp <- do.call(sam, input)


make_reverse_year_map <- function(dat,
                                  years = 1975:2014,
                                  include_other_years = TRUE) {
  dat_years <- colnames(dat$caa)
  if (is.null(dat_years)) {
    stop("dat$caa must have year column names.")
  }

  target_col <- match(as.character(years), dat_years)
  source_col <- match(as.character(rev(years)), dat_years)
  keep <- !is.na(target_col) & !is.na(source_col)

  if (!any(keep)) {
    stop("None of the requested years were found in dat$caa.")
  }

  if (!include_other_years) {
    return(data.frame(
      target_year = years[keep],
      source_year = rev(years)[keep],
      target_col = target_col[keep],
      source_col = source_col[keep]
    ))
  }

  all_target_col <- seq_along(dat_years)
  all_source_col <- all_target_col
  all_source_year <- dat_years
  all_source_col[target_col[keep]] <- source_col[keep]
  all_source_year[target_col[keep]] <- as.character(rev(years)[keep])

  data.frame(
    target_year = dat_years,
    source_year = all_source_year,
    target_col = all_target_col,
    source_col = all_source_col
  )
}

simulate_from_sam <- function(res,
                              seed = NULL,
                              index_error_scale = 2.5,
                              caa_error_scale = 2.5,
                              min_index_sd = 0.25,
                              min_caa_sd = 0.20,
                              value_floor = 1e-12,
                              year_map = NULL,
                              fill_index_na = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  dat <- res$input$dat

  pred_index <- as.matrix(res$pred.index)
  pred_caa <- as.matrix(res$caa)

  sigma_index <- res$sigma
  if (res$input$est.method == "ls" && is.null(res$input$index.key)) {
    sigma_index <- rep(res$sigma, nrow(pred_index))
  }
  sigma_index <- rep_len(as.numeric(sigma_index), nrow(pred_index))
  sigma_index <- pmax(sigma_index * index_error_scale, min_index_sd)

  sigma_caa <- rep_len(as.numeric(res$sigma.logC), nrow(pred_caa))
  sigma_caa <- pmax(sigma_caa * caa_error_scale, min_caa_sd)

  if (is.null(year_map)) {
    year_map <- data.frame(
      target_col = seq_len(ncol(pred_caa)),
      source_col = seq_len(ncol(pred_caa))
    )
  }

  for (a in seq_len(nrow(pred_caa))) {
    dat$caa[a, year_map$target_col] <- exp(rnorm(
      n = nrow(year_map),
      mean = log(pmax(pred_caa[a, year_map$source_col], value_floor)),
      sd = sigma_caa[a]
    ))
  }

  if (isTRUE(res$input$last.catch.zero)) {
    dat$caa[, ncol(dat$caa)] <- 0
  }

  index_target_col <- year_map$target_col
  index_source_col <- year_map$source_col
  index_keep <- index_target_col <= ncol(dat$index) &
    index_source_col <= ncol(pred_index)
  index_target_col <- index_target_col[index_keep]
  index_source_col <- index_source_col[index_keep]

  for (i in seq_len(nrow(pred_index))) {
    if (fill_index_na) {
      ok <- !is.na(pred_index[i, index_source_col])
    } else {
      ok <- !is.na(dat$index[i, index_target_col]) &
        !is.na(pred_index[i, index_source_col])
    }

    if (!any(ok)) next

    dat$index[i, index_target_col[ok]] <- exp(rnorm(
      n = sum(ok),
      mean = log(pmax(pred_index[i, index_source_col[ok]], value_floor)),
      sd = sigma_index[i]
    ))
  }

  dat
}

fit_pseudo_sam <- function(generator_res,
                           pseudo_dat,
                           use_previous_par = TRUE,
                           silent = TRUE) {
  input <- generator_res$input
  input$dat <- pseudo_dat
  input$silent <- silent

  if (use_previous_par) {
    input$p0.list <- generator_res$par_list
  } else {
    input$p0.list <- NULL
  }

  do.call(sam, input)
}

make_iterative_pseudo_sam <- function(initial_res,
                                      n_iter = 10,
                                      seed = 123,
                                      index_error_scale = 2.5,
                                      caa_error_scale = 2.5,
                                      use_previous_par = TRUE,
                                      first_year_map = NULL,
                                      fill_index_na = TRUE) {
  pseudo_dat_list <- vector("list", n_iter)
  pseudo_res_list <- vector("list", n_iter)
  status <- data.frame(
    iter = seq_len(n_iter),
    seed = seed + seq_len(n_iter) - 1L,
    ok = FALSE,
    loglik = NA_real_,
    aic = NA_real_,
    convergence = NA_integer_,
    message = NA_character_,
    stringsAsFactors = FALSE
  )

  generator_res <- initial_res

  for (i in seq_len(n_iter)) {
    iter_seed <- seed + i - 1L

    pseudo_dat <- simulate_from_sam(
      generator_res,
      seed = iter_seed,
      index_error_scale = index_error_scale,
      caa_error_scale = caa_error_scale,
      year_map = if (i == 1L) first_year_map else NULL,
      fill_index_na = fill_index_na
    )
    pseudo_dat_list[[i]] <- pseudo_dat

    fit <- try(
      fit_pseudo_sam(
        generator_res = generator_res,
        pseudo_dat = pseudo_dat,
        use_previous_par = use_previous_par
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      status$message[i] <- as.character(fit)
      next
    }

    pseudo_res_list[[i]] <- fit
    status$ok[i] <- TRUE
    status$loglik[i] <- fit$loglik
    status$aic[i] <- fit$aic
    status$convergence[i] <- fit$opt$convergence

    generator_res <- fit
  }

  list(
    initial_res = initial_res,
    pseudo_dat = pseudo_dat_list,
    pseudo_res = pseudo_res_list,
    status = status
  )
}

reverse_year_map <- make_reverse_year_map(res_tmp$input$dat, years = 1975:2014)

pseudo_masaba <- make_iterative_pseudo_sam(
  initial_res = res_tmp,
  n_iter = 10,
  seed = 20250701,
  index_error_scale = 1,
  caa_error_scale = 1,
  use_previous_par = TRUE,
  first_year_map = reverse_year_map,
  fill_index_na = TRUE
)

print(pseudo_masaba$status)

res_tmp$baa %>% colSums %>% plot


pseudo_masaba$pseudo_res[[3]]$baa %>% colSums %>% plot
pseudo_masaba$pseudo_res[[6]]$baa %>% colSums %>% plot

pseudo_masaba$pseudo_res[[4]]$ssb %>% colSums %>% plot
pseudo_masaba$pseudo_res[[10]]$ssb %>% colSums %>% plot

i <- 6
plot(
  colSums(pseudo_masaba$pseudo_res[[i]]$ssb),
  pseudo_masaba$pseudo_res[[i]]$naa[1,],
  # log = "y"
)

res_tmp$input$dat$caa %>% colSums %>% plot
# pseudo_masaba$pseudo_res[[1]]$input$dat$caa %>% colSums %>% plot

# pseudo_masaba$pseudo_res[[1]]$input$dat$index
# input$min.age

save(
  res_tmp,
  pseudo_masaba,
  file = "tools/pseudo_masaba_iterative.rda"
)

i <- 6 #8番目の結果を使う
sam_ex <- pseudo_masaba$pseudo_res[[i]]
dat_ex <- pseudo_masaba$pseudo_res[[i]]$input$dat

save(sam_ex, file = "data/sam_ex.rda")
save(dat_ex, file = "data/dat_ex.rda")

data("dat_ex")
data(sam_ex)

# p0_list <- sam_ex$par_list

input_ex <- sam_ex$input

# res_test <- do.call(sam, input_ex)
# input_ex$p0.list
res_test  <- sam(
  dat_ex,
  last.catch.zero = TRUE,
  abund = c("N","N","N","SSB","B"),
  min.age=c(0,0,1,0,0),
  max.age = c(0,0,1,6,6),
  rec.age = 0,
  index.key=0:4,
  b.est=FALSE,
  SR = "RW",
  varC = c(0,0,1,1,1,2,2),
  varF = c(0,0,1,1,1,1,1),
  varN = c(0,1,1,1,1,1,1),
  rho.mode=3,
  bias.correct = FALSE,
  silent = TRUE
)

res_test$opt

c(res_test$loglik,sam_ex$loglik)

retro_ex = retro_sam(sam_ex,n=3)
# retro_ex$mohn

save(retro_ex, file = "data/retro_ex.rda")

?do_osa_resid()
tmp = do_osa_resid(sam_ex,subset=101:103)
