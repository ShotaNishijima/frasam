#' Convert vector, matrix, or data.frame to age-year matrix
#'
#' @keywords internal
.as_age_year_matrix <- function(x, target, name = "M") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.matrix(x)) {
    storage.mode(x) <- "numeric"

    if (!identical(dim(x), dim(target))) {
      stop(name, " must have the same dimensions as naa.", call. = FALSE)
    }

    rownames(x) <- rownames(target)
    colnames(x) <- colnames(target)

    return(x)
  }

  if (is.vector(x) && length(x) == nrow(target)) {
    out <- matrix(x, nrow = nrow(target), ncol = ncol(target))
    storage.mode(out) <- "numeric"
    rownames(out) <- rownames(target)
    colnames(out) <- colnames(target)
    return(out)
  }

  stop(
    name,
    " must be either a matrix/data.frame with the same dimensions as naa ",
    "or a vector of length n_age.",
    call. = FALSE
  )
}

#' Generate all permutations of a character vector
#'
#' @keywords internal
.permutations <- function(x) {
  if (length(x) == 1L) {
    return(matrix(x, nrow = 1L))
  }

  out <- lapply(seq_along(x), function(i) {
    cbind(x[i], .permutations(x[-i]))
  })

  do.call(rbind, out)
}


#' Safe ratio
#'
#' @keywords internal
.safe_ratio <- function(num, den, zero_tol = 1e-12) {
  if (anyNA(c(num, den))) {
    return(NA_real_)
  }

  if (abs(den) < zero_tol) {
    if (abs(num) < zero_tol) {
      return(1)
    } else {
      return(NA_real_)
    }
  }

  num / den
}


#' Shapley decomposition for one transition
#'
#' @keywords internal
.shapley_transition <- function(N_from,
                                w_from,
                                w_to,
                                F_from,
                                M_from,
                                process_multiplier,
                                m_from = NULL,
                                m_to = NULL) {
  use_maturity <- !is.null(m_from) && !is.null(m_to)

  values_to_check <- c(N_from, w_from, w_to, F_from, M_from, process_multiplier)
  if (use_maturity) {
    values_to_check <- c(values_to_check, m_from, m_to)
  }

  effect_names <- c("growth", "fishing", "natural", "process")
  if (use_maturity) {
    effect_names <- c("growth", "maturity", "fishing", "natural", "process")
  }

  if (anyNA(values_to_check)) {
    out <- setNames(rep(NA_real_, length(effect_names)), effect_names)
    return(out)
  }

  if (N_from == 0) {
    out <- setNames(rep(0, length(effect_names)), effect_names)
    return(out)
  }

  if (w_from <= 0 || w_to <= 0) {
    stop("Body weights must be positive.", call. = FALSE)
  }

  if (use_maturity && (m_from < 0 || m_to < 0)) {
    stop("Maturity values must be non-negative.", call. = FALSE)
  }

  value <- function(active) {
    w <- if ("growth" %in% active) w_to else w_from

    m <- 1
    if (use_maturity) {
      m <- if ("maturity" %in% active) m_to else m_from
    }

    f_mult <- if ("fishing" %in% active) exp(-F_from) else 1
    m_mult <- if ("natural" %in% active) exp(-M_from) else 1
    p_mult <- if ("process" %in% active) process_multiplier else 1

    N_from * w * m * f_mult * m_mult * p_mult
  }

  perms <- .permutations(effect_names)
  contrib <- setNames(rep(0, length(effect_names)), effect_names)

  for (i in seq_len(nrow(perms))) {
    active <- character(0)

    for (v in perms[i, ]) {
      before <- value(active)
      active <- c(active, v)
      after <- value(active)

      contrib[v] <- contrib[v] + after - before
    }
  }

  contrib / nrow(perms)
}


#' Decompose biomass or spawning-stock-biomass changes by age
#'
#' @param naa Matrix of numbers at age. Rows are ages, columns are years.
#' @param waa Matrix of body weights at age. Rows are ages, columns are years.
#' @param faa Matrix of fishing mortality at age. Rows are ages, columns are years.
#' @param M Matrix or vector of natural mortality at age.
#' @param maa Optional matrix of maturity at age. If NULL, biomass is decomposed.
#'   If supplied, spawning-stock biomass is decomposed.
#' @param plus_group Logical. If TRUE, the last row is treated as a plus group.
#' @param recruitment_age_row Row index for recruitment. Default is 1, assumed to be age 0.
#' @param zero_tol Tolerance for zero denominators.
#'
#' @return A list of matrices with the same dimensions as naa.
#'
#' @export
decompose_biomass_effects <- function(naa,
                                      waa,
                                      faa,
                                      M,
                                      maa = NULL,
                                      plus_group = TRUE,
                                      recruitment_age_row = 1L,
                                      zero_tol = 1e-12) {
  naa <- .as_numeric_matrix(naa, name = "naa")
  waa <- .as_numeric_matrix(waa, name = "waa")
  faa <- .as_numeric_matrix(faa, name = "faa")

  if (!identical(dim(naa), dim(waa)) || !identical(dim(naa), dim(faa))) {
    stop("naa, waa, and faa must have the same dimensions.", call. = FALSE)
  }

  n_age <- nrow(naa)
  n_year <- ncol(naa)

  if (n_age < 2L) {
    stop("naa must have at least two age rows.", call. = FALSE)
  }

  if (n_year < 2L) {
    stop("naa must have at least two year columns.", call. = FALSE)
  }

  if (recruitment_age_row < 1L || recruitment_age_row > n_age) {
    stop("recruitment_age_row is out of range.", call. = FALSE)
  }

  M <- .as_age_year_matrix(M, target = naa, name = "M")

  use_maturity <- !is.null(maa)

  if (use_maturity) {
    maa <- .as_numeric_matrix(maa, name = "maa")

    if (!identical(dim(maa), dim(naa))) {
      stop("maa must have the same dimensions as naa.", call. = FALSE)
    }

    rownames(maa) <- rownames(naa)
    colnames(maa) <- colnames(naa)

    if (any(maa < 0, na.rm = TRUE)) {
      stop("maa must be non-negative.", call. = FALSE)
    }
  }

  if (!identical(dim(naa), dim(waa)) || !identical(dim(naa), dim(faa))) {
    stop("naa, waa, and faa must have the same dimensions.", call. = FALSE)
  }

  make_mat <- function(value = 0) {
    out <- matrix(
      as.numeric(value),
      nrow = n_age,
      ncol = n_year
    )
    rownames(out) <- rownames(naa)
    colnames(out) <- colnames(naa)
    storage.mode(out) <- "numeric"
    out
  }

  recruitment <- make_mat(0)
  growth <- make_mat(0)
  fishing <- make_mat(0)
  natural <- make_mat(0)
  process <- make_mat(0)
  terminal_loss <- make_mat(0)

  maturity <- NULL
  if (use_maturity) {
    maturity <- make_mat(0)
  }

  target_quantity <- if (use_maturity) {
    maa * waa * naa
  } else {
    waa * naa
  }

  recruitment[recruitment_age_row, ] <- target_quantity[recruitment_age_row, ]

  add_contrib <- function(dest_age, year, contrib) {
    growth[dest_age, year]  <<- growth[dest_age, year]  + contrib["growth"]
    fishing[dest_age, year] <<- fishing[dest_age, year] + contrib["fishing"]
    natural[dest_age, year] <<- natural[dest_age, year] + contrib["natural"]
    process[dest_age, year] <<- process[dest_age, year] + contrib["process"]

    if (use_maturity) {
      maturity[dest_age, year] <<-
        maturity[dest_age, year] + contrib["maturity"]
    }
  }

  decompose_one_transition <- function(a_from, a_to, y, process_multiplier) {
    y_prev <- y - 1L

    if (use_maturity) {
      .shapley_transition(
        N_from = naa[a_from, y_prev],
        w_from = waa[a_from, y_prev],
        w_to = waa[a_to, y],
        F_from = faa[a_from, y_prev],
        M_from = M[a_from, y_prev],
        process_multiplier = process_multiplier,
        m_from = maa[a_from, y_prev],
        m_to = maa[a_to, y]
      )
    } else {
      .shapley_transition(
        N_from = naa[a_from, y_prev],
        w_from = waa[a_from, y_prev],
        w_to = waa[a_to, y],
        F_from = faa[a_from, y_prev],
        M_from = M[a_from, y_prev],
        process_multiplier = process_multiplier
      )
    }
  }

  for (y in 2:n_year) {
    y_prev <- y - 1L

    if (plus_group) {
      normal_dest_ages <- if (n_age > 2L) 2:(n_age - 1L) else integer(0)
    } else {
      normal_dest_ages <- 2:n_age
    }

    for (a_to in normal_dest_ages) {
      a_from <- a_to - 1L

      pred_no_process <-
        naa[a_from, y_prev] * exp(-faa[a_from, y_prev] - M[a_from, y_prev])

      process_multiplier <-
        .safe_ratio(naa[a_to, y], pred_no_process, zero_tol = zero_tol)

      contrib <- decompose_one_transition(
        a_from = a_from,
        a_to = a_to,
        y = y,
        process_multiplier = process_multiplier
      )

      add_contrib(dest_age = a_to, year = y, contrib = contrib)
    }

    if (plus_group) {
      a_plus <- n_age
      a_from_1 <- n_age - 1L
      a_from_2 <- n_age

      pred_1 <-
        naa[a_from_1, y_prev] *
        exp(-faa[a_from_1, y_prev] - M[a_from_1, y_prev])

      pred_2 <-
        naa[a_from_2, y_prev] *
        exp(-faa[a_from_2, y_prev] - M[a_from_2, y_prev])

      pred_total <- pred_1 + pred_2

      process_multiplier <-
        .safe_ratio(naa[a_plus, y], pred_total, zero_tol = zero_tol)

      contrib_1 <- decompose_one_transition(
        a_from = a_from_1,
        a_to = a_plus,
        y = y,
        process_multiplier = process_multiplier
      )

      contrib_2 <- decompose_one_transition(
        a_from = a_from_2,
        a_to = a_plus,
        y = y,
        process_multiplier = process_multiplier
      )

      add_contrib(dest_age = a_plus, year = y, contrib = contrib_1)
      add_contrib(dest_age = a_plus, year = y, contrib = contrib_2)
    }

    if (!plus_group) {
      terminal_loss[n_age, y] <- -target_quantity[n_age, y_prev]
    }
  }

  effect_sum_mat <- recruitment + growth + fishing + natural + process
  if (use_maturity) {
    effect_sum_mat <- effect_sum_mat + maturity
  }

  annual_effect_sum <- colSums(effect_sum_mat, na.rm = FALSE)

  annual_effect_sum_with_terminal <-
    colSums(effect_sum_mat + terminal_loss, na.rm = FALSE)

  annual_change <- rep(NA_real_, n_year)
  annual_change[-1] <-
    colSums(target_quantity, na.rm = FALSE)[-1] -
    colSums(target_quantity, na.rm = FALSE)[-n_year]

  names(annual_change) <- colnames(naa)

  out <- list(
    recruitment = recruitment,
    growth = growth,
    fishing = fishing,
    natural = natural,
    process = process,
    terminal_loss = terminal_loss,
    target_quantity = target_quantity,
    annual_effect_sum = annual_effect_sum,
    annual_effect_sum_with_terminal = annual_effect_sum_with_terminal,
    annual_change = annual_change,
    residual = annual_effect_sum - annual_change,
    residual_with_terminal = annual_effect_sum_with_terminal - annual_change,
    target = if (use_maturity) "ssb" else "biomass"
  )

  if (use_maturity) {
    out$maturity <- maturity
  }

  .add_decomposition_summaries(out)
}

#' Add age-aggregated and percentage summaries to a decomposition
#'
#' @keywords internal
.add_decomposition_summaries <- function(x) {
  effect_names <- c(
    "recruitment", "growth", "fishing", "process", "natural", "maturity"
  )
  template <- x$recruitment
  effect_matrices <- setNames(lapply(effect_names, function(effect) {
    if (!is.null(x[[effect]])) x[[effect]] else template * 0
  }), effect_names)
  age_aggregated <- do.call(rbind, lapply(effect_matrices, function(z) {
    colSums(z, na.rm = FALSE)
  }))
  rownames(age_aggregated) <- effect_names
  colnames(age_aggregated) <- colnames(template)
  previous_total <- c(
    NA_real_,
    colSums(x$target_quantity, na.rm = FALSE)[-ncol(x$target_quantity)]
  )
  names(previous_total) <- colnames(template)
  percent_by_age <- lapply(effect_matrices, function(z) {
    out <- sweep(z, 2L, previous_total, "/") * 100
    out[, !is.finite(previous_total)] <- NA_real_
    out
  })
  percent_aggregated <- sweep(age_aggregated, 2L, previous_total, "/") * 100
  percent_aggregated[, !is.finite(previous_total)] <- NA_real_
  x$age_aggregated <- age_aggregated
  x$percent_by_age <- percent_by_age
  x$percent_aggregated <- percent_aggregated
  x$previous_total <- previous_total
  x
}

#' Decompose biomass changes from a VPA or SAM result
#'
#' Extracts numbers at age, fishing mortality, weight, natural mortality,
#' maturity, and the plus-group setting from a fitted `vpa` or `sam` object,
#' then calls [decompose_biomass_effects()].
#'
#' @param result A fitted object of class `vpa` or `sam`.
#' @param target Quantity to decompose: `"biomass"` or `"ssb"`.
#' @param recruitment_age_row Row containing recruitment.
#' @param zero_tol Tolerance for zero denominators.
#'
#' @return The result of [decompose_biomass_effects()], including age-specific,
#'   age-aggregated, and percentage contributions. Percentages use the total
#'   biomass or SSB in the preceding year as denominator.
#'
#' @export
decompose_biomass_factors <- function(
    result,
    target = c("biomass", "ssb"),
    recruitment_age_row = 1L,
    zero_tol = 1e-12) {
  target <- match.arg(target)
  if (length(intersect(class(result), c("vpa", "sam"))) == 0L) {
    stop("result must be an object of class 'vpa' or 'sam'.", call. = FALSE)
  }
  dat <- result$input$dat
  if (is.null(dat)) stop("result$input$dat is missing.", call. = FALSE)
  required_result <- c("naa", "faa")
  missing_result <- required_result[vapply(
    required_result, function(z) is.null(result[[z]]), logical(1)
  )]
  if (length(missing_result) > 0L) {
    stop("Missing result component(s): ", paste(missing_result, collapse = ", "),
         call. = FALSE)
  }
  required_dat <- c("waa", "M")
  if (target == "ssb") required_dat <- c(required_dat, "maa")
  missing_dat <- required_dat[vapply(
    required_dat, function(z) is.null(dat[[z]]), logical(1)
  )]
  if (length(missing_dat) > 0L) {
    stop("Missing result$input$dat component(s): ",
         paste(missing_dat, collapse = ", "), call. = FALSE)
  }
  waa <- dat$waa
  if (inherits(result, "sam") && !is.null(result$waa_est) &&
      identical(dim(result$waa_est), dim(result$naa))) {
    waa <- result$waa_est
  }
  plus_group <- result$input$plus.group
  if (is.null(plus_group) && !is.null(result$data$maxAgePlusGroup)) {
    plus_group <- as.logical(result$data$maxAgePlusGroup[1])
  }
  if (is.null(plus_group) || length(plus_group) != 1L || is.na(plus_group)) {
    stop("The plus-group setting could not be determined from result.", call. = FALSE)
  }
  decompose_biomass_effects(
    naa = result$naa, waa = waa, faa = result$faa, M = dat$M,
    maa = if (target == "ssb") dat$maa else NULL,
    plus_group = isTRUE(plus_group),
    recruitment_age_row = recruitment_age_row, zero_tol = zero_tol
  )
}

