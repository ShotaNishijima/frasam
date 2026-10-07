#######################################
##### Parametric bootstrap of SAM #####
#######################################

#'
#' @export
#'

boo_sam <- function(Res,n=100,seed=1,est=TRUE,method="p",use_p0=TRUE){

  if(is.numeric(set.seed)) set.seed(set.seed)

  Res$input$bias.correct.sd <- FALSE
  Res0 <- Res
  # if (Res$input$bias.correct) {
  #   Res$input$bias.correct <- FALSE
  #   Res <- do.call(sam, Res$input)
  # }

  if (isTRUE(est)) {
    boot.list <- list()
    pred.index <- Res$pred.index
    catch.sim <- sam_catch_simulation_inputs(Res)
    pred.caa <- catch.sim$mean
    sigma.index <- Res$sigma
    sigma.caa <- catch.sim$sd
    catch.present <- catch.sim$present
    resid.index <- log(as.matrix(Res$input$dat$index))-log(as.matrix(Res$pred.index))
    resid.caa <- log(as.matrix(Res$input$dat$caa))-log(as.matrix(pred.caa))
    if (Res$input$last.catch.zero) resid.caa <- resid.caa[, -ncol(resid.caa), drop = FALSE]
    for (j in 1:n) {
      sim.dat <- Res$input$dat
      for(k in 1:100) {
        if (method == "p") { # parametric bootstrap
          for (i in 1:nrow(pred.index)) {
            sim.index <- exp(rnorm(ncol(pred.index),as.numeric(log(pred.index[i,])),sigma.index[i]))
            sim.dat$index[i,!is.na(sim.dat$index[i,])] <- sim.index[!is.na(sim.dat$index[i,])]
          }
          for (i in 1:nrow(pred.caa)) {
            present <- catch.present[i, ]
            sim.dat$caa[i,present] <- exp(rnorm(sum(present),
              log(pred.caa[i,present]),sigma.caa[i,present]))
          }
          if (Res$input$last.catch.zero) sim.dat$caa[,ncol(pred.caa)] <- 0
        }

        if (method == "n") { #non-parametric bootstrap
          for (i in 1:nrow(resid.index)) {
            sim.dat$index[i,!is.na(sim.dat$index[i,])] <- exp(
              log(pred.index[i,!is.na(resid.index[i,])])
              + as.numeric(sample(resid.index[i,!is.na(resid.index[i,])], replace=TRUE))
            )
          }
          for (i in 1:nrow(resid.caa)) {
            present <- which(!is.na(resid.caa[i, ]))
            if (length(present)) {
              residuals <- resid.caa[i,present]
              sim.dat$caa[i,present] <- exp(log(pred.caa[i,present]) +
                residuals[sample.int(length(residuals), replace=TRUE)])
            }
          }
        }

        sim.res <- Res0
        sim.res$input$dat <- sim.dat
        sim.res$input$silent = TRUE
        sim.res$input$p0.list <- Res0$par_list
        sim.res <- try(do.call(sam, sim.res$input))
        if (class(sim.res) != "try-error") break
      }
      cat(sprintf("-----%s-----\n",j))
      boot.list[[j]] <- sim.res
    }
    boot.list
  } else {
    # if (!isTRUE(Res$input$get.random.vcov)) {
    #   Res$input$get.random.vcov <- TRUE
    #   Res <- do.call(sam, Res$input)
    # }
    if(is.null(Res$rep$unbiased)) {
      mu0 <- Res$rep$value
    } else {mu0 <- Res$rep$unbiased$value}
    mu <- mu0[1:(length(mu0)/2)]
    vcov <- Res$rep$cov[1:(length(mu0)/2),1:(length(mu0)/2)]
    nc <- ncol(Res$naa)
    cn <- colnames(Res$naa)
    rn <- rownames(Res$naa)
    waa <- Res$input$dat$waa
    maa <- Res$input$dat$maa
    lapply(1:n, function(i){
      boot <- MASS::mvrnorm(n=1,mu,vcov)
      naa <- exp(matrix(boot[names(boot)=="logN"],ncol=nc))
      faa <- exp(matrix(boot[names(boot)=="logF"],ncol=nc))
      faa <- rbind(faa,Res$input$alpha*faa[nrow(faa),])
      colnames(naa) <- colnames(faa) <- cn
      rownames(naa) <- rownames(faa) <- rn
      baa <- naa*waa
      ssb <- baa*maa
      list(naa=naa,baa=baa,ssb=ssb,faa=faa)
    })
  }
}

# Catch observation groups and their log-scale errors must match the SAM
# likelihood in both bootstrap and population-simulation data.
sam_catch_simulation_inputs <- function(Res) {
  pred.caa <- as.matrix(Res$caa)
  sigma.caa <- matrix(Res$sigma.logC, nrow(pred.caa), ncol(pred.caa))
  if (isTRUE(Res$input$change_plusgroup)) {
    catch.obs <- Res$data$obs[Res$data$obs[, "fleet"] == 1, , drop = FALSE]
    for (i in seq_len(nrow(catch.obs))) {
      ages <- seq.int(catch.obs[i, "age"], catch.obs[i, "maxage"]) + 1L
      year <- match(as.character(catch.obs[i, "year"]), colnames(pred.caa))
      if (length(ages) > 1L) {
        catches <- Res$caa[ages, year]
        pred.caa[ages[1L], year] <- sum(catches)
        sigma.caa[ages[1L], year] <-
          sqrt(sum((catches * sigma.caa[ages, year])^2)) / sum(catches)
      }
    }
  }
  list(mean = pred.caa, sd = sigma.caa,
       present = !is.na(Res$input$dat$caa))
}
