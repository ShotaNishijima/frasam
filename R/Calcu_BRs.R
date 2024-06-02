
#' Calculate the relative number at age at the equilibrium.
#'
#' \emph{Calcu_N} calculates the relative number at age at the equilibrium (N) with a given fishing mortality value.
#' 
#' @param Fish_mort Fishing mortality to calculate SBR. Should be provided with a numeric value.
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @export

Calcu_N <- function(Fish_mort, M, Sel, A = 6){
  ncura <- rep(as.numeric(NA), A+1)
  ncura[1] <- 1 # a = 0
  for(a in 1:(A-1)){ # 1 =< a =< A-1
    ncura[a+1] <- exp(-sum((M + Sel*Fish_mort)[1:a]))
  }
  ncura[A+1] <- exp(-sum((M + Sel*Fish_mort)[1:A])) / (1 - exp(- M[A+1] - Sel[A+1]*Fish_mort))
  return(ncura)
}


#' Calculate spawning biomass per recruit (SBR) with a given fishing mortality value.
#'
#' \emph{Calcu_SBR} calculates SBR from the statistics derived from the \emph{Calcu_stat} function and input data with a given fishing mortality value.
# #' @usage Calcu_SBR(Fish_mort = 0.2, M = M, Sel = Sel, w = w, g = g, A = 6)
#' @param Fish_mort Fishing mortality to calculate SBR. Should be provided with a numeric value.
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @return A value of spawing biomass par recruit with the given fishing mortality.
#' @export
#' 
#'
 

Calcu_SBR <- function(Fish_mort, M, Sel, w, g, A = 6){
  ncura <- Calcu_N(Fish_mort = Fish_mort, M = M, Sel = Sel, A = A)
  sbr_f <- sum(w * g * ncura)
  return(sbr_f)
}

#' Calculate F\%SPR: the fishing mortality that reduces the spawning biomass per recruit (SBR) to X \% of those without fishing (SPR0).
#'
#' \emph{Calcu_SPR_X} calculates F\%SPR (the fishing mortality that reduce the spawning biomass per recruit (SBR) to X \% of those without fishing (SPR0)).
# #' @usage Calcu_SPR_X(x = c(0.2, 0.3, 0.4), M = M, Sel = Sel, w = w, g = g, A = 6, F_max = 1)
#' @param x The proportion (X) to calculate F\% SPR. This parameter should be provided as a proportion (i.e. within the range from 0 to 1). Defalut is 0.2.
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param F_max 
#' This value should be large enough so that the true F\%SPR value be included (otherwise a warning message appears). Defalut is 1.
#' @return A value of F\% SPR.
#' 
#' @export

Calcu_SPR_X <- function(x, M, Sel, w, g, A = 6, F_max = 2){
  SBR0 <- Calcu_SBR(Fish_mort = 0, M = M, Sel = Sel, w = w, g = g, A = A)
  obj_fun <- function(Fish_mort) {
    SBRF <- Calcu_SBR(Fish_mort = Fish_mort, M = M, Sel = Sel, w = w, g = g, A = A)
    return((SBRF/SBR0 - x)^2)
  }
  opt <- optimize(obj_fun, interval = c(0, F_max))
  return(opt$minimum)
}

#' Calculate yield per recruit with a given fishing mortality value.
#'
#' \emph{Calcu_YPR} calculates yield per recruit from the statistics derived from the \emph{Calcu_stat} function and input data with a given fishing mortality value.
# #' @usage Calcu_YPR(Fish_mort = 0.2, M = M, Sel = Sel, w = w, g = g, A = 6, method = "Baranov)
#' @param Fish_mort Fishing mortality to calculate yield per recruit. Should be provided with a numeric value.
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param method Method of the fishing equation. Either "Baranov" or "Pope" is allowed. Default is "Baranov".
#' @return A value of yield per recruit with the given fishing mortality.
#' @export

Calcu_YPR <- function(Fish_mort, M, Sel, w, g, A = 6, method = "Baranov"){
  ncura <- Calcu_N(Fish_mort = Fish_mort, M = M, Sel = Sel, A = A)
  if(method == "Baranov"){
    YPRF <- sum(Sel*Fish_mort/(M+Sel*Fish_mort) * (1 - exp(-M-Sel*Fish_mort)) * w * ncura)
  } else if(method == "Pope"){
    YPRF <- sum((1-exp(-Sel*Fish_mort)) * w * ncura * exp(-M/2))
  } else {
    stop("Unexpected method is provided. method should be either Baranov or Pope.")
  }
  return(YPRF)
}

#'
#' Calculate F0.1 value.
#'
#' \emph{Calcu_F0.1} calculates the fishing mortality rate at which the slope of the yield-per-recruit curve is 10\% of the slope of the curve at its origin.
#' 
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param method Method of the fishing equation. Either "Baranov" or "Pope" is allowed. Default is "Baranov".
#' @param F_max  Maximum value of F values to be searched for determining F\%SPR. This value should be large enough so that the true F\%SPR value be included. Default is 10
# #' @usage Calcu_F0.1(M = M, Sel = Sel, w = w, g = g, A = 6, method = "Baranov, F_max = 10)
#' 
#' @return A value of the F0.1.
#' 
#' @export

Calcu_F0.1 <- function(M, Sel, w, g, A = 6, method = "Baranov", F_max = 10){

  Marginal_YPR0 <- (Calcu_YPR(Fish_mort = 0.001, M = M, Sel = Sel, w = w, g = g, A = A, method = method) -
                      Calcu_YPR(Fish_mort = 0, M = M, Sel = Sel, w = w, g = g, A = A, method = method)) / 0.001

  obj_fun <- function(x){
    Marginal_YPRF <- (Calcu_YPR(Fish_mort = x+0.001, M = M, Sel = Sel, w = w, g = g, A = A, method = method) -
                        Calcu_YPR(Fish_mort = x, M = M, Sel = Sel, w = w, g = g, A = A, method = method)) / 0.001
    return((Marginal_YPRF/Marginal_YPR0 - 0.1)^2)
  }
  opt <- optimize(obj_fun, interval = c(0, F_max))
  return(opt$minimum)
}

#' Calculate the fishing mortality that maximizes the yield per recruit.
#'
#' \emph{Calcu_Fmax} calculates the fishing mortality that maximizes the yield per recruit.
# #' @usage Calcu_SPR_X(M = M, Sel = Sel, w = w, g = g, A = 6, F_max = 3)
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param method Method of the fishing equation. Either "Baranov" or "Pope" is allowed. Default is "Baranov".
#' @param F_max  Maximum value of F values to be searched for determining F\%SPR.
#' This value should be large enough so that the true F\%SPR value be included. Default is 10
#' @return A value of the F that maximizes the yield per recruit.
#' @export

Calcu_Fmax <- function(M, Sel, w, g, A = 6, F_max = 10, method = "Baranov"){
  obj_fun <- function(x){
    YPRF <- Calcu_YPR(Fish_mort = x, M = M, Sel = Sel, w = w, g = g, A = A, method = method)
    return(-YPRF)
  }
  opt <- optimize(obj_fun, interval = c(0, F_max))
  return(opt$minimum)
}

#' Estimate the Beverton-Holt and hockey-stick stock recruitment relationships from data of spawning stock biomass (SBy) and recruits (Ry) with a give steepness (h) value.
#'
#' \emph{Est_SR} estimates the coefficients (i.e. alpha and beta) of the BH relationship from data of spawning stock biomass (SBy) and recruits (Ry).
# #' @usage Est_SR(SBy = SBy, Ry = Ry, h = 0.2, h_fix = T, ref.year = c(1970:2019), M = M, Sel = Sel, w = w, g = g, method = "BH")
#' @param SBy A vector of spawning biomass at year. The vector should be named with years. The length should be the same as that of \emph{ref.year}
#' @param Ry A vector of recruits at year. The vector should be named with years. The length should be the same as that of \emph{ref.year}
#' @param h A value of steepness to be fixed.
#' @param h_fix A logical value deciding whether the steepness is fixed. Default is TRUE. FALSE is allowed only for the BH relationship at this stage.
#' @param inits A vector of initial values of the parameters search (alpha, beta) for BH. Although default is c(1, 1), searching for the suitable sets is recommended.
#' @param ref.year A vector of years to be referred for the estimation. Default is c(1970:2019)
#' @param SB0_max  Maximum value of SB0 values to be searched for parameter estimation. Default is max(SBy)*100.
#' @param Fish_mort Fishing mortality to calculate SBR. Should be provided with a numeric value.
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
# #' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param method The type of the stock-recruitment relationship. At this stage, "BH" (Beverton-Holt) or "HS" (Hockey-stick) is allowed.
#' @return Values of alpha and beta coefficients for the BH relationships or a value of SB0 for the HS relationships. 
#' As for the BH, when h is not fixed (i.e. h_fix = T), the output from the optim() function is shown.
#' @export

Est_SR <- function(SBy, Ry, h, h_fix = T, inits = c(1,1), ref.year = c(1970:2019), SB0_max = max(SBy)*100,
                   M, w, g, A = 6, method){
  
  SB_obs <- SBy[names(SBy) %in% ref.year]
  R_obs <- Ry[names(Ry) %in% ref.year]
  SBR0 <- Calcu_SBR(Fish_mort = 0, M = M, Sel = rep(1,A+1), w = w, g = g, A = A)
  
  if(method == "BH"){
    if(isTRUE(h_fix)){
      alpha <- 4*h / (SBR0*(1-h))
      obj_fun <- function(SB0){ # Estimate SB0.
        beta <- (5*h-1) / ((1-h)*SB0)
        R_fit <- alpha*SB_obs / (1 + beta*SB_obs)
        return(sum((log(R_fit) - log(R_obs))^2))
      }
      opt <- optimize(obj_fun, interval = c(0, SB0_max))
      output <- c(alpha, (5*h-1) / ((1-h)*opt$minimum))
      names(output) <- c("alpha", "beta")
      SB0 <- opt$minimum
      output <- c(output,sigma=sqrt(opt$objective/length(SB_obs)),SBR0=SBR0,SB0=SB0,R0=SB0/SBR0,h=h)
      return(output)
    } else if(!isTRUE(h_fix)) {
      obj_fun <- function(par){ # Estimate alpha and beta
        R_fit <- par[1]*SB_obs / (1 + par[2]*SB_obs)
        return(sum((log(R_fit) - log(R_obs))^2))
      }
      opt <- optim(par = inits, fn = obj_fun)
      if(opt$convergence != 0) warning("Stock-recruitment parameter estimation may not converge")
      return(opt$par)
    }
    
  } else if(method == "HS"){
    if(isTRUE(h_fix)){
      stop("Fixing F is not implemented for the Hockey-stick stock recruitment relationship.")
      # obj_fun <- function(SB0){ # Estimate SB0.
      #   R0 <- SB0 / SBR0
      #   R_fit <- R0*SB_obs / ((1-h)*SB0)
      #   R_fit[SB_obs >= (1-h)*SB0] <- R0
      #   return(sum((log(R_fit) - log(R_obs))^2))
      # }
      # opt <- optimize(obj_fun, interval = c(min(SB_obs)/(1-h), max(SB_obs)/(1-h)))
      # beta <- (1-h) * opt$minimum
      # alpha <- opt$minimum / SBR0
      # return(c(alpha, beta))
    } else if(!isTRUE(h_fix)) {
      # obj_fun <- function(par){ # Estimate alpha and beta
      #   R_fit <- par[1] * SB_obs
      #   R_fit[SB_obs <= par[2]] <- par[1] * par[2]
      #   return(sum((log(R_fit) - log(R_obs))^2))
      # }
      obj_fun <- function(beta,out=FALSE){ # Estimate beta
        SB_b <- SB_obs
        SB_b[SB_obs>beta] <- beta
        alpha <- exp(mean(log(R_obs/SB_b)))
        R_fit <- alpha * SB_b
        if (isTRUE(out)) {
          return( c(alpha=alpha,beta=beta) )
        } else {
          return(sum((log(R_fit) - log(R_obs))^2))
        }
      }
      opt <- optimize(obj_fun,interval=range(SB_obs))
      RES <- c(obj_fun(opt$minimum,out=TRUE),sigma=sqrt(opt$objective/length(SB_obs)))
      R0 <- prod(RES[1:2])
      SB0 <- R0*SBR0
      h <- 1-opt$minimum/SB0
      RES <- c(RES,SBR0=SBR0,SB0=SB0,R0=R0,h=h)
      return( RES )
      # opt <- optim(par = inits, fn = obj_fun)
      # if(opt$convergence != 0) warning("Stock-recruitment parameter estimation may not converge")
      # return(c(alpha = opt$par[1], beta = opt$par[2]))
      # 
      # obj_fun <- function(par){ # Estimate h and SB0.
      #   R0 <- par[2] / SBR0
      #   R_fit <- R0*SB_obs / ((1-par[1])*par[2])
      #   R_fit[SB_obs >= (1-par[1])*par[2]] <- R0
      #   return(sum((log(R_fit) - log(R_obs))^2))
      # }
      # opt <- optim(par = inits, fn = obj_fun, method = "L-BFGS-B",
      #              lower = c(0, 0), upper = c(1, max(SB_obs)))
      # if(opt$convergence != 0) warning("Stock-recruitment parameter estimation may not converge")
      # beta <- (1-opt$par[1]) * opt$par[2]
      # alpha <- opt$par[2] / SBR0
      # return(c(alpha, beta))
    }
    
  } else {
    stop("Unexpected method. Confirm that method is either BH or HS")
  }
}

#' Calculate the fishing mortality that maximizes the maximum sustainable yield (Fmsy) with fixed steepness parameters (h) in the Bevertoh-Holt stock-recruitment relationship.
#'
#' \emph{Calcu_Fmsy} calculates the fishing mortality that maximizes the maximum sustainable yield (Fmsy) with fixed steepness parameters (h).
# #' @usage 
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param method Method of the fishing equation. Either "Baranov" or "Pope" is allowed.
#' @param F_max  Maximum value of F values to be searched for determining Fmsy.
#' @param alpha Estimated Parameter alpha in the Beverton-Holt stock recruitment relationship.
#' @param beta Estimated Parameter beta in the Beverton-Holt stock recruitment relationship.
#' @return A value of the F that maximizes the yield per recruit.
#' @export

Calcu_Fmsy <- function(M, Sel, w, g, A = 6, method,
                          F_max = 10, alpha, beta, method_SR){
  if(method_SR == "BH"){
    obj_fun <- function(Fish_mort){
      YPRF <- Calcu_YPR(Fish_mort = Fish_mort, M = M, Sel = Sel, w = w, g = g, A = A, method = method)
      SBRF <- Calcu_SBR(Fish_mort = Fish_mort, M = M, Sel = Sel, w = w, g = g, A = A)
      SYF <- YPRF * (alpha*SBRF - 1) / (beta*SBRF)
      return(-SYF)
    }
    opt <- optimize(obj_fun, interval = c(0, F_max))
    return(opt$minimum)
    
  } else if(method_SR == "HS"){
    obj_fun_Fstar <- function(Fish_mort){
      SBRF <- Calcu_SBR(Fish_mort = Fish_mort, M = M, Sel = Sel, w = w, g = g, A = A)
      return((alpha*SBRF - 1)^2)
    }
    opt_Fstar <- optimize(obj_fun_Fstar, interval = c(0, F_max))
    F_star <- opt_Fstar$minimum
    obj_fun <- function(Fish_mort){
      YPRF <- Calcu_YPR(Fish_mort = Fish_mort, M = M, Sel = Sel, w = w, g = g, A = A, method = method)
      # if(Fish_mort <= F_star){
        SYF <- YPRF*alpha*beta
      # } else {
      #   SYF <- 0
      # }
      return(-SYF)
    }
    opt <- optimize(obj_fun, interval = c(0, F_star))
    # opt <- optimize(obj_fun, interval = c(0, F_max))
    return(opt$minimum)
    
  } else {
    stop("Unexpected method. Confirm that method is either BH or HS")
  }
}

#' Calculate the stock biomass at the maximum sustainable yield (Bmsy).
#'
#' \emph{Calcu_Bmsy} calculates the stock biomass at the maximum sustainable yield (Fmsy) for the Bevertoh-Holt stock-recruitment relationship.
# #' @usage 
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param alpha Estimated Parameter alpha in the Beverton-Holt stock recruitment relationship.
#' @param beta Estimated Parameter beta in the Beverton-Holt stock recruitment relationship.
#' @return A value of the stock biomass at the maximum sustainable yield (Bmsy).
#' 
#' @export
#' 

Calcu_Bmsy <- function(Fmsy, M, Sel, w, g, A = 6,
                       alpha, beta, method_SR){
  ncura_msy <- Calcu_N(Fish_mort = Fmsy, M = M, Sel = Sel, A = A)
  SBR_msy <- Calcu_SBR(Fish_mort = Fmsy, M = M, Sel = Sel, w = w, g = g, A = A)
  if(method_SR == "BH"){
    R <- (alpha*SBR_msy - 1) / (beta*SBR_msy)
    Bmsy <- sum(w * R * ncura_msy)
  } else if(method_SR == "HS"){
    Bmsy <- sum(w*alpha*beta*ncura_msy)
  } else {
    stop("Unexpected method. Confirm that method is either BH or HS")
  }
  return(Bmsy)
}

#' Calculate the spawning stock biomass at the maximum sustainable yield (SSBmsy).
#'
#' \emph{Calcu_SBmsy} calculates the spawning stock biomass at the maximum sustainable yield (Fmsy) for the Bevertoh-Holt stock-recruitment relationship.
# #' @usage 
#' @param M A vector of natural mortality rate at age a. The length should be A+1.
#' @param Sel A vector of selectivity at age a. The length should be A+1.
#' @param w A vector of weight at age a. The length should be A+1.
#' @param g A vector of maturity at age a. The length should be A+1.
#' @param A Plus group age. Default is 6.
#' @param alpha Estimated Parameter alpha in the Beverton-Holt stock recruitment relationship.
#' @param beta Estimated Parameter beta in the Beverton-Holt stock recruitment relationship.
#' @return A value of the spawning stock biomass at the maximum sustainable yield (SSBmsy).

Calcu_SBmsy <- function(Fmsy, M, Sel, w, g, A = 6,
                        alpha, beta, method_SR){
  ncura_msy <- Calcu_N(Fish_mort = Fmsy, M = M, Sel = Sel, A = A)
  SBR_msy <- Calcu_SBR(Fish_mort = Fmsy, M = M, Sel = Sel, w = w, g = g, A = A)
  if(method_SR == "BH"){
    SBmsy <- (alpha*SBR_msy - 1) / beta
  } else if(method_SR == "HS"){
    SBmsy <- alpha*beta*SBR_msy
  } else {
    stop("Unexpected method. Confirm that method is either BH or HS")
  }
  return(SBmsy)
}

