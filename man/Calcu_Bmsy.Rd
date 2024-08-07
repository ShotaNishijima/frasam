% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calcu_BRs.R
\name{Calcu_Bmsy}
\alias{Calcu_Bmsy}
\title{Calculate the stock biomass at the maximum sustainable yield (Bmsy).}
\usage{
Calcu_Bmsy(Fmsy, M, Sel, w, g, A = 6, alpha, beta, method_SR)
}
\arguments{
\item{M}{A vector of natural mortality rate at age a. The length should be A+1.}

\item{Sel}{A vector of selectivity at age a. The length should be A+1.}

\item{w}{A vector of weight at age a. The length should be A+1.}

\item{g}{A vector of maturity at age a. The length should be A+1.}

\item{A}{Plus group age. Default is 6.}

\item{alpha}{Estimated Parameter alpha in the Beverton-Holt stock recruitment relationship.}

\item{beta}{Estimated Parameter beta in the Beverton-Holt stock recruitment relationship.}
}
\value{
A value of the stock biomass at the maximum sustainable yield (Bmsy).
}
\description{
\emph{Calcu_Bmsy} calculates the stock biomass at the maximum sustainable yield (Fmsy) for the Bevertoh-Holt stock-recruitment relationship.
}
