% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.R
\name{contAdjustBin}
\alias{contAdjustBin}
\title{Continuity adjusted binomial limits}
\usage{
contAdjustBin(limit, n, target)
}
\arguments{
\item{limit}{a numeric vector containing the (1-limit)100\% values for the control limits.}

\item{n}{number of clusters (institutions) being compared.}

\item{target}{institution target value (a proportion)}
}
\description{
Adjust binomial distribution control limits for lack of continuity. See #4 of A.1.1 of Spiegelhalter (2005).
}
