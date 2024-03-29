% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poho_residual.R
\name{ph_residual}
\alias{ph_residual}
\title{Post-Hoc Residual Test}
\usage{
ph_residual(
  field1,
  field2,
  categories1 = NULL,
  categories2 = NULL,
  residual = "adjusted"
)
}
\arguments{
\item{field1}{list or dataframe with the first categorical field}

\item{field2}{list or dataframe with the second categorical field}

\item{categories1}{optional list with order and/or selection for categories of field1}

\item{categories2}{optional list with order and/or selection for categories of field2}

\item{residual}{optional methdod for residual to test. Either "adjusted" (default) or "standardized".}
}
\description{
Post-Hoc Residual Test
}
