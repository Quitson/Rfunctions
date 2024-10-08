% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{bs_CI}
\alias{bs_CI}
\title{bs_CI}
\usage{
bs_CI(data, conf_level = 0.95)
}
\arguments{
\item{data}{vector}

\item{conf_level}{confidence level}
}
\value{
confidence interval (vector)
}
\description{
Bootstrap is used e.g. for very skewed distributions when CLT is inaccurate. Bootstrap is not advised for large data sets, almost regardless of distribution
since the CIs are virtually the same (as CLT CI) in this case and bootstrap is calculation intense. Bootstrap t CI (number of resamples = 10^4)
}
