\name{ca_verify}
\alias{ca_verify}
\alias{print.ca_verify}
\title{Verify strength-t coverage of a covering array}
\usage{
ca_verify(x, t, v = NULL, threads = NULL, report = 10L)
}
\arguments{
  \item{x}{an integer matrix or data frame of integers; rows are runs,
    columns are factors. \code{NA} entries are wildcards ("flexible
    values") and count as every symbol.}
  \item{t}{interaction strength to verify (positive integer).}
  \item{v}{number of symbols per column; if \code{NULL}, inferred from
    the data range. Symbols may be coded \code{0..v-1} or \code{1..v};
    1-based input is detected and shifted automatically.}
  \item{threads}{number of threads. \code{NULL} (the default) picks
    automatically: half the machine's logical cores, single-threaded
    for small jobs, capped during CRAN checks. Override with
    \code{options(caverify.threads = n)} or this argument. Long runs
    are interruptible without aborting R.}
  \item{report}{maximum number of missing-tuple examples to collect.}
}
\value{
  An object of class \code{ca_verify}: a list with \code{covered}
  (logical), \code{colsets}, \code{gaps}, \code{missing_tuples},
  \code{examples} (matrix; first \code{t} entries per row are 1-based
  column indices, the rest the missing value combination in the input
  coding), and \code{t}, \code{v}, \code{N}, \code{k}.
}
\description{
  Checks that an \eqn{N \times k} array over \eqn{v} symbols covers
  every \eqn{t}-way interaction: for every choice of \eqn{t} columns
  and every one of the \eqn{v^t} value combinations, at least one row
  exhibits that combination. This is the certificate check for a
  covering array, run in compiled C at speeds suitable for large
  arrays and strengths.
}
\examples{
ca <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
ca_verify(ca, t = 2)          # a covering array of strength 2
ca_verify(ca[-1, ], t = 2)    # removing a run breaks coverage
}
