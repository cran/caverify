\name{ca_verify}
\alias{ca_verify}
\alias{print.ca_verify}
\title{Verify strength-t coverage of a covering array}
\usage{
ca_verify(x, t, v = NULL, threads = NULL, report = 10L)
}
\arguments{
  \item{x}{an integer matrix or data frame of integers; rows are runs,
    columns are factors. \code{NA} marks a flexible ("don't care")
    entry: a row contributes nothing to a projection in which it has
    an \code{NA}, so a verified array is covering no matter how its
    \code{NA} entries are later filled.}
  \item{t}{interaction strength to verify (positive integer).}
  \item{v}{number of symbols. One of: \code{NULL} (the default) to
    infer per-column symbol counts from each column's own maximum
    (\code{"auto"} is an explicit alias); a single integer
    \code{>= 2} declaring a uniform array; or an integer vector of
    length \code{ncol(x)} giving each column its own number of
    symbols. Versions 0.1.x inferred one uniform value from the
    global data range instead; changed on a recommendation by Ulrike
    Groemping, a default should not assume a uniform CA. Declare
    \code{v} explicitly when the check should also catch a column
    that fails to reach its intended number of symbols. Symbols may
    be coded \code{0..v-1} or \code{1..v} (per column for mixed
    levels); 1-based input is detected and shifted automatically.}
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
  coding), and \code{t}, \code{v} (a single integer for uniform
  arrays, an integer vector of per-column symbol counts for
  mixed-level arrays), \code{N}, \code{k}.
}
\description{
  Checks that an \eqn{N \times k} array covers every \eqn{t}-way
  interaction: for every choice of \eqn{t} columns and every
  combination of their symbols, at least one row exhibits that
  combination. This is the certificate check for a covering array, run
  in compiled C at speeds suitable for large arrays and strengths.

  Both uniform arrays (every column has \eqn{v} symbols, \eqn{v^t}
  combinations per column set) and mixed-level arrays (column \eqn{i}
  has its own \eqn{v_i} symbols, \eqn{\prod v_i} combinations per
  column set) are supported; a uniform array is the degenerate
  mixed-radix case in which all radices coincide.
}
\examples{
ca <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
ca_verify(ca, t = 2)          # a covering array of strength 2
ca_verify(ca[-1, ], t = 2)    # removing a run breaks coverage

# a mixed-level array: one 3-level and two 2-level columns
mca <- as.matrix(expand.grid(0:2, 0:1, 0:1))
ca_verify(mca, t = 2, v = c(3, 2, 2))
ca_verify(mca, t = 2)   # the default infers each column's count
}
