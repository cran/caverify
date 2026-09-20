#' Verify strength-t coverage of a covering array
#'
#' Checks that an N x k array covers every t-way interaction: for every
#' choice of t columns and every combination of their symbols, at least
#' one row exhibits that combination. This is the certificate check for
#' a covering array, performed in C at speeds suitable for large arrays
#' and strengths.
#'
#' Both uniform arrays (every column has v symbols) and mixed-level
#' arrays (column i has its own number of symbols v[i]) are supported.
#' Tuple counts and indexing are mixed-radix; a uniform array is the
#' degenerate case in which all radices coincide.
#'
#' Missing values (NA) mark flexible ("don't care") entries. A row
#' contributes nothing to a projection in which it has an NA, so a
#' verified array is covering no matter how its NA entries are later
#' filled. (Versions 0.1.x instead let an NA count as every symbol at
#' once, which could certify arrays that no single choice of values
#' would make covering; changed in 0.2.0 following a counterexample
#' by Ulrike Groemping.)
#'
#' Symbols may be coded 0..v-1 or 1..v; 1-based input is detected and
#' shifted automatically (reported examples are given in the input
#' coding).
#'
#' @param x an integer matrix or data frame of integers, rows = runs,
#'   columns = factors.
#' @param t interaction strength to verify (a positive integer).
#' @param v number of symbols. One of: NULL (default) to infer
#'   per-column symbol counts from each column's own maximum ("auto"
#'   is an explicit alias for the same behaviour); a single integer
#'   >= 2 declaring a uniform array; or an integer vector of length
#'   ncol(x) giving each column its own number of symbols. Versions
#'   0.1.x instead inferred one uniform value from the global data
#'   range; changed on a recommendation by Ulrike Groemping, since a
#'   default should not assume a uniform CA. Declare v explicitly
#'   when the check should also catch a column that fails to reach
#'   its intended number of symbols.
#' @param threads number of threads. NULL (default) picks automatically:
#'   half the machine's logical cores (data.table-style politeness),
#'   single-threaded for small jobs, capped during CRAN checks, and
#'   overridable via options(caverify.threads = n) or this argument.
#'   Effective when compiled with OpenMP; otherwise single-threaded.
#'   Long runs are interruptible (Escape or Ctrl-C).
#' @param report maximum number of missing-tuple examples to collect.
#' @return An object of class \code{ca_verify}: a list with elements
#'   \code{covered} (logical), \code{colsets} (number of column sets
#'   checked), \code{gaps} (column sets with at least one missing
#'   tuple), \code{missing_tuples} (total missing tuples),
#'   \code{examples} (matrix of up to \code{report} missing examples:
#'   first t columns give 1-based column indices, remaining t columns
#'   the missing value combination in the input coding), and the call
#'   parameters \code{t}, \code{v} (a single integer for uniform
#'   arrays, an integer vector of per-column symbol counts for
#'   mixed-level arrays), \code{N}, \code{k}.
#' @examples
#' ## a strength-2 covering array on 3 two-level factors in 4 runs
#' ca <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
#' ca_verify(ca, t = 2)
#'
#' ## removing a run breaks coverage
#' ca_verify(ca[-1, ], t = 2)
#'
#' ## a mixed-level array: one 3-level and two 2-level columns,
#' ## full factorial, so any strength-2 projection is covered
#' mca <- as.matrix(expand.grid(0:2, 0:1, 0:1))
#' ca_verify(mca, t = 2, v = c(3, 2, 2))
#' ca_verify(mca, t = 2)   # the default infers each column's count
#' @export
resolve_threads <- function(threads, k, t) {
    if (!is.null(threads)) return(max(1L, as.integer(threads)))
    opt <- getOption("caverify.threads", NULL)
    if (!is.null(opt)) return(max(1L, as.integer(opt)))
    nc <- tryCatch(parallel::detectCores(logical = TRUE),
                   error = function(e) 1L)
    if (is.na(nc) || nc < 1L) nc <- 1L
    auto <- max(1L, nc %/% 2L)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && !identical(tolower(chk), "false"))
        auto <- min(auto, 2L)
    if (choose(k, t) < 4096) auto <- 1L
    auto
}

ca_verify <- function(x, t, v = NULL, threads = NULL, report = 10L) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x)) stop("'x' must be a matrix or data frame")
    storage.mode(x) <- "integer"
    t <- as.integer(t)
    if (length(t) != 1L || is.na(t) || t < 1L) stop("'t' must be a positive integer")
    k <- ncol(x)

    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1L])) stop("'x' contains no non-missing values")
    shift <- 0L

    if (is.null(v) || identical(v, "auto")) {
        ## default: per-column inference from each column's own maximum.
        ## NULL behaves this way on U. Groemping's recommendation, since
        ## a default should not assume a uniform CA ("auto" is kept as
        ## an explicit alias). Note the tradeoff: an intended-uniform
        ## array with a defective column (top symbol absent, or
        ## constant) passes at the reduced per-column counts; declare a
        ## scalar v where that defect-catching matters.
        colmax <- suppressWarnings(apply(x, 2L, max, na.rm = TRUE))
        if (any(!is.finite(colmax)))
            stop("cannot infer per-column symbol counts: column(s) ",
                 paste(which(!is.finite(colmax)), collapse = ", "),
                 " are all NA; supply 'v'")
        if (rng[1L] == 0L) vs <- as.integer(colmax) + 1L
        else if (rng[1L] == 1L) { vs <- as.integer(colmax); shift <- 1L }
        else stop("cannot infer symbol counts: symbols start at ", rng[1L],
                  " (expected 0- or 1-based); supply 'v' and 0-based symbols")
    } else {
        if (!is.numeric(v) || anyNA(v) || any(v != as.integer(v)))
            stop("'v' must be NULL, \"auto\", a single integer, or an integer vector of length ncol(x)")
        v <- as.integer(v)
        if (length(v) == 1L) {
            if (v < 2L) stop("'v' must be >= 2")
            if (rng[1L] >= 1L && rng[2L] == v) shift <- 1L
            else if (rng[1L] >= 0L && rng[2L] <= v - 1L) shift <- 0L
            else stop("symbols out of range for v = ", v,
                      " (saw ", rng[1L], "..", rng[2L], ")")
            vs <- rep.int(v, k)
        } else {
            if (length(v) != k)
                stop("'v' has length ", length(v),
                     " but 'x' has ", k, " columns")
            if (any(v < 1L)) stop("all elements of 'v' must be >= 1")
            colmax <- suppressWarnings(apply(x, 2L, max, na.rm = TRUE))
            colmax[!is.finite(colmax)] <- -1L  ## all-NA column fits any coding
            if (rng[1L] >= 1L && all(colmax <= v) && any(colmax == v)) shift <- 1L
            else if (rng[1L] >= 0L && all(colmax <= v - 1L)) shift <- 0L
            else stop("symbols out of range for the given per-column 'v' ",
                      "(neither 0..v[i]-1 nor 1..v[i] fits every column)")
            vs <- v
        }
    }
    if (shift) x <- x - 1L

    nt <- resolve_threads(threads, k, t)
    res <- .Call(C_ca_verify, x, t, vs, nt, as.integer(report))
    if (shift && nrow(res$examples) > 0L) {
        tt <- ncol(res$examples) / 2L
        res$examples[, (tt + 1L):(2L * tt)] <-
            res$examples[, (tt + 1L):(2L * tt), drop = FALSE] + 1L
    }
    res$t <- t
    res$v <- if (length(unique(vs)) == 1L) vs[1L] else vs
    res$N <- nrow(x); res$k <- k
    class(res) <- "ca_verify"
    res
}

## "4^2 3 2^3"-style exponent notation for a per-column levels vector
levels_notation <- function(vs) {
    r <- rle(as.integer(vs))
    paste(ifelse(r$lengths > 1L,
                 paste0(r$values, "^", r$lengths),
                 as.character(r$values)),
          collapse = " ")
}

#' @export
print.ca_verify <- function(x, ...) {
    mixed <- length(x$v) > 1L
    if (isTRUE(x$covered)) {
        if (mixed) {
            cat(sprintf(
                "VERIFIED: MCA(N=%d; t=%d, k=%d, levels %s) - all %s column sets cover all their tuples\n",
                x$N, x$t, x$k, levels_notation(x$v),
                format(x$colsets, big.mark = ",")))
        } else {
            cat(sprintf(
                "VERIFIED: CA(N=%d; t=%d, k=%d, v=%d) - all %s column sets cover all %s tuples\n",
                x$N, x$t, x$k, x$v,
                format(x$colsets, big.mark = ","),
                format(x$v^x$t, big.mark = ",")))
        }
    } else {
        cat(sprintf(
            "NOT a covering array at t=%d: %s of %s column sets have gaps (%s missing tuples)\n",
            x$t, format(x$gaps, big.mark = ","),
            format(x$colsets, big.mark = ","),
            format(x$missing_tuples, big.mark = ",")))
        if (nrow(x$examples) > 0L) {
            tt <- ncol(x$examples) / 2L
            cat("first missing examples (columns -> value combination):\n")
            for (i in seq_len(nrow(x$examples))) {
                cat("  (", paste(x$examples[i, 1:tt], collapse = ","),
                    ") -> (", paste(x$examples[i, (tt + 1):(2 * tt)],
                                    collapse = ","), ")\n", sep = "")
            }
        }
    }
    invisible(x)
}
