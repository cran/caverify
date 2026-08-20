#' Verify strength-t coverage of a covering array
#'
#' Checks that an N x k array over v symbols covers every t-way
#' interaction: for every choice of t columns and every one of the v^t
#' value combinations, at least one row exhibits that combination. This
#' is the certificate check for a covering array, performed in C at
#' speeds suitable for large arrays and strengths.
#'
#' Missing values (NA) are treated as wildcards ("flexible values"):
#' an NA entry counts as every symbol simultaneously.
#'
#' Symbols may be coded 0..v-1 or 1..v; 1-based input is detected and
#' shifted automatically (reported examples are given in the input
#' coding).
#'
#' @param x an integer matrix or data frame of integers, rows = runs,
#'   columns = factors.
#' @param t interaction strength to verify (a positive integer).
#' @param v number of symbols per column. If NULL (default), inferred
#'   from the data range.
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
#'   parameters \code{t}, \code{v}, \code{N}, \code{k}.
#' @examples
#' ## a strength-2 covering array on 3 two-level factors in 4 runs
#' ca <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
#' ca_verify(ca, t = 2)
#'
#' ## removing a run breaks coverage
#' ca_verify(ca[-1, ], t = 2)
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

    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1L])) stop("'x' contains no non-missing values")
    shift <- 0L
    if (is.null(v)) {
        if (rng[1L] == 0L) v <- rng[2L] + 1L
        else if (rng[1L] == 1L) { v <- rng[2L]; shift <- 1L }
        else stop("cannot infer 'v': symbols start at ", rng[1L],
                  " (expected 0- or 1-based); supply 'v' and 0-based symbols")
    } else {
        v <- as.integer(v)
        if (rng[1L] >= 1L && rng[2L] == v) shift <- 1L
        else if (rng[1L] >= 0L && rng[2L] <= v - 1L) shift <- 0L
        else stop("symbols out of range for v = ", v,
                  " (saw ", rng[1L], "..", rng[2L], ")")
    }
    if (shift) x <- x - 1L

    nt <- resolve_threads(threads, ncol(x), t)
    res <- .Call(C_ca_verify, x, t, v, nt, as.integer(report))
    if (shift && nrow(res$examples) > 0L) {
        tt <- ncol(res$examples) / 2L
        res$examples[, (tt + 1L):(2L * tt)] <-
            res$examples[, (tt + 1L):(2L * tt), drop = FALSE] + 1L
    }
    res$t <- t; res$v <- v; res$N <- nrow(x); res$k <- ncol(x)
    class(res) <- "ca_verify"
    res
}

#' @export
print.ca_verify <- function(x, ...) {
    if (isTRUE(x$covered)) {
        cat(sprintf(
            "VERIFIED: CA(N=%d; t=%d, k=%d, v=%d) - all %s column sets cover all %s tuples\n",
            x$N, x$t, x$k, x$v,
            format(x$colsets, big.mark = ","),
            format(x$v^x$t, big.mark = ",")))
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
