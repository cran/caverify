library(caverify)

## pure-R brute-force oracle (shares no logic with the C code)
brute <- function(x, t, v) {
    k <- ncol(x); miss <- 0; gaps <- 0
    combs <- utils::combn(k, t)
    tuples <- as.matrix(expand.grid(rep(list(0:(v - 1)), t)))
    for (ci in seq_len(ncol(combs))) {
        cols <- combs[, ci]
        sub <- x[, cols, drop = FALSE]
        hit <- logical(nrow(tuples))
        for (r in seq_len(nrow(sub))) {
            row <- sub[r, ]
            if (anyNA(row)) next   ## don't-care: contributes nothing here
            ok <- rep(TRUE, nrow(tuples))
            for (j in seq_len(t)) ok <- ok & (tuples[, j] == row[j])
            hit <- hit | ok
        }
        if (!all(hit)) { gaps <- gaps + 1; miss <- miss + sum(!hit) }
    }
    list(covered = gaps == 0, gaps = gaps, missing = miss)
}

## 1. known CA(4; 2, 3, 2)
ca <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
r <- ca_verify(ca, 2)
stopifnot(isTRUE(r$covered), r$colsets == 3, r$gaps == 0)

## 2. full factorial 2^3, strengths 2 and 3
ff <- as.matrix(expand.grid(0:1, 0:1, 0:1))
stopifnot(isTRUE(ca_verify(ff, 2)$covered), isTRUE(ca_verify(ff, 3)$covered))

## 3. broken array detected, counts match brute force
br <- ca[-1, ]
r <- ca_verify(br, 2)
b <- brute(br, 2, 2)
stopifnot(!r$covered, r$gaps == b$gaps, r$missing_tuples == b$missing)
stopifnot(nrow(r$examples) >= 1)

## 4. don't-care semantics: an NA row contributes nothing, so adding
## one cannot restore coverage (0.2.0 semantics change)
wc <- rbind(br, c(NA, NA, NA))
rwc <- ca_verify(wc, 2)
rbr <- ca_verify(br, 2)
stopifnot(!rwc$covered, rwc$gaps == rbr$gaps,
          rwc$missing_tuples == rbr$missing_tuples)

## 4b. regression: the Groemping counterexample. One concrete column
## plus 20 all-NA columns. Under the per-column default the all-NA
## columns are uninferable and the call must error informatively;
## with v declared, every column set must report uncovered.
gx <- cbind(1:4, matrix(NA_integer_, 4, 20))
stopifnot(inherits(try(ca_verify(gx, 4), silent = TRUE), "try-error"))
rg <- ca_verify(gx, 4, v = rep(4L, 21))
stopifnot(!rg$covered, rg$gaps == rg$colsets, rg$colsets == choose(21, 4))

## 4c. the default is per-column inference (0.2.0, on U. Groemping's
## recommendation): a mixed-level array is judged correctly with no v
plan0 <- cbind(rep(0:2, each = 3), c(rep(0:1, 4), 1), c(rep(0:1, each = 4), 0))
rd <- ca_verify(plan0, 2)
re <- ca_verify(plan0, 2, v = c(3, 2, 2))
stopifnot(identical(rd$covered, re$covered), rd$gaps == re$gaps,
          rd$missing_tuples == re$missing_tuples,
          identical(rd$v, c(3L, 2L, 2L)))

## 5. 1-based auto-shift equivalence
r0 <- ca_verify(ca, 2)
r1 <- ca_verify(ca + 1L, 2)
stopifnot(identical(r0$covered, r1$covered), r0$colsets == r1$colsets)

## 6. randomized cross-validation against the oracle
set.seed(20260810)
for (i in 1:25) {
    N <- sample(4:9, 1); k <- sample(3:5, 1); v <- sample(2:3, 1)
    x <- matrix(sample(0:(v - 1), N * k, replace = TRUE), N, k)
    if (i %% 5 == 0) x[sample(length(x), 2)] <- NA
    r <- ca_verify(x, 2, v = v)
    b <- brute(x, 2, v)
    stopifnot(identical(isTRUE(r$covered), b$covered),
              r$gaps == b$gaps, r$missing_tuples == b$missing)
}

## 7. threads argument accepted (single- or multi- depending on OpenMP)
stopifnot(isTRUE(ca_verify(ca, 2, threads = 4)$covered))

cat("all caverify tests passed\n")

## 8. on-the-fly argument expressions (reported by U. Groemping under
## RStudio; must work in plain R - regression guard)
plan <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
r <- ca_verify(plan[-1, ], 2)
stopifnot(!r$covered)
r <- ca_verify(rbind(plan, plan)[-1, ], 2)
stopifnot(isTRUE(r$covered))

## 9. interruptibility smoke: a run large enough to span several batches
## completes normally (interrupt behavior itself is manual-test only)
set.seed(1)
big <- matrix(sample(0:1, 40 * 18, replace = TRUE), 40, 18)
r <- ca_verify(big, 4)
stopifnot(r$colsets == choose(18, 4))

## ---- mixed-level (0.2.0) ----

## pure-R brute-force oracle for per-column symbol counts
## (shares no logic with the C code)
brute_mixed <- function(x, t, vs) {
    k <- ncol(x); miss <- 0; gaps <- 0
    combs <- utils::combn(k, t)
    for (ci in seq_len(ncol(combs))) {
        cols <- combs[, ci]
        tuples <- as.matrix(expand.grid(lapply(cols, function(cc) 0:(vs[cc] - 1))))
        sub <- x[, cols, drop = FALSE]
        hit <- logical(nrow(tuples))
        for (r in seq_len(nrow(sub))) {
            row <- sub[r, ]
            if (anyNA(row)) next   ## don't-care: contributes nothing here
            ok <- rep(TRUE, nrow(tuples))
            for (j in seq_len(t)) ok <- ok & (tuples[, j] == row[j])
            hit <- hit | ok
        }
        if (!all(hit)) { gaps <- gaps + 1; miss <- miss + sum(!hit) }
    }
    list(covered = gaps == 0, gaps = gaps, missing = miss)
}

## 10. mixed full factorial 3 x 2 x 2 is covered at t = 2 and t = 3
mca <- as.matrix(expand.grid(0:2, 0:1, 0:1))
stopifnot(isTRUE(ca_verify(mca, 2, v = c(3, 2, 2))$covered),
          isTRUE(ca_verify(mca, 3, v = c(3, 2, 2))$covered))

## 11. the mixed-level example from the CAs::coverage documentation:
## one of the three 2-column projections misses exactly one tuple
plan <- cbind(rep(0:2, each = 3), c(rep(0:1, 4), 1), c(rep(0:1, each = 4), 0))
r <- ca_verify(plan, 2, v = c(3, 2, 2))
b <- brute_mixed(plan, 2, c(3, 2, 2))
stopifnot(!r$covered, r$gaps == 1, r$missing_tuples == 1,
          r$gaps == b$gaps, r$missing_tuples == b$missing)
## the reported example must really be missing from the array
ex <- r$examples[1, ]
cols <- ex[1:2]; val <- ex[3:4]
stopifnot(!any(plan[, cols[1]] == val[1] & plan[, cols[2]] == val[2]))

## 12. "auto" inference and explicit vector agree; result carries v vector
r_auto <- ca_verify(plan, 2, v = "auto")
stopifnot(identical(r_auto$covered, r$covered),
          r_auto$gaps == r$gaps,
          r_auto$missing_tuples == r$missing_tuples,
          identical(r_auto$v, c(3L, 2L, 2L)))

## 13. uniform array via vector v matches scalar v exactly
ca2 <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0))
rs <- ca_verify(ca2[-1, ], 2)
rv <- ca_verify(ca2[-1, ], 2, v = c(2, 2, 2))
stopifnot(identical(rs$covered, rv$covered), rs$gaps == rv$gaps,
          rs$missing_tuples == rv$missing_tuples,
          identical(rs$examples, rv$examples), identical(rv$v, 2L))

## 14. 1-based mixed input auto-shifts; examples come back 1-based
r1b <- ca_verify(plan + 1L, 2, v = c(3, 2, 2))
stopifnot(r1b$gaps == r$gaps, r1b$missing_tuples == r$missing_tuples,
          all(r1b$examples[, 3:4] == r$examples[, 3:4] + 1L))

## 15. mixed array: an NA row contributes nothing (don't-care
## semantics), so it cannot restore broken coverage
gap <- plan[-(1:2), ]                       ## break coverage
rgap <- ca_verify(gap, 2, v = c(3, 2, 2))
stopifnot(!rgap$covered)
fix <- rbind(gap, c(NA, NA, NA))
rfix <- ca_verify(fix, 2, v = c(3, 2, 2))
stopifnot(!rfix$covered, rfix$gaps == rgap$gaps,
          rfix$missing_tuples == rgap$missing_tuples)

## 16. randomized mixed-level cross-validation against the oracle,
## with and without NA wildcards, strengths 2 and 3
set.seed(20260823)
for (i in 1:40) {
    k <- sample(3:5, 1); t <- sample(2:min(3, k), 1)
    vs <- sample(2:4, k, replace = TRUE)
    N <- sample(4:12, 1)
    x <- sapply(vs, function(v) sample(0:(v - 1), N, replace = TRUE))
    x[1, 1] <- 0L   ## pin the coding to 0-based (shift detection is global)
    if (i %% 4 == 0) x[sample(length(x), 2)] <- NA
    r <- ca_verify(x, t, v = vs)
    b <- brute_mixed(x, t, vs)
    stopifnot(identical(isTRUE(r$covered), b$covered),
              r$gaps == b$gaps, r$missing_tuples == b$missing)
}

## 17. v vector length mismatch and bad values are rejected
stopifnot(inherits(try(ca_verify(plan, 2, v = c(3, 2)), silent = TRUE), "try-error"))
stopifnot(inherits(try(ca_verify(plan, 2, v = c(3, 2, 0)), silent = TRUE), "try-error"))

cat("all mixed-level caverify tests passed\n")
