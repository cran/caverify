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
            ok <- rep(TRUE, nrow(tuples))
            for (j in seq_len(t)) {
                if (!is.na(row[j])) ok <- ok & (tuples[, j] == row[j])
            }
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

## 4. wildcard row restores coverage
wc <- rbind(br, c(NA, NA, NA))
stopifnot(isTRUE(ca_verify(wc, 2)$covered))

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
