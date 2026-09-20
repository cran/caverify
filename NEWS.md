# caverify 0.2.0

* The default `v = NULL` now infers per-column symbol counts from each
  column's own maximum, identical to `v = "auto"`, which remains as an
  explicit alias. Versions 0.1.x inferred one uniform value from the
  global data range, which misjudges mixed-level arrays; changed on a
  recommendation by Ulrike Groemping, since a default should not
  assume a uniform CA. Declare a scalar `v` where the check should
  also catch a column that fails to reach its intended symbol count.
* NA semantics corrected, a breaking change from 0.1.x. An NA marks a
  flexible ("don't care") entry: a row now contributes nothing to a
  projection in which it has an NA, so a verified array is covering no
  matter how its NA entries are later filled. The 0.1.x behaviour let
  one NA count as every symbol at once, per tuple independently, and
  could certify arrays that no single choice of values would make
  covering (a 4-run array passed as a strength-4 covering array on 21
  four-level columns). Thanks to Ulrike Groemping for the
  counterexample. The new behaviour also agrees exactly with
  `CAs::coverage()` on arrays with NAs.
* Mixed-level covering arrays (MCAs): `v` now also accepts an integer
  vector of length `ncol(x)` giving each column its own number of
  symbols, or the string `"auto"` to infer per-column symbol counts
  from each column's own maximum. Tuple counting and indexing in the C
  kernel are mixed-radix; a uniform array is the degenerate case and
  takes the same code path.
* A scalar `v` behaves exactly as in 0.1.x. The only change for
  existing callers is the `v = NULL` default described above, which
  now reads each column's own symbol count. A uniform array with every
  symbol present in every column verifies identically under both
  readings.
* `print` shows mixed-level results with the levels profile in
  exponent notation (e.g. `levels 4^2 3 2^3`).
* Tests: mixed-level brute-force oracle plus randomized mixed-level
  cross-validation added.

# caverify 0.1.3

* Fixed an installation failure on R-devel with clang 22, seen on the
  CRAN flavour r-devel-linux-x86_64-fedora-clang. The R headers remap
  `match` to `Rf_match` unless `R_NO_REMAP` is defined, and clang 22's
  `omp.h` uses `match` as a clause of `#pragma omp declare variant`, so
  an `omp.h` included after the R headers no longer parsed. `omp.h` is
  now included first and the C code compiles with `R_NO_REMAP` and the
  `Rf_` prefixed API. No user-visible change.

# caverify 0.1.2

* Automatic thread selection: by default the checker now uses half the
  machine's logical cores (single-threaded for small jobs, capped to 2
  during CRAN checks). No setup needed; `options(caverify.threads = n)`
  or the `threads` argument still override.

# caverify 0.1.1

* Long verifications are now interruptible (Escape / Ctrl-C) without
  aborting the R session: combinations are processed in batches with an
  interrupt check between batches, and all working memory is R-managed
  so interruption cannot leak.
* `threads` now defaults to `getOption("caverify.threads", 1L)`.
* DESCRIPTION metadata cleanup (URL, BugReports).

# caverify 0.1.0

* First version: registered .Call interface, NA as wildcard ("flexible
  value"), optional OpenMP threading, 0- or 1-based symbol detection,
  tests including randomized cross-validation against a pure R oracle.
