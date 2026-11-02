# caverify 0.3.0

* A second exact method, used by default where it fits in memory. Every
  set of t columns is written as its t-1 smallest columns plus its
  largest column. For each set of t-1 columns the rows are grouped by
  their values on those columns, and one pass over each group, one bit
  per column, decides the coverage of all the later columns at once. On
  the strength-4 binary arrays used in testing it was about 40 to 50
  times faster than 0.2.0 at a few hundred columns, and the gain grows
  with the number of columns. The method is the one of `ca_verify.c`
  (version 2) in ca-tools, described in `tools/RESEARCH_ca_verify.md`
  there; this is a portable rewrite of it for R, with mixed levels and
  the NA rule of 0.2.0.
* The method of 0.2.0 is kept and is used where the new one would need
  more than `getOption("caverify.maxmem", 2^30)` bytes or a column has
  more than 64 symbols. `options(caverify.kernel = "direct")` or
  `"prefix"` forces one of them, for testing and comparison. The result
  names the method used in a new element `kernel`.
* Results are unchanged: `covered`, `colsets`, `gaps` and
  `missing_tuples` agree with 0.2.0 on every array tested (random
  uniform, mixed-level and binary arrays with and without NAs, strengths
  1 to 6, and covering arrays from 'CAs' with a row removed, an entry
  changed or an entry set to NA).
* The examples are now the first column sets with a gap in
  lexicographic order, each with its smallest missing value combination,
  for any number of threads. In 0.2.0 this held for one thread; with
  several threads the examples chosen depended on scheduling.

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
