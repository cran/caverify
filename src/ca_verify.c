/* ca_verify.c - fast strength-t coverage checker for covering arrays,
 * R .Call interface.
 *
 * Port of the standalone ca_verify CLI tool (C. Smolen, ca-tools) to a
 * native R entry point: the array arrives as an R integer matrix (no
 * files, no system calls), R's NA_integer_ marks a flexible
 * ("don't care") entry, and results return as an R list.
 *
 * Mixed-level support (0.2.0): the number of symbols is a per-column
 * vector vs[0..k-1]; tuple ranks within a column set are mixed-radix
 * numbers (last column of the set changing fastest, radix vs[c[j]] in
 * position j). A uniform array is the degenerate case vs[j] == v for
 * all j and takes the same code path.
 *
 * Threading: optional OpenMP (portable through R's SHLIB_OPENMP_CFLAGS
 * mechanism, including Windows/Rtools); falls back to single-threaded
 * cleanly when OpenMP is unavailable. No R API calls occur inside the
 * parallel region (R is not thread-safe); all SEXP construction happens
 * after the join.
 *
 * Symbols are expected as 0..v-1 (the R wrapper auto-shifts 1..v input).
 */
/* Include order matters here. R's headers remap much of the API into the
 * global namespace unless R_NO_REMAP is defined, and one of those macros
 * is match -> Rf_match. Clang 22's omp.h uses
 *   #pragma omp declare variant ... match(...)
 * so an omp.h seen after the R headers has its match clause rewritten and
 * no longer parses (r-devel-linux-x86_64-fedora-clang, 2026-08-20).
 * R_NO_REMAP removes the remap entirely; including omp.h first is kept as
 * a second line of defence against any other collision of the same kind. */
#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define MAX_T 20

typedef struct {
    long long combos;
    long long gaps;
    long long missing;
    int n_ex;
    int max_ex;
    int *ex_cols;          /* max_ex x T, 0-based column indices */
    long long *ex_tuple;   /* max_ex tuple ranks                 */
} Acc;

/* advance combination c[0..t-1] over k columns; return 0 when exhausted */
static int next_comb(int *c, int k, int t) {
    int i = t - 1;
    while (i >= 0 && c[i] == k - t + i) i--;
    if (i < 0) return 0;
    c[i]++;
    for (int j = i + 1; j < t; j++) c[j] = c[j - 1] + 1;
    return 1;
}

/* mark the tuple exhibited by row r on columns c[].
 * NA ("don't care" / flexible value) semantics, 0.2.0: a row with an
 * NA inside the projection contributes NOTHING to that projection.
 * Coverage must hold on concrete entries alone, so the array stays
 * covering no matter how the NAs are later filled. (The 0.1.x
 * per-tuple wildcard credited one NA with every symbol at once,
 * which certifies arrays no single filling could make covering;
 * counterexample due to Ulrike Groemping, 2026-08-23.) */
static void mark_row(uint64_t *bm, const int *m, int nrow, int r,
                     const int *c, int t, const long long *mult) {
    long long base = 0;
    for (int j = 0; j < t; j++) {
        int s = m[r + (long long)c[j] * nrow];
        if (s == NA_INTEGER) return;   /* contributes nothing here */
        base += (long long)s * mult[j];
    }
    bm[base >> 6] |= 1ULL << (base & 63);
}

SEXP C_ca_verify(SEXP mat, SEXP t_, SEXP v_, SEXP nthreads_, SEXP maxreport_)
{
    if (!Rf_isInteger(mat) || !Rf_isMatrix(mat))
        Rf_error("internal: matrix of integers expected");
    const int *m = INTEGER(mat);
    const int nrow = Rf_nrows(mat), ncol = Rf_ncols(mat);
    const int t = Rf_asInteger(t_);
    int nt = Rf_asInteger(nthreads_);
    const int max_ex = Rf_asInteger(maxreport_);

    if (t < 1 || t > MAX_T) Rf_error("t must be in 1..%d", MAX_T);
    if (ncol < t) Rf_error("array has fewer columns (%d) than t (%d)", ncol, t);
    if (!Rf_isInteger(v_) || Rf_length(v_) != ncol)
        Rf_error("internal: 'v' must be an integer vector of length ncol");
    const int *vs = INTEGER(v_);
    for (int j = 0; j < ncol; j++) {
        if (vs[j] == NA_INTEGER || vs[j] < 1)
            Rf_error("number of symbols for column %d must be >= 1", j + 1);
    }
    if (nt < 1) nt = 1;
    if (nt > 64) nt = 64;

    /* worst-case tuples per column set: product of the t largest vs */
    int *vsrt = (int *) R_alloc((size_t) ncol, sizeof(int));
    memcpy(vsrt, vs, (size_t) ncol * sizeof(int));
    for (int i = 0; i < ncol - 1; i++)          /* descending selection sort */
        for (int j = i + 1; j < ncol; j++)
            if (vsrt[j] > vsrt[i]) { int h = vsrt[i]; vsrt[i] = vsrt[j]; vsrt[j] = h; }
    long long vt_max = 1;
    for (int i = 0; i < t; i++) {
        vt_max *= vsrt[i];
        if (vt_max > (1LL << 31))
            Rf_error("tuples per column set too many (limit 2^31)");
    }

    /* validate symbol range once, outside the parallel region */
    for (int j = 0; j < ncol; j++) {
        const int *col = m + (long long) j * nrow;
        for (int r = 0; r < nrow; r++) {
            int s = col[r];
            if (s != NA_INTEGER && (s < 0 || s >= vs[j]))
                Rf_error("symbol %d in column %d out of range 0..%d (after any auto-shift)",
                         s, j + 1, vs[j] - 1);
        }
    }

#ifndef _OPENMP
    nt = 1;
#endif

    /* All allocations via R_alloc: reclaimed automatically if an
     * interrupt (R_CheckUserInterrupt longjmp) fires between batches,
     * so interruption is leak-free. (v0.1.1 interruptibility fix.) */
    Acc *acc = (Acc *) R_alloc((size_t) nt, sizeof(Acc));
    memset(acc, 0, (size_t) nt * sizeof(Acc));
    for (int i = 0; i < nt; i++) {
        acc[i].max_ex = max_ex;
        acc[i].ex_cols = (int *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1) * t,
                                         sizeof(int));
        acc[i].ex_tuple = (long long *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1),
                                                sizeof(long long));
    }
    const long long words = (vt_max + 63) >> 6;
    uint64_t *bms = (uint64_t *) R_alloc((size_t) nt * words, sizeof(uint64_t));

    /* Batched enumeration: combinations are collected serially into a
     * buffer, processed in parallel, and R_CheckUserInterrupt() runs
     * BETWEEN batches (never inside the parallel region - the R API is
     * not thread-safe). Escape/Ctrl-C therefore works without aborting
     * the R session. */
#define CA_BATCH 4096
    int *buf = (int *) R_alloc((size_t) CA_BATCH * t, sizeof(int));
    int c[MAX_T];
    for (int j = 0; j < t; j++) c[j] = j;
    int done = 0;
    while (!done) {
        int nb = 0;
        while (nb < CA_BATCH) {
            memcpy(buf + (size_t) nb * t, c, t * sizeof(int));
            nb++;
            if (!next_comb(c, ncol, t)) { done = 1; break; }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic)
#endif
        for (int bi = 0; bi < nb; bi++) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            uint64_t *bm = bms + (size_t) tid * words;
            const int *cc = buf + (size_t) bi * t;
            /* mixed-radix multipliers for this column set:
             * last column changes fastest */
            long long mult[MAX_T], mm = 1, vt = 1;
            for (int j = t - 1; j >= 0; j--) { mult[j] = mm; mm *= vs[cc[j]]; }
            vt = mm;
            const long long cwords = (vt + 63) >> 6;
            memset(bm, 0, (size_t) cwords * sizeof(uint64_t));
            for (int r = 0; r < nrow; r++)
                mark_row(bm, m, nrow, r, cc, t, mult);
            long long covered = 0;
            for (long long i = 0; i < cwords; i++)
                covered += __builtin_popcountll(bm[i]);
            Acc *w = &acc[tid];
            w->combos++;
            if (covered != vt) {
                w->gaps++;
                w->missing += vt - covered;
                if (w->n_ex < w->max_ex) {
                    for (long long idx = 0; idx < vt; idx++) {
                        if (!(bm[idx >> 6] & (1ULL << (idx & 63)))) {
                            memcpy(w->ex_cols + (size_t) w->n_ex * t,
                                   cc, t * sizeof(int));
                            w->ex_tuple[w->n_ex] = idx;
                            w->n_ex++;
                            break;
                        }
                    }
                }
            }
        }
        R_CheckUserInterrupt();
    }

    long long combos = 0, gaps = 0, missing = 0;
    int n_ex = 0;
    for (int i = 0; i < nt; i++) {
        combos += acc[i].combos;
        gaps += acc[i].gaps;
        missing += acc[i].missing;
        n_ex += acc[i].n_ex;
    }
    if (n_ex > max_ex) n_ex = max_ex;

    /* examples: n_ex x (t cols 1-based, then t tuple symbols 0-based) */
    SEXP ex = PROTECT(Rf_allocMatrix(INTSXP, n_ex, 2 * t));
    int *e = INTEGER(ex);
    int row = 0;
    for (int i = 0; i < nt && row < n_ex; i++) {
        for (int j = 0; j < acc[i].n_ex && row < n_ex; j++, row++) {
            for (int q = 0; q < t; q++)
                e[row + (long long)q * n_ex] =
                    acc[i].ex_cols[(size_t)j * t + q] + 1;
            long long idx = acc[i].ex_tuple[j];
            int digs[MAX_T];
            for (int q = t - 1; q >= 0; q--) {
                int vq = vs[acc[i].ex_cols[(size_t)j * t + q]];
                digs[q] = (int)(idx % vq); idx /= vq;
            }
            for (int q = 0; q < t; q++)
                e[row + (long long)(t + q) * n_ex] = digs[q];
        }
    }
    /* R_alloc memory is released by R automatically. */

    const char *names[] = {"covered", "colsets", "gaps", "missing_tuples",
                           "examples", ""};
    SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(gaps == 0));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal((double)combos));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal((double)gaps));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal((double)missing));
    SET_VECTOR_ELT(out, 4, ex);
    UNPROTECT(2);
    return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_ca_verify", (DL_FUNC) &C_ca_verify, 5},
    {NULL, NULL, 0}
};

void R_init_caverify(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
