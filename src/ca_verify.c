/* ca_verify.c - fast strength-t coverage checker for covering arrays,
 * R .Call interface.
 *
 * The array arrives as an R integer matrix (no files, no system calls),
 * R's NA_integer_ marks a flexible ("don't care") entry, and results
 * return as an R list.
 *
 * Mixed-level support (0.2.0): the number of symbols is a per-column
 * vector vs[0..k-1]; tuple ranks within a column set are mixed-radix
 * numbers (last column of the set changing fastest, radix vs[c[j]] in
 * position j). A uniform array is the degenerate case.
 *
 * Two exact kernels (0.3.0). Both give identical results; the first is
 * kept as a fallback and as an independent reference for the tests.
 *
 *  "direct"  (the kernel of 0.2.0): every t-set of columns is checked
 *            on its own, one bitmap of tuples per column set.
 *
 *  "prefix"  (new in 0.3.0, the method of ca_verify version 2 in
 *            ca-tools): every t-set is written uniquely as a (t-1)-set P
 *            (its t-1 smallest columns) plus its largest column d. For a
 *            fixed P the rows are grouped by their value tuple on P into
 *            classes. The t-set P + {d} covers tuple (s, u) iff some row
 *            of class s has value u in column d. With one bit plane per
 *            row and value (bit d set iff the row has value u in column
 *            d), the OR of the planes of the rows of class s has bit d
 *            clear exactly where P + {d} misses (s, u). This checks all
 *            candidate columns d at once, 64 per machine word, so the
 *            work per (t-1)-set is about N * v * k / 64 word operations
 *            instead of N * t operations for each of its k t-sets.
 *            Binary arrays without NAs use one plane per row with an OR
 *            and an AND accumulator.
 *
 * NA ("don't care") semantics, both kernels: a row with an NA inside a
 * projection contributes nothing to that projection. In the prefix
 * kernel a row with an NA in P belongs to no class, and a row with an
 * NA in column d has no bit for d in any plane.
 *
 * Results: colsets, gaps (column sets with at least one missing tuple),
 * missing_tuples, and up to `report` examples. The examples are the
 * first column sets with a gap in lexicographic order, each with its
 * smallest missing tuple, whatever the number of threads. (In 0.2.0
 * this held for one thread; with several threads the selection
 * depended on scheduling.)
 *
 * Threading: optional OpenMP (portable through R's SHLIB_OPENMP_CFLAGS
 * mechanism); single-threaded when OpenMP is unavailable. No R API call
 * occurs inside a parallel region; R_CheckUserInterrupt() runs between
 * batches, and all memory comes from R_alloc, so an interrupt is
 * leak-free.
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

static int popc(uint64_t x) { return __builtin_popcountll(x); }
static int ctz(uint64_t x) { return __builtin_ctzll(x); }

/* ================================================================== */
/* direct kernel (0.2.0)                                               */
/* ================================================================== */

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

/* mark the tuple exhibited by row r on columns c[]; a row with an NA in
 * the projection contributes nothing to it */
static void mark_row(uint64_t *bm, const int *m, int nrow, int r,
                     const int *c, int t, const long long *mult) {
    long long base = 0;
    for (int j = 0; j < t; j++) {
        int s = m[r + (long long)c[j] * nrow];
        if (s == NA_INTEGER) return;
        base += (long long)s * mult[j];
    }
    bm[base >> 6] |= 1ULL << (base & 63);
}

/* Output of either kernel: totals, and examples as (cols[t], tuple rank)
 * in lexicographic order of the column sets. */
typedef struct {
    double combos, gaps, missing;
    int n_ex;
    int *ex_cols;          /* max_ex x t */
    long long *ex_tuple;
} Result;

/* The direct kernel processes column sets in batches of CA_BATCH in
 * lexicographic order; within a batch the examples of each column set
 * are stored by batch position, so the merged list is in lexicographic
 * order for any thread count. */
#define CA_BATCH 4096

static void kernel_direct(const int *m, int nrow, int ncol, int t,
                          const int *vs, long long vt_max, int nt,
                          int max_ex, Result *res)
{
    const long long words = (vt_max + 63) >> 6;
    uint64_t *bms = (uint64_t *) R_alloc((size_t) nt * words, sizeof(uint64_t));
    int *buf = (int *) R_alloc((size_t) CA_BATCH * t, sizeof(int));
    long long *firstmiss = (long long *) R_alloc(CA_BATCH, sizeof(long long));
    double *bmiss = (double *) R_alloc(CA_BATCH, sizeof(double));
    int c[MAX_T];
    for (int j = 0; j < t; j++) c[j] = j;
    int done = 0;
    res->combos = res->gaps = res->missing = 0;
    res->n_ex = 0;
    while (!done) {
        int nb = 0;
        while (nb < CA_BATCH) {
            memcpy(buf + (size_t) nb * t, c, t * sizeof(int));
            nb++;
            if (!next_comb(c, ncol, t)) { done = 1; break; }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
#endif
        for (int bi = 0; bi < nb; bi++) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            uint64_t *bm = bms + (size_t) tid * words;
            const int *cc = buf + (size_t) bi * t;
            long long mult[MAX_T], mm = 1;
            for (int j = t - 1; j >= 0; j--) { mult[j] = mm; mm *= vs[cc[j]]; }
            const long long vt = mm;
            const long long cwords = (vt + 63) >> 6;
            memset(bm, 0, (size_t) cwords * sizeof(uint64_t));
            for (int r = 0; r < nrow; r++)
                mark_row(bm, m, nrow, r, cc, t, mult);
            long long covered = 0;
            for (long long i = 0; i < cwords; i++) covered += popc(bm[i]);
            bmiss[bi] = (double)(vt - covered);
            firstmiss[bi] = -1;
            if (covered != vt) {
                for (long long idx = 0; idx < vt; idx++)
                    if (!(bm[idx >> 6] & (1ULL << (idx & 63)))) {
                        firstmiss[bi] = idx; break;
                    }
            }
        }
        for (int bi = 0; bi < nb; bi++) {
            res->combos += 1;
            if (firstmiss[bi] >= 0) {
                res->gaps += 1;
                res->missing += bmiss[bi];
                if (res->n_ex < max_ex) {
                    memcpy(res->ex_cols + (size_t) res->n_ex * t,
                           buf + (size_t) bi * t, t * sizeof(int));
                    res->ex_tuple[res->n_ex] = firstmiss[bi];
                    res->n_ex++;
                }
            }
        }
        R_CheckUserInterrupt();
    }
}

/* ================================================================== */
/* prefix kernel (0.3.0)                                               */
/* ================================================================== */

#define CH 4                /* words per column chunk (256 columns) */

typedef struct {
    /* shared, read-only */
    const int *m; int N, K, t, vmax, Kw, NCH, binary;
    const int *vs;
    const uint64_t *PL;     /* planes, chunk-major: ((u*NCH+ci)*N + r)*CH + x */
    const uint64_t *vmask;  /* vmax x Kw: bit d set iff vs[d] > u */
    const uint64_t *kmask;  /* Kw: bit d set iff d < K */
} Shared;

typedef struct {
    double combos, gaps, missing;
    int *cid;               /* (t) x N class ids per prefix depth, -1 = none */
    uint64_t *keys;         /* N sort keys */
    int *rows;              /* N row indices grouped by class */
    long long *cls;         /* class id of each group */
    int *start;             /* N+1 group starts */
    int *cnt;               /* counting-sort buckets, 4N+2 */
    uint64_t *gapmask, *found, *dm;   /* Kw each */
    uint64_t *acc;          /* (vmax+1) x CH */
    int *pd; long long *ptuple; int np;   /* first gaps of this prefix */
    /* examples of the current unit */
    int max_ex, n_ex;
    int *ex_cols; long long *ex_tuple;
    int p[MAX_T];
} PWork;

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

/* d gets its first missing tuple (s, u) */
static void note_first(PWork *w, int d, long long rank) {
    w->pd[w->np] = d; w->ptuple[w->np] = rank; w->np++;
}

/* class s has misses umiss[u][x] (already masked to valid d) for the
 * words wb .. wb+nw-1 */
static void fold_misses(const Shared *S, PWork *w, long long s, int wb, int nw,
                        const uint64_t *umiss) {
    for (int x = 0; x < nw; x++) {
        const int wi = wb + x;
        uint64_t any = 0;
        for (int u = 0; u < S->vmax; u++) {
            uint64_t mu = umiss[u * CH + x];
            any |= mu;
            w->missing += popc(mu);
        }
        if (!any) continue;
        w->gapmask[wi] |= any;
        uint64_t nb = any & ~w->found[wi];
        if (!nb) continue;
        w->found[wi] |= nb;
        if (w->n_ex >= w->max_ex) continue;
        while (nb) {
            const int b = ctz(nb); nb &= nb - 1;
            const int d = wi * 64 + b;
            int u = 0;
            while (!((umiss[u * CH + x] >> b) & 1)) u++;
            note_first(w, d, s * S->vs[d] + u);
        }
    }
}

/* an empty class s: every value of every valid d is missing */
static void empty_class(const Shared *S, PWork *w, long long s, int w0) {
    for (int wi = w0; wi < S->Kw; wi++) {
        uint64_t un = w->dm[wi];
        if (!un) continue;
        for (int u = 0; u < S->vmax; u++)
            w->missing += popc(un & S->vmask[(size_t)u * S->Kw + wi]);
        w->gapmask[wi] |= un;
        uint64_t nb = un & ~w->found[wi];
        if (!nb) continue;
        w->found[wi] |= nb;
        if (w->n_ex >= w->max_ex) continue;
        while (nb) {
            const int b = ctz(nb); nb &= nb - 1;
            const int d = wi * 64 + b;
            note_first(w, d, s * S->vs[d]);
        }
    }
}

/* one nonempty class s with rows L[0..ml-1]: words from chunk c0 on */
static void class_rows(const Shared *S, PWork *w, long long s,
                       const int *L, int ml, int c0) {
    const int N = S->N, NCH = S->NCH, Kw = S->Kw, vmax = S->vmax;
    const size_t cs = (size_t) N * CH;
    uint64_t um[64 * CH];   /* vmax <= 64 for this kernel */
    for (int ci = c0; ci < NCH; ci++) {
        const int wb = ci * CH;
        const int nw = Kw - wb < CH ? Kw - wb : CH;
        uint64_t badany = 0;
        if (S->binary) {
            uint64_t o[CH] = {0, 0, 0, 0};
            uint64_t a[CH] = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
            const uint64_t *pb = S->PL + (size_t) ci * cs;
            for (int i = 0; i < ml; i++) {
                const uint64_t *p = pb + (size_t) L[i] * CH;
                for (int x = 0; x < CH; x++) { o[x] |= p[x]; a[x] &= p[x]; }
            }
            for (int x = 0; x < nw; x++) {
                const uint64_t dmx = w->dm[wb + x];
                um[0 * CH + x] = a[x] & dmx;    /* value 0 missing: all rows 1 */
                um[1 * CH + x] = ~o[x] & dmx;   /* value 1 missing: no row 1 */
                badany |= um[x] | um[CH + x];
            }
        } else {
            for (int u = 0; u < vmax; u++) {
                uint64_t o[CH] = {0, 0, 0, 0};
                const uint64_t *pb = S->PL + ((size_t) u * NCH + ci) * cs;
                for (int i = 0; i < ml; i++) {
                    const uint64_t *p = pb + (size_t) L[i] * CH;
                    for (int x = 0; x < CH; x++) o[x] |= p[x];
                }
                for (int x = 0; x < nw; x++) {
                    const uint64_t mu = ~o[x] & w->dm[wb + x]
                                      & S->vmask[(size_t) u * Kw + wb + x];
                    um[u * CH + x] = mu;
                    badany |= mu;
                }
            }
        }
        if (badany) fold_misses(S, w, s, wb, nw, um);
    }
}

/* all t-sets P + {d}, d > c = max(P); class ids in w->cid + (t-1)*N */
static void process_prefix(const Shared *S, PWork *w, int c) {
    const int N = S->N, Kw = S->Kw, t = S->t;
    const int *cid = w->cid + (size_t)(t - 1) * N;
    const int w0 = (c + 1) >> 6;
    long long ncls = 1;
    for (int j = 0; j < t - 1; j++) ncls *= S->vs[w->p[j]];

    /* candidate columns d > c */
    for (int wi = 0; wi < Kw; wi++) {
        uint64_t mk = S->kmask[wi];
        if (wi < w0) mk = 0;
        else if (wi == w0 && ((c + 1) & 63)) mk &= ~0ULL << ((c + 1) & 63);
        w->dm[wi] = mk;
        w->gapmask[wi] = 0;
        w->found[wi] = 0;
    }
    w->combos += (double)(S->K - 1 - c);
    w->np = 0;

    /* group rows by class id, classes ascending: counting sort when the
     * number of classes is small against N, otherwise a sort of keys */
    int ng = 0, n = 0;
    if (ncls <= 4LL * N + 1) {
        int *cnt = w->cnt;
        memset(cnt, 0, (size_t)(ncls + 1) * sizeof(int));
        for (int r = 0; r < N; r++) if (cid[r] >= 0) { cnt[cid[r] + 1]++; n++; }
        for (long long s = 0; s < ncls; s++) {
            if (cnt[s + 1]) { w->cls[ng] = s; w->start[ng] = cnt[s]; ng++; }
            cnt[s + 1] += cnt[s];
        }
        for (int r = 0; r < N; r++) if (cid[r] >= 0) w->rows[cnt[cid[r]]++] = r;
    } else {
        for (int r = 0; r < N; r++)
            if (cid[r] >= 0) w->keys[n++] = ((uint64_t) cid[r] << 32) | (uint64_t) r;
        qsort(w->keys, (size_t) n, sizeof(uint64_t), cmp_u64);
        for (int i = 0; i < n; i++) {
            const long long k = (long long)(w->keys[i] >> 32);
            if (i == 0 || k != w->cls[ng - 1]) { w->cls[ng] = k; w->start[ng] = i; ng++; }
            w->rows[i] = (int)(w->keys[i] & 0xffffffffULL);
        }
    }
    w->start[ng] = n;

    /* classes in ascending order; the empty ones that come before a
     * nonempty one matter only through the smallest of them (for the
     * examples) and their number (for the counts) */
    const int c0 = w0 / CH;
    long long next_s = 0;           /* smallest class id not yet visited */
    int empty_done = 0;
    for (int g = 0; g <= ng; g++) {
        const long long s = g < ng ? w->cls[g] : ncls;
        if (s > next_s) {           /* classes next_s .. s-1 are empty */
            long long more = s - next_s;
            if (!empty_done) { empty_class(S, w, next_s, w0); empty_done = 1; more--; }
            /* counts for the other empty classes of this run */
            if (more > 0) {
                double per = 0;
                for (int wi = w0; wi < Kw; wi++)
                    for (int u = 0; u < S->vmax; u++)
                        per += popc(w->dm[wi] & S->vmask[(size_t) u * Kw + wi]);
                w->missing += per * (double) more;
            }
        }
        if (g == ng) break;
        class_rows(S, w, s, w->rows + w->start[g], w->start[g + 1] - w->start[g], c0);
        next_s = s + 1;
    }

    double gg = 0;
    for (int wi = w0; wi < Kw; wi++) gg += popc(w->gapmask[wi]);
    w->gaps += gg;

    /* examples of this prefix, in order of d */
    if (w->np && w->n_ex < w->max_ex) {
        for (int i = 1; i < w->np; i++) {
            int d = w->pd[i]; long long tp = w->ptuple[i]; int j = i - 1;
            while (j >= 0 && w->pd[j] > d) {
                w->pd[j + 1] = w->pd[j]; w->ptuple[j + 1] = w->ptuple[j]; j--;
            }
            w->pd[j + 1] = d; w->ptuple[j + 1] = tp;
        }
        for (int i = 0; i < w->np && w->n_ex < w->max_ex; i++) {
            int *ec = w->ex_cols + (size_t) w->n_ex * t;
            for (int j = 0; j < t - 1; j++) ec[j] = w->p[j];
            ec[t - 1] = w->pd[i];
            w->ex_tuple[w->n_ex] = w->ptuple[i];
            w->n_ex++;
        }
    }
    w->np = 0;
}

/* depth-first over prefixes: p[0..j-1] fixed, class ids at depth j in
 * w->cid + j*N */
static void dfs(const Shared *S, PWork *w, int j) {
    const int N = S->N, K = S->K, t = S->t;
    if (j == t - 1) { process_prefix(S, w, w->p[t - 2]); return; }
    const int lo = w->p[j - 1] + 1, hi = K - t + j;
    const int *prev = w->cid + (size_t) j * N;
    int *cur = w->cid + (size_t)(j + 1) * N;
    for (int c = lo; c <= hi; c++) {
        w->p[j] = c;
        const int *col = S->m + (size_t) c * N;
        const int vc = S->vs[c];
        for (int r = 0; r < N; r++) {
            const int a = col[r];
            cur[r] = (prev[r] < 0 || a == NA_INTEGER) ? -1 : prev[r] * vc + a;
        }
        dfs(S, w, j + 1);
    }
}

/* one unit = all prefixes with smallest column c0 (t >= 2), or the
 * single empty prefix (t == 1, c0 == -1) */
static void run_unit(const Shared *S, PWork *w, int c0) {
    const int N = S->N, t = S->t;
    if (t == 1) {
        for (int r = 0; r < N; r++) w->cid[r] = 0;
        process_prefix(S, w, -1);
        return;
    }
    w->p[0] = c0;
    const int *col = S->m + (size_t) c0 * N;
    int *cur = w->cid + (size_t) N;
    for (int r = 0; r < N; r++) cur[r] = col[r] == NA_INTEGER ? -1 : col[r];
    dfs(S, w, 1);
}

/* bytes the prefix kernel would need (shared planes plus per thread) */
static double prefix_bytes(int N, int K, int t, int vmax, int nt) {
    const double Kw = (K + 63) / 64, NCH = (Kw + CH - 1) / CH;
    const double planes = (double)(vmax == 2 ? 1 : vmax) * NCH * N * CH * 8;
    const double per = (double) t * N * 4 + N * 8.0 + N * 4.0 + (N + 1) * 12.0
                     + Kw * 8 * 3 + (double) K * 12;
    return planes + per * nt;
}

static void kernel_prefix(const int *m, int nrow, int ncol, int t,
                          const int *vs, int vmax, int nt, int max_ex,
                          Result *res)
{
    Shared S;
    S.m = m; S.N = nrow; S.K = ncol; S.t = t; S.vs = vs; S.vmax = vmax;
    S.Kw = (ncol + 63) / 64;
    S.NCH = (S.Kw + CH - 1) / CH;
    int hasna = 0;
    for (long long i = 0; i < (long long) nrow * ncol; i++)
        if (m[i] == NA_INTEGER) { hasna = 1; break; }
    S.binary = (vmax == 2 && !hasna);
    for (int d = 0; d < ncol; d++) if (vs[d] != 2) S.binary = 0;

    const int Kw = S.Kw, NCH = S.NCH, N = nrow;
    const size_t cs = (size_t) N * CH;
    const int np_ = S.binary ? 1 : vmax;
    uint64_t *PL = (uint64_t *) R_alloc((size_t) np_ * NCH * cs, sizeof(uint64_t));
    memset(PL, 0, (size_t) np_ * NCH * cs * sizeof(uint64_t));
    for (int d = 0; d < ncol; d++) {
        const int ci = (d >> 6) / CH, x = (d >> 6) % CH;
        const uint64_t bit = 1ULL << (d & 63);
        const int *col = m + (size_t) d * N;
        for (int r = 0; r < N; r++) {
            const int a = col[r];
            if (a == NA_INTEGER) continue;
            if (S.binary) { if (a) PL[(size_t) ci * cs + (size_t) r * CH + x] |= bit; }
            else PL[((size_t) a * NCH + ci) * cs + (size_t) r * CH + x] |= bit;
        }
    }
    S.PL = PL;
    uint64_t *vmask = (uint64_t *) R_alloc((size_t) vmax * Kw, sizeof(uint64_t));
    memset(vmask, 0, (size_t) vmax * Kw * sizeof(uint64_t));
    for (int d = 0; d < ncol; d++)
        for (int u = 0; u < vs[d]; u++)
            vmask[(size_t) u * Kw + (d >> 6)] |= 1ULL << (d & 63);
    S.vmask = vmask;
    uint64_t *kmask = (uint64_t *) R_alloc((size_t) Kw, sizeof(uint64_t));
    for (int wi = 0; wi < Kw; wi++) {
        const int lo = wi * 64, rem = ncol - lo;
        kmask[wi] = rem >= 64 ? ~0ULL : ((1ULL << rem) - 1);
    }
    S.kmask = kmask;

    /* units: smallest column c0 = 0 .. K - t (t >= 2), or one unit */
    const int nunits = t == 1 ? 1 : ncol - t + 1;
    int B = 2 * nt; if (B < 1) B = 1;
    if (B > nunits) B = nunits;
    PWork *W = (PWork *) R_alloc((size_t) B, sizeof(PWork));
    PWork *slot = W;    /* one work area per unit of the batch */
    for (int i = 0; i < B; i++) {
        PWork *w = &slot[i];
        memset(w, 0, sizeof(PWork));
        w->cid = (int *) R_alloc((size_t) t * N + 1, sizeof(int));
        w->keys = (uint64_t *) R_alloc((size_t) N + 1, sizeof(uint64_t));
        w->rows = (int *) R_alloc((size_t) N + 1, sizeof(int));
        w->cls = (long long *) R_alloc((size_t) N + 1, sizeof(long long));
        w->start = (int *) R_alloc((size_t) N + 2, sizeof(int));
        w->cnt = (int *) R_alloc((size_t) 4 * N + 3, sizeof(int));
        w->gapmask = (uint64_t *) R_alloc((size_t) Kw, sizeof(uint64_t));
        w->found = (uint64_t *) R_alloc((size_t) Kw, sizeof(uint64_t));
        w->dm = (uint64_t *) R_alloc((size_t) Kw, sizeof(uint64_t));
        w->pd = (int *) R_alloc((size_t) ncol + 1, sizeof(int));
        w->ptuple = (long long *) R_alloc((size_t) ncol + 1, sizeof(long long));
        w->max_ex = max_ex;
        w->ex_cols = (int *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1) * t, sizeof(int));
        w->ex_tuple = (long long *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1), sizeof(long long));
    }

    res->combos = res->gaps = res->missing = 0;
    res->n_ex = 0;
    for (int u0 = 0; u0 < nunits; u0 += B) {
        const int nb = nunits - u0 < B ? nunits - u0 : B;
        const int need = max_ex - res->n_ex;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
#endif
        for (int bi = 0; bi < nb; bi++) {
            PWork *w = &slot[bi];
            w->combos = w->gaps = w->missing = 0;
            w->n_ex = 0;
            w->max_ex = need > 0 ? need : 0;
            run_unit(&S, w, t == 1 ? -1 : u0 + bi);
        }
        /* merge in unit order = lexicographic order of the t-sets */
        for (int bi = 0; bi < nb; bi++) {
            PWork *w = &slot[bi];
            res->combos += w->combos;
            res->gaps += w->gaps;
            res->missing += w->missing;
            for (int i = 0; i < w->n_ex && res->n_ex < max_ex; i++) {
                memcpy(res->ex_cols + (size_t) res->n_ex * t,
                       w->ex_cols + (size_t) i * t, t * sizeof(int));
                res->ex_tuple[res->n_ex] = w->ex_tuple[i];
                res->n_ex++;
            }
        }
        R_CheckUserInterrupt();
    }
}

/* ================================================================== */
/* entry point                                                         */
/* ================================================================== */

/* kernel_: 0 = automatic, 1 = direct, 2 = prefix.
 * maxmem_: memory budget in bytes for the prefix kernel (automatic
 * choice only). */
SEXP C_ca_verify(SEXP mat, SEXP t_, SEXP v_, SEXP nthreads_, SEXP maxreport_,
                 SEXP kernel_, SEXP maxmem_)
{
    if (!Rf_isInteger(mat) || !Rf_isMatrix(mat))
        Rf_error("internal: matrix of integers expected");
    const int *m = INTEGER(mat);
    const int nrow = Rf_nrows(mat), ncol = Rf_ncols(mat);
    const int t = Rf_asInteger(t_);
    int nt = Rf_asInteger(nthreads_);
    int max_ex = Rf_asInteger(maxreport_);
    int kernel = Rf_asInteger(kernel_);
    const double maxmem = Rf_asReal(maxmem_);
    if (max_ex < 0 || max_ex == NA_INTEGER) max_ex = 0;

    if (t < 1 || t > MAX_T) Rf_error("t must be in 1..%d", MAX_T);
    if (ncol < t) Rf_error("array has fewer columns (%d) than t (%d)", ncol, t);
    if (!Rf_isInteger(v_) || Rf_length(v_) != ncol)
        Rf_error("internal: 'v' must be an integer vector of length ncol");
    const int *vs = INTEGER(v_);
    int vmax = 1;
    for (int j = 0; j < ncol; j++) {
        if (vs[j] == NA_INTEGER || vs[j] < 1)
            Rf_error("number of symbols for column %d must be >= 1", j + 1);
        if (vs[j] > vmax) vmax = vs[j];
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

    /* validate symbol range once, outside any parallel region */
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

    /* The prefix kernel needs vmax <= 64 (one word of value flags per
     * class and column in the gap path) and row indices below 2^32. */
    const int prefix_ok = vmax <= 64 && nrow < 2147483647;
    if (kernel == 2 && !prefix_ok)
        Rf_error("the prefix kernel needs at most 64 symbols per column");
    if (kernel == 0) {
        kernel = 1;
        if (prefix_ok && prefix_bytes(nrow, ncol, t, vmax, nt) <= maxmem)
            kernel = 2;
    }

    Result res;
    res.ex_cols = (int *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1) * t, sizeof(int));
    res.ex_tuple = (long long *) R_alloc((size_t)(max_ex > 0 ? max_ex : 1),
                                         sizeof(long long));
    if (kernel == 2)
        kernel_prefix(m, nrow, ncol, t, vs, vmax, nt, max_ex, &res);
    else
        kernel_direct(m, nrow, ncol, t, vs, vt_max, nt, max_ex, &res);

    /* examples: n_ex x (t cols 1-based, then t tuple symbols 0-based) */
    const int n_ex = res.n_ex;
    SEXP ex = PROTECT(Rf_allocMatrix(INTSXP, n_ex, 2 * t));
    int *e = INTEGER(ex);
    for (int row = 0; row < n_ex; row++) {
        const int *cc = res.ex_cols + (size_t) row * t;
        for (int q = 0; q < t; q++)
            e[row + (long long) q * n_ex] = cc[q] + 1;
        long long idx = res.ex_tuple[row];
        int digs[MAX_T];
        for (int q = t - 1; q >= 0; q--) {
            int vq = vs[cc[q]];
            digs[q] = (int)(idx % vq); idx /= vq;
        }
        for (int q = 0; q < t; q++)
            e[row + (long long)(t + q) * n_ex] = digs[q];
    }

    const char *names[] = {"covered", "colsets", "gaps", "missing_tuples",
                           "examples", "kernel", ""};
    SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(res.gaps == 0));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(res.combos));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(res.gaps));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(res.missing));
    SET_VECTOR_ELT(out, 4, ex);
    SET_VECTOR_ELT(out, 5, Rf_mkString(kernel == 2 ? "prefix" : "direct"));
    UNPROTECT(2);
    return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_ca_verify", (DL_FUNC) &C_ca_verify, 7},
    {NULL, NULL, 0}
};

void R_init_caverify(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
