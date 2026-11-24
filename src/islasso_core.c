/*
 * islasso: numerical core of the induced smoothed lasso.
 *
 * At fixed smoothing scales s (standard errors) and mixture weights c the
 * IS estimating equation is the gradient of the smooth, strictly convex
 *
 *   F(theta) = 1/2 ||zt - Xt theta||^2 + sum_j lam_j { alpha G_j(theta_j) + (1 - alpha) theta_j^2 / 2 },
 *   G_j(b)   = c_j g(b; s_j) + (1 - c_j) g(b; eps s_j),   g(b; s) = E|b + s Z|,
 *
 * whose gradient contains P_j = c_j erf(t/sqrt2) + (1 - c_j) erf(t/(eps sqrt2)), t = b/s_j
 * (the smoothed sign of the IS-lasso) and whose Hessian H = Xt'Xt + diag(D) is the
 * slope matrix U'(theta) of the sandwich formula.  F is minimised by Newton's method
 * with an Armijo line search; at the solution the sandwich diagonal
 *   v = diag(H^-1 Xt'Xt H^-1)   and the pseudo hat values   h = diag(H^-1 Xt'Xt)
 * are returned.
 *
 * Two linear-algebra back ends:
 *   "dense": Xt'Xt (p x p) is given, Cholesky of H, cost O(p^3) per step;
 *   "wood" : Xt (n x p) is given, Woodbury/push-through identities, cost O(n^2 p) per step
 *            (used when p > n).  Unpenalised columns are handled exactly through a Schur
 *            complement.
 *
 * Compile with -DISL_STANDALONE for testing outside R.
 */

#include <math.h>
#include <string.h>
#include <stddef.h>

#ifdef ISL_STANDALONE
#include <stdlib.h>
#include <stdio.h>
#define FCONE , (size_t) 1
#define F77_CALL(x) x ## _
void dpotrf_(const char *, const int *, double *, const int *, int *, size_t);
void dpotrs_(const char *, const int *, const int *, const double *, const int *, double *,
             const int *, int *, size_t);
void dpotri_(const char *, const int *, double *, const int *, int *, size_t);
void dsyrk_(const char *, const char *, const int *, const int *, const double *, const double *,
            const int *, const double *, double *, const int *, size_t, size_t);
void dgemm_(const char *, const char *, const int *, const int *, const int *, const double *,
            const double *, const int *, const double *, const int *, const double *, double *,
            const int *, size_t, size_t);
void dgemv_(const char *, const int *, const int *, const double *, const double *, const int *,
            const double *, const int *, const double *, double *, const int *, size_t);
void dsymv_(const char *, const int *, const double *, const double *, const int *,
            const double *, const int *, const double *, double *, const int *, size_t);
static void *isl_alloc(size_t n, size_t sz) { void *p = calloc(n > 0 ? n : 1, sz); if (!p) abort(); return p; }
#define ALLOC_D(n) ((double *) isl_alloc((size_t) (n), sizeof(double)))
#define ALLOC_I(n) ((int *) isl_alloc((size_t) (n), sizeof(int)))
#define FREE(x) free(x)
#define CHECK_INTERRUPT()
#else
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#ifndef FCONE
#define FCONE
#endif
static double *isl_alloc_d(size_t n) { return (double *) R_alloc(n > 0 ? n : 1, sizeof(double)); }
static int *isl_alloc_i(size_t n) { return (int *) R_alloc(n > 0 ? n : 1, sizeof(int)); }
#define ALLOC_D(n) isl_alloc_d((size_t) (n))
#define ALLOC_I(n) isl_alloc_i((size_t) (n))
#define FREE(x) ((void) (x))
#define CHECK_INTERRUPT() R_CheckUserInterrupt()
#endif

#define ISL_EPS 1e-5            /* sd of the "zero mass" mixture component */
#define ISL_SQRT2 1.4142135623730950488
#define ISL_DNORM0 0.39894228040143267794
#define ISL_DELTA 1e-6          /* below this |theta| the smoothed sign is linearised */

/* ------------------------------------------------------------------ */
/* penalty pieces                                                      */
/* ------------------------------------------------------------------ */

/* g(b; s) = E|b + sZ| */
static double isl_g(double b, double s)
{
    double t = b / s;
    return b * erf(t / ISL_SQRT2) + 2.0 * s * ISL_DNORM0 * exp(-0.5 * t * t);
}

static double isl_P(double b, double s, double c)
{
    double t = b / s;
    return c * erf(t / ISL_SQRT2) + (1.0 - c) * erf(t / (ISL_EPS * ISL_SQRT2));
}

/*
 * Penalty value, gradient and curvatures.
 * For |theta_j| < ISL_DELTA the smoothed sign is linearised, P(b) ~ b P(delta)/delta, as in
 * the fixed-point iterations of islasso <= 1.6.2 (threshold 1e-6 on the standardised scale):
 * the objective stays convex and C^1.  Dn = curvature of this objective (Newton),
 * Ds = lam (alpha P'(theta) + 1 - alpha), the exact derivative used in the sandwich formula.
 */
static double isl_pen(int p, const double *th, const double *se, const double *c,
                      const double *lam, double alpha, double *gr, double *Dn, double *Ds)
{
    double val = 0.0;
    for (int j = 0; j < p; j++) {
        if (lam[j] <= 0.0) {
            if (gr) gr[j] = 0.0;
            if (Dn) Dn[j] = 0.0;
            if (Ds) Ds[j] = 0.0;
            continue;
        }
        double s = se[j], b = th[j], cj = c[j], ab = fabs(b);
        double lin = ab < ISL_DELTA;
        double Pd = lin ? isl_P(ISL_DELTA, s, cj) / ISL_DELTA : 0.0;
        if (gr) {
            double P = lin ? Pd * b : isl_P(b, s, cj);
            gr[j] = lam[j] * (alpha * P + (1.0 - alpha) * b);
        }
        if (Dn || Ds) {
            double t = b / s, t2 = t / ISL_EPS;
            double dP = 2.0 * ISL_DNORM0 * (cj * exp(-0.5 * t * t) +
                                            (1.0 - cj) * exp(-0.5 * t2 * t2) / ISL_EPS) / s;
            if (Ds) Ds[j] = lam[j] * (alpha * dP + (1.0 - alpha));
            if (Dn) Dn[j] = lam[j] * (alpha * (lin ? Pd : dP) + (1.0 - alpha));
        }
        double G;
        if (lin) {
            G = cj * isl_g(ISL_DELTA, s) + (1.0 - cj) * isl_g(ISL_DELTA, ISL_EPS * s) +
                0.5 * Pd * (b * b - ISL_DELTA * ISL_DELTA);
        } else {
            G = cj * isl_g(b, s) + (1.0 - cj) * isl_g(b, ISL_EPS * s);
        }
        val += lam[j] * (alpha * G + 0.5 * (1.0 - alpha) * b * b);
    }
    return val;
}

/* curvature of the MM (local quadratic) majoriser: lam (alpha P(b)/b + 1 - alpha), b floored
 * at ISL_DELTA (exactly the matrix of the fixed-point iterations of islasso <= 1.6.2) */
static void isl_mm_curv(int p, const double *th, const double *se, const double *c,
                        const double *lam, double alpha, double *D)
{
    for (int j = 0; j < p; j++) {
        if (lam[j] <= 0.0) { D[j] = 0.0; continue; }
        double ab = fabs(th[j]) < ISL_DELTA ? ISL_DELTA : fabs(th[j]);
        D[j] = lam[j] * (alpha * isl_P(ab, se[j], c[j]) / ab + (1.0 - alpha));
    }
}

/* ------------------------------------------------------------------ */
/* dense back end: G = Xt'Xt (p x p), bz = Xt'zt                        */
/* ------------------------------------------------------------------ */

typedef struct {
    int p;
    const double *G, *bz;
    double *L;      /* Cholesky factor of H (lower) */
    double *Gth;    /* G theta */
} dense_t;

/* F up to the constant ||zt||^2/2; also stores G theta */
static double dense_F(dense_t *ds, const double *th, double pen)
{
    int p = ds->p, one = 1;
    double d1 = 1.0, d0 = 0.0, q = 0.0;
    F77_CALL(dsymv)("L", &p, &d1, ds->G, &p, th, &one, &d0, ds->Gth, &one FCONE);
    for (int j = 0; j < p; j++) q += th[j] * (0.5 * ds->Gth[j] - ds->bz[j]);
    return q + pen;
}

/* factor H = G + diag(D); returns LAPACK info */
static int dense_factor(dense_t *ds, const double *D)
{
    int p = ds->p, info;
    size_t pp = (size_t) p * p;
    memcpy(ds->L, ds->G, pp * sizeof(double));
    for (int j = 0; j < p; j++) ds->L[j + (size_t) j * p] += D[j];
    F77_CALL(dpotrf)("L", &p, ds->L, &p, &info FCONE);
    return info;
}

static void dense_solve(dense_t *ds, double *x)
{
    int p = ds->p, one = 1, info;
    F77_CALL(dpotrs)("L", &p, &one, ds->L, &p, x, &p, &info FCONE);
}

/* sandwich diagonal v and hat values h from the factor of H (D = its penalty part).
 * Uses diag(H^-1 G H^-1) = diag(H^-1) - sum_k D_k (H^-1)_jk^2 and h = 1 - D diag(H^-1). */
static void dense_sandwich(dense_t *ds, const double *D, double *v, double *h)
{
    int p = ds->p, info;
    double *Hi = ds->L;
    F77_CALL(dpotri)("L", &p, Hi, &p, &info FCONE);
    for (int j = 0; j < p; j++) {           /* fill the upper triangle */
        for (int k = j + 1; k < p; k++) Hi[j + (size_t) k * p] = Hi[k + (size_t) j * p];
    }
    for (int j = 0; j < p; j++) {
        double s = 0.0, hjj = Hi[j + (size_t) j * p];
        const double *col = Hi + (size_t) j * p;
        for (int k = 0; k < p; k++) s += D[k] * col[k] * col[k];
        v[j] = hjj - s;
        if (v[j] < 0.0) v[j] = 0.0;
        h[j] = 1.0 - D[j] * hjj;
    }
}

/* ------------------------------------------------------------------ */
/* Woodbury back end: Xt (n x p) given                                 */
/* H = diag(D) + Xt'Xt, U = unpenalised columns, P = the others        */
/*   A = I + Xt_P D_P^-1 Xt_P',  S = D_U + Xt_U' A^-1 Xt_U              */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, p, nu, ucap;
    const double *X, *z, *cs;     /* cs = column sums of squares of Xt */
    const double *lam;
    int *U, *isU;                 /* indices of the exactly handled columns, membership */
    double *A;                    /* chol(A), n x n */
    double *XD;                   /* Xt scaled by D^-1/2 (U columns zero), n x p */
    double *XU, *AiXU, *S;        /* n x ucap, n x ucap, ucap x ucap */
    double *Dw;                   /* working curvature */
    double *r, *t, *tmp, *xt;     /* workspace: n, n, max(n,p), p */
    double *key; int *ord;        /* workspace for the selection of U */
} wood_t;

static int wood_init(wood_t *ws, int n, int p, const double *X, const double *z,
                     const double *lam)
{
    ws->n = n; ws->p = p; ws->X = X; ws->z = z; ws->lam = lam;
    double *cs = ALLOC_D(p);
    int nu0 = 0;
    for (int j = 0; j < p; j++) {
        double s = 0.0; const double *xj = X + (size_t) j * n;
        for (int i = 0; i < n; i++) s += xj[i] * xj[i];
        cs[j] = s > 0.0 ? s : 1.0;
        if (lam[j] <= 0.0) nu0++;
    }
    ws->cs = cs;
    ws->ucap = n < p ? n : p;
    if (nu0 > ws->ucap) ws->ucap = nu0;
    ws->U = ALLOC_I(ws->ucap); ws->isU = ALLOC_I(p);
    ws->A = ALLOC_D((size_t) n * n); ws->XD = ALLOC_D((size_t) n * p);
    ws->XU = ALLOC_D((size_t) n * ws->ucap); ws->AiXU = ALLOC_D((size_t) n * ws->ucap);
    ws->S = ALLOC_D((size_t) ws->ucap * ws->ucap); ws->Dw = ALLOC_D(p);
    ws->r = ALLOC_D(n); ws->t = ALLOC_D(n); ws->tmp = ALLOC_D(n > p ? n : p); ws->xt = ALLOC_D(p);
    ws->key = ALLOC_D(p); ws->ord = ALLOC_I(p);
    ws->nu = 0;
    return nu0;
}

static void wood_free(wood_t *ws)
{
    FREE((void *) ws->cs); FREE(ws->U); FREE(ws->isU); FREE(ws->A); FREE(ws->XD); FREE(ws->XU);
    FREE(ws->AiXU); FREE(ws->S); FREE(ws->Dw); FREE(ws->r); FREE(ws->t); FREE(ws->tmp);
    FREE(ws->xt); FREE(ws->key); FREE(ws->ord);
}

/* insertion of index j in the ordered (by key) list of length m, capacity cap */
static int ord_insert(int *ord, const double *key, int m, int cap, int j)
{
    int pos = m;
    while (pos > 0 && key[ord[pos - 1]] > key[j]) pos--;
    if (pos >= cap) return m;
    int last = m < cap ? m : cap - 1;
    for (int k = last; k > pos; k--) ord[k] = ord[k - 1];
    ord[pos] = j;
    return m < cap ? m + 1 : cap;
}

static double wood_F(wood_t *ws, const double *th, double pen)
{
    int n = ws->n, p = ws->p, one = 1;
    double d1 = 1.0, dm1 = -1.0, q = 0.0;
    memcpy(ws->r, ws->z, (size_t) n * sizeof(double));
    F77_CALL(dgemv)("N", &n, &p, &dm1, ws->X, &n, th, &one, &d1, ws->r, &one FCONE);
    for (int i = 0; i < n; i++) q += ws->r[i] * ws->r[i];
    return 0.5 * q + pen;
}

/* gradient of the quadratic part: -Xt'(zt - Xt th); requires ws->r from wood_F */
static void wood_lossgrad(wood_t *ws, double *g)
{
    int n = ws->n, p = ws->p, one = 1;
    double dm1 = -1.0, d0 = 0.0;
    F77_CALL(dgemv)("T", &n, &p, &dm1, ws->X, &n, ws->r, &one, &d0, g, &one FCONE);
}

/* Factorisation for H = diag(D) + Xt'Xt.  Columns with tiny curvature relative to
 * their scale (and the unpenalised ones) are handled exactly through the Schur
 * complement S; at most ucap of them, the others are floored at ISL_TAU_U * cs. */
#define ISL_TAU_U 1e-6
static int wood_factor(wood_t *ws, const double *D)
{
    int n = ws->n, p = ws->p, info = 0, m = 0, cap = ws->ucap;
    const double *lam = ws->lam;
    double d1 = 1.0, d0 = 0.0;
    for (int j = 0; j < p; j++) {
        ws->isU[j] = 0;
        ws->key[j] = lam[j] <= 0.0 ? -1.0 : D[j] / ws->cs[j];
        if (lam[j] <= 0.0 || !(D[j] > ISL_TAU_U * ws->cs[j]))
            m = ord_insert(ws->ord, ws->key, m, cap, j);
    }
    for (int k = 0; k < m; k++) { ws->U[k] = ws->ord[k]; ws->isU[ws->ord[k]] = 1; }
    ws->nu = m;
    int nu = m;
    for (int j = 0; j < p; j++) {
        double dj = D[j];
        if (lam[j] > 0.0 && !ws->isU[j] && !(dj > ISL_TAU_U * ws->cs[j])) dj = ISL_TAU_U * ws->cs[j];
        if (lam[j] > 0.0 && ws->isU[j]) {
            double fl = 1e-12 * ws->cs[j];
            if (!(dj > fl)) dj = fl;
        }
        ws->Dw[j] = dj;
        double sc = ws->isU[j] ? 0.0 : 1.0 / sqrt(dj);
        const double *xj = ws->X + (size_t) j * n;
        double *yj = ws->XD + (size_t) j * n;
        for (int i = 0; i < n; i++) yj[i] = sc * xj[i];
    }
    memset(ws->A, 0, (size_t) n * n * sizeof(double));
    for (int i = 0; i < n; i++) ws->A[i + (size_t) i * n] = 1.0;
    F77_CALL(dsyrk)("L", "N", &n, &p, &d1, ws->XD, &n, &d1, ws->A, &n FCONE FCONE);
    F77_CALL(dpotrf)("L", &n, ws->A, &n, &info FCONE);
    if (info) return info;
    if (nu > 0) {
        for (int k = 0; k < nu; k++)
            memcpy(ws->XU + (size_t) k * n, ws->X + (size_t) ws->U[k] * n, (size_t) n * sizeof(double));
        memcpy(ws->AiXU, ws->XU, (size_t) n * nu * sizeof(double));
        F77_CALL(dpotrs)("L", &n, &nu, ws->A, &n, ws->AiXU, &n, &info FCONE);
        F77_CALL(dgemm)("T", "N", &nu, &nu, &n, &d1, ws->XU, &n, ws->AiXU, &n, &d0, ws->S, &nu
                        FCONE FCONE);
        for (int k = 0; k < nu; k++) ws->S[k + (size_t) k * nu] += ws->Dw[ws->U[k]];
        F77_CALL(dpotrf)("L", &nu, ws->S, &nu, &info FCONE);
        if (info) return -info;
    }
    return 0;
}

/* x <- H^-1 x */
static void wood_solve(wood_t *ws, double *x)
{
    int n = ws->n, p = ws->p, nu = ws->nu, one = 1, info;
    double d1 = 1.0, d0 = 0.0, dm1 = -1.0;
    double *rp = ws->tmp, *t = ws->t, *xt = ws->xt;
    for (int j = 0; j < p; j++) rp[j] = ws->isU[j] ? 0.0 : x[j] / ws->Dw[j];
    /* t = A^-1 Xt rP */
    F77_CALL(dgemv)("N", &n, &p, &d1, ws->X, &n, rp, &one, &d0, t, &one FCONE);
    F77_CALL(dpotrs)("L", &n, &one, ws->A, &n, t, &n, &info FCONE);
    double *a = ws->key;   /* reuse (length p >= nu) */
    if (nu > 0) {
        for (int k = 0; k < nu; k++) a[k] = x[ws->U[k]];
        F77_CALL(dgemv)("T", &n, &nu, &dm1, ws->XU, &n, t, &one, &d1, a, &one FCONE);
        F77_CALL(dpotrs)("L", &nu, &one, ws->S, &nu, a, &nu, &info FCONE);
        /* t <- t + A^-1 XU a, so that bP = D^-1 (rP - Xt_P' t) */
        F77_CALL(dgemv)("N", &n, &nu, &d1, ws->AiXU, &n, a, &one, &d1, t, &one FCONE);
    }
    /* bP = rP - D^-1 Xt' t */
    F77_CALL(dgemv)("T", &n, &p, &d1, ws->X, &n, t, &one, &d0, xt, &one FCONE);
    for (int j = 0; j < p; j++) if (!ws->isU[j]) x[j] = rp[j] - xt[j] / ws->Dw[j];
    for (int k = 0; k < nu; k++) x[ws->U[k]] = a[k];
}

/* M = H^-1 Xt' (p x n, column major with leading dimension p) */
static void wood_HinvXt(wood_t *ws, double *M)
{
    int n = ws->n, p = ws->p, nu = ws->nu, info;
    double d1 = 1.0, d0 = 0.0, dm1 = -1.0;
    size_t nn = (size_t) n * n;
    double *Q = ALLOC_D(nn);
    memcpy(Q, ws->A, nn * sizeof(double));
    F77_CALL(dpotri)("L", &n, Q, &n, &info FCONE);
    for (int j = 0; j < n; j++) for (int k = j + 1; k < n; k++) Q[j + (size_t) k * n] = Q[k + (size_t) j * n];
    double *SiB = NULL;
    if (nu > 0) {
        /* SiB = S^-1 (A^-1 XU)'  (nu x n);  Q <- A^-1 - A^-1 XU SiB */
        SiB = ALLOC_D((size_t) nu * n);
        for (int i = 0; i < n; i++) for (int k = 0; k < nu; k++)
            SiB[k + (size_t) i * nu] = ws->AiXU[i + (size_t) k * n];
        F77_CALL(dpotrs)("L", &nu, &n, ws->S, &nu, SiB, &nu, &info FCONE);
        F77_CALL(dgemm)("N", "N", &n, &n, &nu, &dm1, ws->AiXU, &n, SiB, &nu, &d1, Q, &n FCONE FCONE);
    }
    /* M = Xt' Q, then rows of P scaled by 1/D, rows of U replaced by SiB */
    F77_CALL(dgemm)("T", "N", &p, &n, &n, &d1, ws->X, &n, Q, &n, &d0, M, &p FCONE FCONE);
    for (int j = 0; j < p; j++) {
        if (!ws->isU[j]) {
            double sc = 1.0 / ws->Dw[j];
            for (int i = 0; i < n; i++) M[j + (size_t) i * p] *= sc;
        }
    }
    for (int k = 0; k < nu; k++) {
        int j = ws->U[k];
        for (int i = 0; i < n; i++) M[j + (size_t) i * p] = SiB[k + (size_t) i * nu];
    }
    FREE(Q);
    if (SiB) { FREE(SiB); }
}

static void wood_sandwich(wood_t *ws, double *v, double *h, double *M)
{
    int n = ws->n, p = ws->p;
    wood_HinvXt(ws, M);
    for (int j = 0; j < p; j++) { v[j] = 0.0; h[j] = 0.0; }
    for (int i = 0; i < n; i++) {
        const double *mi = M + (size_t) i * p;
        for (int j = 0; j < p; j++) {
            double m = mi[j];
            v[j] += m * m;
            h[j] += m * ws->X[i + (size_t) j * n];
        }
    }
}

/* ------------------------------------------------------------------ */
/* inner solver                                                        */
/* ------------------------------------------------------------------ */

/*
 * Minimises F at fixed (se, c).
 * mode 0 = dense (G, bz), 1 = wood (X, z).
 * method 0 = MM (local quadratic majoriser, monotone; the undamped version of the
 *            fixed-point iteration of islasso <= 1.6.2), 1 = Newton with Armijo search.
 * On exit: theta updated; v, h = sandwich diagonal and hat values at the solution;
 * D = exact curvature (penalty part of U'). Returns the number of steps;
 * *info: 0 ok, 1 = max iterations, 2 = factorisation failure, 3 = line search stalled.
 */
int isl_inner(int mode, int method, int n, int p, const double *G, const double *bz,
              const double *X, const double *z, double *theta, const double *se,
              const double *c, const double *lam, double alpha, double tol, int maxit,
              double *v, double *h, double *D, int *info)
{
    double *g = ALLOC_D(p), *d = ALLOC_D(p), *thn = ALLOC_D(p), *Ds = ALLOC_D(p);
    double *gl = ALLOC_D(p);
    dense_t ds; wood_t ws;
    int it = 0, fact;
    *info = 0;
    if (mode == 0) {
        ds.p = p; ds.G = G; ds.bz = bz;
        ds.L = ALLOC_D((size_t) p * p); ds.Gth = ALLOC_D(p);
    } else {
        wood_init(&ws, n, p, X, z, lam);
    }

    double *Dn = ALLOC_D(p);
    double pen = isl_pen(p, theta, se, c, lam, alpha, g, Dn, NULL);
    double f = mode == 0 ? dense_F(&ds, theta, pen) : wood_F(&ws, theta, pen);

    int usemm = (method == 0);
    for (it = 1; it <= maxit; it++) {
        CHECK_INTERRUPT();
        if (mode == 0) {
            for (int j = 0; j < p; j++) g[j] += ds.Gth[j] - bz[j];
        } else {
            wood_lossgrad(&ws, gl);
            for (int j = 0; j < p; j++) g[j] += gl[j];
        }
        if (usemm) isl_mm_curv(p, theta, se, c, lam, alpha, Ds);
        else memcpy(Ds, Dn, (size_t) p * sizeof(double));
        fact = mode == 0 ? dense_factor(&ds, Ds) : wood_factor(&ws, Ds);
        if (fact && !usemm) {                 /* Newton: fall back to the MM majoriser */
            usemm = 1;
            isl_mm_curv(p, theta, se, c, lam, alpha, Ds);
            fact = mode == 0 ? dense_factor(&ds, Ds) : wood_factor(&ws, Ds);
        }
        if (fact) { *info = 2; break; }
        double gd = 0.0, dmax = 0.0;
        for (int j = 0; j < p; j++) d[j] = -g[j];
        if (mode == 0) dense_solve(&ds, d); else wood_solve(&ws, d);
        for (int j = 0; j < p; j++) gd += g[j] * d[j];
        if (!(gd < 0.0)) {
            if (gd == 0.0) break;             /* stationary point */
            gd = 0.0;                          /* not a descent direction: steepest descent */
            for (int j = 0; j < p; j++) { d[j] = -g[j]; gd -= g[j] * g[j]; }
        }
        /* MM: full step decreases F in exact arithmetic; backtracking as a safeguard */
        double step = 1.0, fn = f, sig = usemm ? 0.0 : 1e-4;
        int ok = 0;
        for (int ls = 0; ls < 60; ls++) {
            for (int j = 0; j < p; j++) thn[j] = theta[j] + step * d[j];
            double pn = isl_pen(p, thn, se, c, lam, alpha, NULL, NULL, NULL);
            fn = mode == 0 ? dense_F(&ds, thn, pn) : wood_F(&ws, thn, pn);
            if (isfinite(fn) && fn <= f + sig * step * gd + 1e-14 * fabs(f)) { ok = 1; break; }
            step *= 0.5;
        }
        if (!ok) {                             /* restore cached quantities at theta */
            pen = isl_pen(p, theta, se, c, lam, alpha, g, Dn, NULL);
            f = mode == 0 ? dense_F(&ds, theta, pen) : wood_F(&ws, theta, pen);
            if (!usemm) { usemm = 1; continue; }   /* Newton stalled: switch to MM steps */
            *info = 3;
            break;
        }
        /* convergence is judged on the full direction (not on the damped step) */
        for (int j = 0; j < p; j++) {
            double dj = fabs(d[j]);
            if (dj > dmax) dmax = dj;
            theta[j] = thn[j];
        }
        f = fn;
        pen = isl_pen(p, theta, se, c, lam, alpha, g, Dn, NULL);
        if (dmax <= tol) break;
    }
    if (it > maxit) { it = maxit; *info = 1; }

    /* sandwich at the solution (same se): exact curvature D */
    isl_pen(p, theta, se, c, lam, alpha, NULL, NULL, D);
    if (mode == 0) {
        if (dense_factor(&ds, D)) *info = 2;
        else dense_sandwich(&ds, D, v, h);
        FREE(ds.L); FREE(ds.Gth);
    } else {
        if (wood_factor(&ws, D)) *info = 2;
        else {
            double *M = ALLOC_D((size_t) p * n);
            wood_sandwich(&ws, v, h, M);
            FREE(M);
        }
        wood_free(&ws);
    }
    FREE(g); FREE(d); FREE(thn); FREE(Ds); FREE(gl); FREE(Dn);
    return it;
}

/*
 * Full unscaled sandwich V = H^-1 Xt'Xt H^-1 (p x p) and, optionally, the observation
 * leverages hobs_i = x_i' H^-1 x_i (X required).  mode 0 needs G; mode 1 needs X.
 * Returns 0 on success.
 */
int isl_vcov(int mode, int n, int p, const double *G, const double *X, const double *theta,
             const double *se, const double *c, const double *lam, double alpha,
             double *V, double *hobs)
{
    double *D = ALLOC_D(p);
    double d1 = 1.0, d0 = 0.0;
    int info = 0;
    isl_pen(p, theta, se, c, lam, alpha, NULL, NULL, D);
    if (mode == 0) {
        dense_t ds; ds.p = p; ds.G = G; ds.L = ALLOC_D((size_t) p * p);
        if (dense_factor(&ds, D)) { FREE(ds.L); FREE(D); return 2; }
        double *Hi = ds.L;
        F77_CALL(dpotri)("L", &p, Hi, &p, &info FCONE);
        for (int j = 0; j < p; j++) for (int k = j + 1; k < p; k++) Hi[j + (size_t) k * p] = Hi[k + (size_t) j * p];
        double *B = ALLOC_D((size_t) p * p);
        F77_CALL(dgemm)("N", "N", &p, &p, &p, &d1, Hi, &p, G, &p, &d0, B, &p FCONE FCONE);
        F77_CALL(dgemm)("N", "N", &p, &p, &p, &d1, B, &p, Hi, &p, &d0, V, &p FCONE FCONE);
        for (int j = 0; j < p; j++) for (int k = j + 1; k < p; k++) {
            double m = 0.5 * (V[j + (size_t) k * p] + V[k + (size_t) j * p]);
            V[j + (size_t) k * p] = V[k + (size_t) j * p] = m;
        }
        if (hobs && X) {
            double *XH = ALLOC_D((size_t) n * p);
            F77_CALL(dgemm)("N", "N", &n, &p, &p, &d1, X, &n, Hi, &p, &d0, XH, &n FCONE FCONE);
            for (int i = 0; i < n; i++) hobs[i] = 0.0;
            for (int j = 0; j < p; j++) for (int i = 0; i < n; i++)
                hobs[i] += X[i + (size_t) j * n] * XH[i + (size_t) j * n];
            FREE(XH);
        }
        FREE(B); FREE(ds.L);
    } else {
        wood_t ws;
        wood_init(&ws, n, p, X, NULL, lam);
        if (wood_factor(&ws, D)) info = 2;
        else {
            double *M = ALLOC_D((size_t) p * n);
            wood_HinvXt(&ws, M);
            F77_CALL(dsyrk)("L", "N", &p, &n, &d1, M, &p, &d0, V, &p FCONE FCONE);
            for (int j = 0; j < p; j++) for (int k = j + 1; k < p; k++) V[j + (size_t) k * p] = V[k + (size_t) j * p];
            if (hobs) {
                for (int i = 0; i < n; i++) {
                    double s = 0.0;
                    for (int j = 0; j < p; j++) s += X[i + (size_t) j * n] * M[j + (size_t) i * p];
                    hobs[i] = s;
                }
            }
            FREE(M);
        }
        wood_free(&ws);
    }
    FREE(D);
    return info;
}

/* quadratic form b' H^-1 G H^-1 b of the unscaled sandwich (e.g. variance of a linear
 * combination), without forming the p x p matrix.  Returns NAN on failure. */
double isl_quad(int mode, int n, int p, const double *G, const double *X, const double *theta,
                const double *se, const double *c, const double *lam, double alpha,
                const double *b)
{
    double *D = ALLOC_D(p), *x = ALLOC_D(p), q = 0.0, d1 = 1.0, d0 = 0.0;
    int one = 1;
    isl_pen(p, theta, se, c, lam, alpha, NULL, NULL, D);
    memcpy(x, b, (size_t) p * sizeof(double));
    if (mode == 0) {
        dense_t ds; ds.p = p; ds.G = G; ds.L = ALLOC_D((size_t) p * p); ds.Gth = ALLOC_D(p);
        if (dense_factor(&ds, D)) q = NAN;
        else {
            dense_solve(&ds, x);
            F77_CALL(dsymv)("L", &p, &d1, G, &p, x, &one, &d0, ds.Gth, &one FCONE);
            for (int j = 0; j < p; j++) q += x[j] * ds.Gth[j];
        }
        FREE(ds.L); FREE(ds.Gth);
    } else {
        wood_t ws;
        wood_init(&ws, n, p, X, NULL, lam);
        if (wood_factor(&ws, D)) q = NAN;
        else {
            wood_solve(&ws, x);
            F77_CALL(dgemv)("N", &n, &p, &d1, X, &n, x, &one, &d0, ws.r, &one FCONE);
            for (int i = 0; i < n; i++) q += ws.r[i] * ws.r[i];
        }
        wood_free(&ws);
    }
    FREE(D); FREE(x);
    return q;
}

/* penalty, gradient and curvature at given theta (exported for testing and for R) */
double isl_penalty(int p, const double *th, const double *se, const double *c,
                   const double *lam, double alpha, double *gr, double *D)
{
    return isl_pen(p, th, se, c, lam, alpha, gr, NULL, D);
}

#ifndef ISL_STANDALONE
/* ------------------------------------------------------------------ */
/* .Call interface                                                     */
/* ------------------------------------------------------------------ */

static int get_int(SEXP x) { return asInteger(x); }

SEXP C_is_inner(SEXP mode_, SEXP method_, SEXP G_, SEXP bz_, SEXP X_, SEXP z_, SEXP theta_,
                SEXP se_, SEXP c_, SEXP lam_, SEXP alpha_, SEXP tol_, SEXP maxit_)
{
    int mode = get_int(mode_), method = get_int(method_), p = length(theta_), n = 0;
    const double *G = NULL, *bz = NULL, *X = NULL, *z = NULL;
    if (mode == 0) {
        if (!isReal(G_) || !isReal(bz_) || XLENGTH(G_) != (R_xlen_t) p * p || length(bz_) != p)
            error("invalid dense input");
        G = REAL(G_); bz = REAL(bz_);
    } else {
        if (!isReal(X_) || !isReal(z_)) error("invalid input");
        n = length(z_);
        if (XLENGTH(X_) != (R_xlen_t) n * p) error("dimension mismatch");
        X = REAL(X_); z = REAL(z_);
    }
    if (p < 1) error("empty coefficient vector");
    if (!isReal(theta_) || !isReal(se_) || !isReal(c_) || !isReal(lam_) ||
        length(se_) != p || length(c_) != p || length(lam_) != p) error("invalid input");
    SEXP out = PROTECT(allocVector(VECSXP, 6));
    SEXP th = PROTECT(duplicate(theta_));
    SEXP v = PROTECT(allocVector(REALSXP, p));
    SEXP h = PROTECT(allocVector(REALSXP, p));
    SEXP D = PROTECT(allocVector(REALSXP, p));
    int info = 0;
    int it = isl_inner(mode, method, n, p, G, bz, X, z, REAL(th), REAL(se_), REAL(c_), REAL(lam_),
                        asReal(alpha_), asReal(tol_), get_int(maxit_), REAL(v), REAL(h),
                        REAL(D), &info);
    SET_VECTOR_ELT(out, 0, th);
    SET_VECTOR_ELT(out, 1, v);
    SET_VECTOR_ELT(out, 2, h);
    SET_VECTOR_ELT(out, 3, D);
    SET_VECTOR_ELT(out, 4, ScalarInteger(it));
    SET_VECTOR_ELT(out, 5, ScalarInteger(info));
    SEXP nms = PROTECT(allocVector(STRSXP, 6));
    const char *nn[6] = {"theta", "v", "h", "D", "iter", "info"};
    for (int k = 0; k < 6; k++) SET_STRING_ELT(nms, k, mkChar(nn[k]));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
}

SEXP C_is_vcov(SEXP mode_, SEXP G_, SEXP X_, SEXP theta_, SEXP se_, SEXP c_, SEXP lam_,
               SEXP alpha_, SEXP hobs_)
{
    int mode = get_int(mode_), p = length(theta_), n = 0, want_h = asLogical(hobs_);
    const double *G = NULL, *X = NULL;
    if (p < 1) error("empty coefficient vector");
    if (!isReal(theta_) || !isReal(se_) || !isReal(c_) || !isReal(lam_) ||
        length(se_) != p || length(c_) != p || length(lam_) != p) error("invalid input");
    if (mode == 0) {
        if (!isReal(G_) || XLENGTH(G_) != (R_xlen_t) p * p) error("invalid dense input");
        G = REAL(G_);
    }
    if (isReal(X_)) {
        X = REAL(X_);
        n = (int) (XLENGTH(X_) / p);
    } else if (mode == 1 || want_h) error("'X' is required");
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP V = PROTECT(allocMatrix(REALSXP, p, p));
    SEXP hb = PROTECT(allocVector(REALSXP, want_h ? n : 0));
    int info = isl_vcov(mode, n, p, G, X, REAL(theta_), REAL(se_), REAL(c_), REAL(lam_),
                        asReal(alpha_), REAL(V), want_h ? REAL(hb) : NULL);
    SET_VECTOR_ELT(out, 0, V);
    SET_VECTOR_ELT(out, 1, hb);
    SET_VECTOR_ELT(out, 2, ScalarInteger(info));
    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, mkChar("V"));
    SET_STRING_ELT(nms, 1, mkChar("hobs"));
    SET_STRING_ELT(nms, 2, mkChar("info"));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
}

SEXP C_is_quad(SEXP mode_, SEXP G_, SEXP X_, SEXP theta_, SEXP se_, SEXP c_, SEXP lam_,
               SEXP alpha_, SEXP b_)
{
    int mode = get_int(mode_), p = length(theta_), n = 0;
    if (p < 1 || length(b_) != p || !isReal(b_)) error("invalid input");
    if (!isReal(theta_) || !isReal(se_) || !isReal(c_) || !isReal(lam_) ||
        length(se_) != p || length(c_) != p || length(lam_) != p) error("invalid input");
    const double *G = NULL, *X = NULL;
    if (mode == 0) {
        if (!isReal(G_) || XLENGTH(G_) != (R_xlen_t) p * p) error("invalid dense input");
        G = REAL(G_);
    } else {
        if (!isReal(X_)) error("'X' is required");
        X = REAL(X_);
        n = (int) (XLENGTH(X_) / p);
    }
    return ScalarReal(isl_quad(mode, n, p, G, X, REAL(theta_), REAL(se_), REAL(c_), REAL(lam_),
                               asReal(alpha_), REAL(b_)));
}

SEXP C_is_penalty(SEXP theta_, SEXP se_, SEXP c_, SEXP lam_, SEXP alpha_)
{
    int p = length(theta_);
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP gr = PROTECT(allocVector(REALSXP, p));
    SEXP D = PROTECT(allocVector(REALSXP, p));
    double val = isl_pen(p, REAL(theta_), REAL(se_), REAL(c_), REAL(lam_), asReal(alpha_),
                         REAL(gr), NULL, REAL(D));
    SET_VECTOR_ELT(out, 0, ScalarReal(val));
    SET_VECTOR_ELT(out, 1, gr);
    SET_VECTOR_ELT(out, 2, D);
    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, mkChar("value"));
    SET_STRING_ELT(nms, 1, mkChar("gradient"));
    SET_STRING_ELT(nms, 2, mkChar("curvature"));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
}
#endif
