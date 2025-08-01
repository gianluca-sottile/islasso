#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

void F77_SUB(inv_lapack)(int *n, double *a, double *ainv, int *info, int *ipiv, double *work);
void F77_SUB(standardize)(double *X, double *xm, double *xse, int *n, int *p, int *intercept);
void F77_SUB(check_out)(double *theta, double *cov, double *xm, double *xse, int *p, int *intercept);


void F77_SUB(family)(int *fam, int *link, int *func, double *x, int *n, double *y);

void F77_SUB(binomial_variance)(double *x, int *n, double *varmu);

void F77_SUB(logitlink)(double *x, int *n, double *mu);
void F77_SUB(logitlinkinv)(double *x, int *n, double *eta);
void F77_SUB(logitmueta)(double *x, int *n, double *eta);

void F77_SUB(probitlink)(double *x, int *n, double *mu);
void F77_SUB(probitlinkinv)(double *x, int *n, double *eta);
void F77_SUB(probitmueta)(double *x, int *n, double *eta);

void F77_SUB(poisson_variance)(double *x, int *n, double *varmu);

void F77_SUB(loglink)(double *x, int *n, double *mu);
void F77_SUB(loglinkinv)(double *x, int *n, double *eta);
void F77_SUB(logmueta)(double *x, int *n, double *eta);

void F77_SUB(gamma_variance)(double *x, int *n, double *varmu);

void F77_SUB(inverselink)(double *x, int *n, double *mu);
void F77_SUB(inverselinkinv)(double *x, int *n, double *eta);
void F77_SUB(inversemueta)(double *x, int *n, double *eta);

void F77_SUB(identitylink)(double *x, int *n, double *mu);
void F77_SUB(identitylinkinv)(double *x, int *n, double *eta);
void F77_SUB(identitymueta)(double *x, int *n, double *eta);


void F77_SUB(penalty)(double *theta, double *se, double *pi, int *p, double *pen, double *alpha);
void F77_SUB(gradient)(double *theta, double *se, double *lambda, double *xtw,
             double *res, double *pi, int *n, int *p, double *grad, double *alpha);
void F77_SUB(hessian_theta)(double *theta, double *se, double *lambda, double *xtx,
	     double *pi, int *p, double *hess, double *alpha);
void F77_SUB(hessian)(double *theta, double *se, double *lambda, double *xtx,
	     double *pi, int *p, double *hess, double *alpha);


void F77_SUB(islasso)(double *X, double *y, int *n, int *p, double *theta, double *se,
             double *cov, double *lambda, double *alpha, double *pi, int *estpi,
             int *itmax, int *itmaxse, double *tol, double *sigma2, int *trace,
             int *adaptive, double *offset, int *conv, int *stand, int *intercept,
             double *eta, double *mu, double *varmu, double *mu_eta_val, double *w,
             double *res, double *dev, double *weights, double *hi, double *edf,
             double *xtw, double *xtx, double *grad, double *hess, double *invH, double *pen);

void F77_SUB(islasso_glm)(double *X, double *y, int *n, int *p, double *theta, double *se,
             double *cov, double *lambda, double *alpha, double *pi, int *estpi,
             int *itmax, int *itmaxse, double *tol, double *sigma2, int *trace,
             int *adaptive, double *offset, int *conv, int *stand, int *intercept,
             double *eta, double *mu, double *varmu, double *mu_eta_val, double *w,
             double *res, double *dev, double *weights, double *hi, double *edf,
             double *xtw, double *xtx, double *grad, double *hess, double *invH, double *pen,
             int *fam, int *link);
