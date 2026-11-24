#include <R.h>

void F77_SUB(islasso_trace2_3)(int *i) {
    if(*i == 1) Rprintf("\nFisher scoring step:\n");
    Rprintf(".");
    if((*i % 100) == 0) Rprintf(" %5d\n", *i);
}

void F77_SUB(islasso_trace1_2_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     sup|BETAn - BETAo|1 = %2.8f\n", *ind2);
    Rprintf("     sup|SEn - SEo| = %2.8f\n", *ind);
}

void F77_SUB(islasso_trace1_7_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    if(*i == 1)Rprintf("\nIS-lasso algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, DEV = %10.4f (%5.2f df), phi = %7.4f, sup|BETAn - BETAo| = %2.8f e sup|SEn - SEo| = %2.8f\n", *i, *dev, *redf, *s2, *ind2, *ind);
}

void F77_SUB(islasso_trace1_8)(int *i) {
    Rprintf("\n===================================");
    if(*i == 1)Rprintf("\nConvergence criterion is met!\n\n");
}

void F77_SUB(islasso_trace2_2_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    Rprintf("\n=============================================================\n");
    Rprintf("IS-lasso (GLM) algorithm step = %d\n", *i);
    Rprintf("  Choosen lambda value = %7.3f\n", *lmb);
    Rprintf("  Residual deviance = %10.6f on %5.2f degrees of freedom\n", *dev, *redf);
    Rprintf("  Estimated dispersion parameter = %7.4f\n", *s2);
    Rprintf("  Checking convergence criterion (threshold = %g):\n", *eps);
    Rprintf("     sup|BETAn - BETAo| = %2.8f\n", *ind2);
    Rprintf("     sup|SEn - SEo| = %2.8f\n", *ind);
}

void F77_SUB(islasso_trace2_7_2)(double *eps, int *i, double *lmb, double *dev, double *redf, double *s2, double *ind, double *ind2) {
    if(*i == 1)Rprintf("\nIS-lasso (GLM) algorithm (choosen lambda = %7.3f, threshold = %g)\n\n", *lmb, *eps);
    Rprintf("Step = %4d, DEV = %10.4f (%5.2f df), phi = %7.4f, sup|BETAn - BETAo| = %2.8f e sup|SEn - SEo| = %2.8f\n", *i, *dev, *redf, *s2, *ind2, *ind);
}

void F77_SUB(islasso_trace2_6)(int *i) {
    Rprintf(" %5d. Convergence criterion is met!\n\n", *i);
}
