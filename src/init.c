
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "islasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"inv_lapack", (DL_FUNC) &F77_SUB(inv_lapack), 6},
    {"standardize", (DL_FUNC) &F77_SUB(standardize), 6},
    {"check_out", (DL_FUNC) &F77_SUB(check_out), 6},

    {"family", (DL_FUNC) &F77_SUB(family), 6},

    {"binomial_variance", (DL_FUNC) &F77_SUB(binomial_variance), 3},
    {"logitlink", (DL_FUNC) &F77_SUB(logitlink), 3},
    {"logitlinkinv", (DL_FUNC) &F77_SUB(logitlinkinv), 3},
    {"logitmueta", (DL_FUNC) &F77_SUB(logitmueta), 3},
    {"probitlink", (DL_FUNC) &F77_SUB(probitlink), 3},
    {"probitlinkinv", (DL_FUNC) &F77_SUB(probitlinkinv), 3},
    {"probitmueta", (DL_FUNC) &F77_SUB(probitmueta), 3},

    {"poisson_variance", (DL_FUNC) &F77_SUB(poisson_variance), 3},
    {"loglink", (DL_FUNC) &F77_SUB(loglink), 3},
    {"loglinkinv", (DL_FUNC) &F77_SUB(loglinkinv), 3},
    {"logmueta", (DL_FUNC) &F77_SUB(logmueta), 3},

    {"gamma_variance", (DL_FUNC) &F77_SUB(gamma_variance), 3},
    {"inverselink", (DL_FUNC) &F77_SUB(inverselink), 3},
    {"inverselinkinv", (DL_FUNC) &F77_SUB(inverselinkinv), 3},
    {"inversemueta", (DL_FUNC) &F77_SUB(inversemueta), 3},
    {"identitylink", (DL_FUNC) &F77_SUB(identitylink), 3},
    {"identitylinkinv", (DL_FUNC) &F77_SUB(identitylinkinv), 3},
    {"identitymueta", (DL_FUNC) &F77_SUB(identitymueta), 3},

    {"penalty", (DL_FUNC) &F77_SUB(penalty), 5},
    {"gradient", (DL_FUNC) &F77_SUB(gradient), 10},
    {"hessian", (DL_FUNC) &F77_SUB(hessian), 8},
    {"hessian_theta", (DL_FUNC) &F77_SUB(hessian_theta), 8},

    {"islasso", (DL_FUNC) &F77_SUB(islasso), 37},
    {"islasso_glm", (DL_FUNC) &F77_SUB(islasso_glm), 39},

    {NULL, NULL, 0}
};

void attribute_visible R_init_islasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
