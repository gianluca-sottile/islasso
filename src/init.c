#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP C_is_inner(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_is_vcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_is_penalty(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_is_quad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_is_gauss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"is_inner", (DL_FUNC) &C_is_inner, 13},
    {"is_vcov", (DL_FUNC) &C_is_vcov, 9},
    {"is_penalty", (DL_FUNC) &C_is_penalty, 5},
    {"is_quad", (DL_FUNC) &C_is_quad, 9},
    {"is_gauss", (DL_FUNC) &C_is_gauss, 18},
    {NULL, NULL, 0}
};

void attribute_visible R_init_islasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
