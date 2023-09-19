/* prokmeans.c    2023-08-15 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `prokmeans'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* declare functions here to register them below */

SEXP C_kmeans_dnorm(SEXP X, SEXP cluster0, SEXP PARA);
SEXP kmeansConstr_Call(SEXP x, SEXP CLS, SEXP PARS);

static R_CallMethodDef Call_entries[] = {
    {"C_kmeans_dnorm", (DL_FUNC) &C_kmeans_dnorm, 3},
    {"kmeansConstr_Call", (DL_FUNC) &kmeansConstr_Call, 3},
    {NULL, NULL, 0}
};

void R_init_prokmeans(DllInfo *info)
{
    R_registerRoutines(info, NULL, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
