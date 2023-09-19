/* kmeans_dnorm.c    2023-08-16 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `prokmeans'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

#define INDEX(i, j) i + j * n

int kmeansConstr(double *x, int n, int p, int *cls0, int *cls,
		 int K, int iterlim, int restartlim, int diag,
		 int verbose)
{
    int c, i, j, k, *N, *N0, pK = p * K, Nreclass, Nrclprv, iter, restart = 0, Ndiag, propc;
    double *sums, *sums0, *means, d, dmin;

    /* assume no NA's, so the counts are the same for all variables */

    N0 = (int*)R_alloc(K, sizeof(int));
    N = (int*)R_alloc(K, sizeof(int));

    sums0 = (double*)R_alloc(pK, sizeof(double));
    sums = (double*)R_alloc(pK, sizeof(double));
    means = (double*)R_alloc(pK, sizeof(double));
    /* in these arrays, the rows are the p variables (columns) of x,
       and the columns are the K groups */

    memset(N0, 0, K * sizeof(int));
    memset(sums0, 0, pK * sizeof(double));

    /* compute the sums and counts only for the observations that
       are constrained and cannot be reassigned */
    for (i = 0; i < n; i++) {
	if (!cls0[i]) continue;
	c = cls[i] - 1;
	(N0[c])++;
	for (j = 0, k = p * c; j < p; j++, k++)
	    sums0[k] += x[INDEX(i, j)];
    }

 start:
    for (i = 0; i < n; i++) {
	if (cls0[i]) continue;
	GetRNGstate();
	cls[i] = (int)ceil(unif_rand() * K);
	PutRNGstate();
    }
    Ndiag = 0;

    if (verbose) Rprintf("Restart %d:\n", restart);

    iter = 1;
    Nrclprv = n;
    Nreclass = 0; /* avoids a warning during compilation */
    while (iter < iterlim) {
	/* initialise the sums and counts from the constrained observations */
	memcpy(sums, sums0, pK * sizeof(double));
	memcpy(N, N0, K * sizeof(int));

	/* do the sums and counts */
	for (i = 0; i < n; i++) {
	    if (cls0[i]) continue;
	    c = cls[i] - 1;
	    (N[c])++;
	    for (j = 0, k = p * c; j < p; j++, k++)
		sums[k] += x[INDEX(i, j)];
	}
	/* compute the means */
	for (c = 0; c < K; c++) {
	    if (!N[c]) continue;
	    for (j = 0, k = p * c; j < p; j++, k++)
		means[k] = sums[k] / N[c];
	}

	for (i = 0; i < n; i++) {
	    if (cls0[i]) continue;
	    dmin = R_PosInf;
	    propc = -1;
	    for (c = 0; c < K; c++) {
		if (!N[c]) continue; /* skip the empty classes */
		d = 0;
		for (j = 0, k = p * c; j < p; j++, k++)
		    d += pow(x[INDEX(i, j)] - means[k], 2);
		/* this is the Euclidean distance, but no need to
		   compute the square root because we compare sums of
		   squared differences */
		if (d < dmin) {
		    dmin = d;
		    propc = c;
		}
	    }
	    if (propc == -1) continue;
	    propc++;
	    if (propc != cls[i]) {
		cls[i] = propc;
		Nreclass++;
	    }
	}
	if (verbose)
	    Rprintf("  iteration %d -> %d reclassified\n", iter, Nreclass);
	if (!Nreclass) break;
	if (Nreclass >= Nrclprv) Ndiag++; else Ndiag = 0;
	if (Ndiag >= diag) {
	    restart++;
	    if (restart > restartlim) {
		if (verbose)
		    Rprintf("No convergence after %d restarts; results likely dubious.\n", restart);
		break;
	    }
	    goto start;
	}
	Nrclprv = Nreclass;
	Nreclass = 0; /* avoids a warning during compilation */
	iter++;
    }
    return Nreclass;
}

SEXP kmeansConstr_Call(SEXP x, SEXP CLS, SEXP PARS)
{
    int *cls0, *cls, n, p, *pars, v;
    SEXP res;

    PROTECT(x = coerceVector(x, REALSXP));
    PROTECT(CLS = coerceVector(CLS, INTSXP));
    PROTECT(PARS = coerceVector(PARS, INTSXP));
    n = nrows(x);
    p = ncols(x);
    cls0 = INTEGER(CLS);
    pars = INTEGER(PARS);

    /* the values in CLS are not erased by copying into res
       (cls points to res, not to CLS) */
    PROTECT(res = allocVector(INTSXP, n));
    cls = INTEGER(res);
    memcpy(cls, cls0, n * sizeof(int));

    /* both cls0 and cls are passed to kmeansConstr() */
    v = kmeansConstr(REAL(x), n, p, cls0, cls, pars[0], pars[1], pars[2], pars[3], pars[4]);
    if (v && pars[4]) Rprintf("Reached iteration number limit.\n");

    UNPROTECT(4);
    return res;
}
