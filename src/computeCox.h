#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>

bool strata_ok(SEXP x);
int printInitialMessage(int verbose);
int printMessage(int stepCounter, int num_events, int currentPercent, int verbose);
SEXP compute_cox_SEXP _strata, SEXP _X, SEXP _w,
                SEXP _deltal, SEXP _wx, SEXP _n,
                SEXP _p, SEXP _num_events, SEXP _verbose);
