#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>

bool strata_ok(SEXP x);
int printInitialMessage();
int printMessage(int stepCounter, int numEvents, int currentPercent);
SEXP ComputeCox(SEXP _strata, SEXP _X, SEXP _w, SEXP _deltal, SEXP _WX, SEXP _n, SEXP _p, SEXP _numEvents);
