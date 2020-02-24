/*
 * SEXP method2,
 * SEXP eps2,
 * method  = asInteger(method2);
 * eps     = asReal(eps2);
 */

#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include "computeCox.h"

bool strata_ok(SEXP x)
{
  bool ok = true;
  for (int i = 0; i < length(x); i++) {
    SEXP name = getAttrib(VECTOR_ELT(x, i), R_NamesSymbol);
    ok = ok && length(name) == 8 &&
      strcmp("start",  CHAR(STRING_ELT(name, 0))) == 0 &&
      strcmp("end",    CHAR(STRING_ELT(name, 1))) == 0 &&
      strcmp("label",  CHAR(STRING_ELT(name, 2))) == 0 &&
      strcmp("J",      CHAR(STRING_ELT(name, 3))) == 0 &&
      strcmp("nfails", CHAR(STRING_ELT(name, 4))) == 0 &&
      strcmp("start0", CHAR(STRING_ELT(name, 5))) == 0 &&
      strcmp("start1", CHAR(STRING_ELT(name, 6))) == 0 &&
      strcmp("stop1",  CHAR(STRING_ELT(name, 7))) == 0;
  }
  return(ok);
}

int printInitialMessage() {
  Rprintf("Processing W*X               :   0%%|                    |\r");
  fflush(stdout);
  return(0);
}

int printMessage(int stepCounter, int numEvents, int currentPercent)
{
  int newPercent = 100 * stepCounter / (float)numEvents;
   int stars      = 20 * stepCounter  / (float)numEvents;
  if (newPercent > currentPercent) {
    Rprintf("Processing W*X               : %3d%%|", newPercent);
    for (int i = 0; i < stars; i++ ) { Rprintf("#"); }
    for (int i = stars; i < 20; i++) { Rprintf(" "); }
    Rprintf("|\r");
    if (stepCounter == numEvents) {
      Rprintf("\n\n");
    }
  }
  fflush(stdout);
  return(newPercent);
}


SEXP ComputeCox(SEXP _strata, SEXP _X, SEXP _w, SEXP _deltal, SEXP _WX, SEXP _n, SEXP _p, SEXP _numEvents)
{
  bool ok = strata_ok(_strata);
  if (!ok) {
    Rprintf("ERROR IN STRATA\n");
    SEXP ans;
    ans = allocList(0);
    return(ans);
  }
  int n          = INTEGER(_n)[0];  // use asInteger(_n);
  int p          = INTEGER(_p)[0];  // use asInteger(_p);
  int numEvents  = INTEGER(_numEvents)[0]; // use asInteger(_numEvents);
  double *X      = REAL(_X);
  double *w      = REAL(_w);
  double *deltal = REAL(_deltal);
  double *WX     = REAL(_WX);

  double *wz    = (double*)malloc(sizeof(double) * n);
  double *YWX   = (double*)malloc(sizeof(double) * p);
  double *ZWX   = (double*)malloc(sizeof(double) * p);

  if (wz == NULL || YWX == NULL || ZWX == NULL) {
    Rprintf("Error allocating memory\n");
    free(wz);
    free(YWX);
    free(ZWX);
    SEXP ans;
    ans = allocList(0);
    return(ans);
  }

  int currentPercent = printInitialMessage();
  int stepCounter = 0;
  //Rprintf("n = %d\np = %d\n", n, p);
  //Rprintf("length(strata) = %d\n", length(_strata));
  for (int i = 0; i < length(_strata); i++) {
    //Rprintf("i = %d\n", i);
    SEXP level   = VECTOR_ELT(_strata, i);
    int end     = INTEGER(VECTOR_ELT(level, 1))[0];
    int J       = INTEGER(VECTOR_ELT(level, 3))[0];
    double *nfails = REAL(VECTOR_ELT(level, 4));
    double *start0 = REAL(VECTOR_ELT(level, 5));
    double *start1 = REAL(VECTOR_ELT(level, 6));
    double *stop1  = REAL(VECTOR_ELT(level, 7));
    //Rprintf("J = %d\n", J);
    if (J > 0) {
      //Rprintf("J > 0\n");
      for (int j = 0; j < J; j++) {
        //Rprintf("%d / %d\n", j, J);
        int nj = nfails[j];
        int yIndexStart = start0[j] - 1;
        int yIndexEnd   = end;
        int zIndexStart = start1[j] - 1;
        int zIndexEnd   = stop1[j];
        //Rprintf("%d %d %d %d %d\n", nj, yIndexStart, yIndexEnd, zIndexStart, zIndexEnd);
        double Aj1 = 0.;
        double Aj2 = 0.;
        for (int yindex = yIndexStart; yindex < yIndexEnd; yindex++) {
          if (zIndexStart <= yindex && yindex < zIndexEnd) {
            wz[yindex] = w[yindex] / nj;
            Aj2 += wz[yindex];
          } else {
            wz[yindex] = 0.;
          }
          Aj1 += w[yindex];
        }
        double h1, h2, h3, h4, h5;
        h1 = h2 = h3 = h4 = h5 = 0.;
        for (int r = 0; r < nj; r++) {
          double Ajr = Aj1 - r * Aj2;
          double denom2 = 1 / Ajr;
          double denom1 = denom2 / Ajr;
          h1 += denom1;
          h2 += r * denom1;
          h3 += r * r * denom1;
          h4 += denom2;
          h5 += r * denom2;
        }
        for (int pct = 0; pct < p; pct++) {
          YWX[pct] = ZWX[pct] = 0.;
        }
        for (long row = yIndexStart; row < yIndexEnd; row++) {
          for (long col = 0; col < p; col++) {
            long idx = col * n + row;
            YWX[col] += w[row]  * X[idx];
            ZWX[col] += wz[row] * X[idx];
          }
          deltal[row] += h5 * wz[row] - h4 * w[row];
        }
        for (long col = 0; col < p; col++) {
          for (long row = yIndexStart; row < yIndexEnd; row++) {
            long idx = col * n + row;
            WX[idx] +=
              (h2 * wz[row] - h1 *  w[row]) * YWX[col] +
              (h2 *  w[row] - h3 * wz[row]) * ZWX[col] +
              (h4 *  w[row] - h5 * wz[row]) * X[idx];
          }
        }
        stepCounter += nj;
        //Rprintf("%d %d %d\n", stepCounter, numEvents, currentPercent);
        currentPercent = printMessage(stepCounter, numEvents, currentPercent);
      }
    }
  }
  fflush(stdout);
  free(wz);
  free(YWX);
  free(ZWX);
  SEXP ans;
  ans = allocList(0);
  return(ans);
}
