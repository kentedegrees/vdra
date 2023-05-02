#include <R_ext/Rdynload.h>
#include "compute_cox.h"

void R_init_vdra(DllInfo *dll);

static const R_CallMethodDef Callentries[] = {
  {"ComputeCox",       (DL_FUNC) &compute_cox,      9},
  {NULL, NULL, 0}
};

void R_init_vdra(DllInfo *dll) {
  R_registerRoutines(dll, NULL, Callentries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}




