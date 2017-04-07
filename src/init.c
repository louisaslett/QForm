#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void machine_eps(double*);

static R_NativePrimitiveArgType machine_eps_t[] = {
  REALSXP
};

static R_CMethodDef cMethods[] = {
  {"machine_eps", (DL_FUNC) &machine_eps, 1, machine_eps_t},
  {NULL, NULL, 0}
};

void R_init_QForm(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
