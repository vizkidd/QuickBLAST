#ifndef RCPP_WORKAROUND_HPP
#define RCPP_WORKAROUND_HPP

// --- R 4.7 / R-devel Workaround ---
// Forward-declare the underlying struct. 
struct SEXPREC;

extern "C" {
  // MUST use SEXPREC* here, because SEXP doesn't exist yet!
  extern SEXPREC* R_NamespaceRegistry;
}
// ----------------------------------

#endif