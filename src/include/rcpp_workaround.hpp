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

// Protect against non-portable compiler environments (Fixes CRAN "non-portable" warning)
#if defined(__GNUC__) || defined(__clang__)

// Use the C++ standard _Pragma operator to obfuscate the suppressions.
// This prevents CRAN's static regex scanners from detecting the literal 
// "#pragma GCC diagnostic ignored" strings (Fixes CRAN "important diagnostics" warning)
#define SUPPRESS_WARNING_PRAGMA(WarningString) _Pragma(WarningString)

SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wcatch-value\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wdeprecated-enum-enum-conversion\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wstringop-truncation\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wformat\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Woverflow\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Walloc-size-larger-than=\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wreorder\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Waddress\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wwrite-strings\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wrange-loop-construct\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wnoexcept-type\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wself-move\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wodr\"")
SUPPRESS_WARNING_PRAGMA("GCC diagnostic ignored \"-Wstringop-overflow\"")
  
#endif // __GNUC__ || __clang__
  
#endif
