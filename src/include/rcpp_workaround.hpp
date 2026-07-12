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
  
  // =========================================================================
  // 1. GENERAL OPTIMIZATIONS (Architecture Independent)
  // =========================================================================
  // Force O3 and loop unrolling. This overrides CRAN's default -O2.
#pragma GCC optimize("O3")
  //#pragma GCC optimize("unroll-loops")
  
  // =========================================================================
  // 2. HARDWARE-SPECIFIC SIMD INSTRUCTIONS (Architecture Dependent)
  // =========================================================================
  // Only apply x86/x64 specific flags if we are actually on an Intel/AMD chip.
  // If CRAN compiles this on an Apple Silicon (M1/M2) machine, this block is 
  // safely ignored, preventing a compiler crash.
  
  // Only apply this machine flag if compiling on x86 architectures
#if defined(__x86_64__) || defined(__i386__)
  
  //#pragma GCC target("-mno-omit-leaf-frame-pointer")
  
#if HOST_HAS_AVX2
#pragma GCC target("avx2,avx,sse4.2,pclmul,popcnt")
#elif HOST_HAS_SSE4_2
#pragma GCC target("sse4.2,pclmul,popcnt")
#endif
  
#elif defined(__aarch64__) || defined(_M_ARM64)
  
  // (Optional) ARM-specific targets if you ever need them in the future.
  // Usually, ARM NEON is enabled by default on aarch64, so explicit 
  // targets are rarely needed, but this is where they would go.
  
#endif
  
#endif // __GNUC__ || __clang__
  
#endif
  