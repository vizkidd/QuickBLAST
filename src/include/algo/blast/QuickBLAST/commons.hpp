#pragma once
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override

#include <string>
#include <vector>
#include <thread>
#include <utility>
#include <algorithm>   // std::max
#include <mutex>
#include <typeinfo>
#include <iconv.h>
#include <string_view>
#include <cstring>
#include <cerrno>
#include <stdexcept>

// #include <R_ext/Print.h>
// #include <R.h>

// #include <Rinternals.h>
// --- R 4.7 / R-devel Workaround ---
// R Core hid this symbol, but older Rcpp headers still expect it.
// We manually declare it here to satisfy GCC 15's strict template checks.
// Forward-declare the underlying struct. 
struct SEXPREC;

extern "C" {
  // MUST use SEXPREC* here, because SEXP doesn't exist yet!
  extern SEXPREC* R_NamespaceRegistry;
}
// ----------------------------------

#include <RcppCommon.h>
#include <Rcpp.h>
#include <progress.hpp> //RcppProgress
#include <progress_bar.hpp> //RcppProgress
#include <RcppThread.h> //RcppThread
// #include <R_ext/Rdynload.h>

#if defined(_OPENMP)
#include "omp.h"
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppThread)]]
#if defined(_WIN32) || defined(__MINGW32__) || defined(MINGW32) || defined(WIN32)
// #define _WIN32_WINNT _WIN32_WINNT_WIN7
#ifdef QBLIBRARY_EXPORTS
#define QBLIBRARY_API __declspec(dllexport)
#endif
#else
#define QBLIBRARY_API //__declspec(dllimport)
#endif

struct FastaSequenceData
{
    int rec_no = 1;
    std::string header;
    std::string seq;
};
struct BLASTHitData
{
    int rec_no = 1;
    std::vector<std::string_view> col_names();
    std::vector<std::string_view> col_values();
};

struct safe_jthread {
  std::thread t;

  // 1. Default constructor (needed for std::vector sizing)
  safe_jthread() noexcept = default;

  // 2. Accept an already-constructed thread
  explicit safe_jthread(std::thread&& t_) noexcept : t(std::move(t_)) {}

  // 3. Accept a lambda/function with arguments
  template<class Function, class... Args>
  explicit safe_jthread(Function&& f, Args&&... args)
    : t(std::forward<Function>(f), std::forward<Args>(args)...) {}

  // 4. Destructor guarantees join
  ~safe_jthread() {
    if (t.joinable()) {
      t.join();
    }
  }

  // 5. Delete copy operations (Threads CANNOT be copied)
  safe_jthread(const safe_jthread&) = delete;
  safe_jthread& operator=(const safe_jthread&) = delete;

  // 6. Move operations MUST be marked noexcept for std::vector compatibility!
  safe_jthread(safe_jthread&& other) noexcept : t(std::move(other.t)) {}

  safe_jthread& operator=(safe_jthread&& other) noexcept {
    if (this != &other) {
      // If the current thread is running, we must join it before overwriting it
      if (t.joinable()) t.join();
      t = std::move(other.t);
    }
    return *this;
  }

  std::thread& get() { return t; }
};
