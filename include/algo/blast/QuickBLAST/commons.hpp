#pragma once
//#pragma message("Including commons => " __FILE__)
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override

#include <string>
#include <vector>
#include <thread>
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
#include <RcppCommon.h>
#include <Rcpp.h>
#include <progress.hpp> //RcppProgress
#include <progress_bar.hpp> //RcppProgress
#include <RcppThread.h> //RcppThread
#include <R_ext/Rdynload.h>

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
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
