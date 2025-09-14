#pragma once

#include <string>
#include <vector>
#include <thread>
#include <algorithm>   // std::max

// #include <R_ext/Print.h>
#include <R.h>
#include <RcppCommon.h>
#include <Rcpp.h>
#include <progress.hpp> //RcppProgress
#include <R_ext/Rdynload.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]

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
