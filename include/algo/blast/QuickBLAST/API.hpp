#pragma once

#include <chrono>
#include <iostream>
#include <string_view>
#include <map>
#include <filesystem>
#include <functional>
#include <memory>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#include "omp.h"
#endif

std::map<unsigned int, Rcpp::XPtr<QuickBLAST>> cppObj_list = {};
struct QuickBLASTHandle
{
    std::shared_ptr<QuickBLAST> ptr;
    // QuickBLAST *ptr;
    int id;
};
/*struct ArrowRBHandle {
    std::shared_ptr<arrow::RecordBatch> ptr;
};*/
struct ArrowRBVHandle
{
    std::shared_ptr<arrow::RecordBatchVector> ptr;
    // arrow::RecordBatchVector *ptr;
};
