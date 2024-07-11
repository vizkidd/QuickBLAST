#pragma once

#include <chrono>
#include <iostream>
#include <string_view>
#include <map>
#include <filesystem>
#include <functional>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>

std::map<unsigned int, std::shared_ptr<QuickBLAST>> obj_list;

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
