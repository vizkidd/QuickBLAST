#pragma once

#include <algo/blast/QuickBLAST/commons.hpp>

#include <string>
#include <vector>
#include <thread>
#include <algorithm>   // std::max

#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

using namespace ncbi;
using namespace objects;
using namespace blast;
using namespace std::chrono_literals;
using namespace ncbi::objects;

class DebugHelper
{
public:
  DebugHelper();
  ~DebugHelper();
  static std::string SafeGetNamedScoreAsString(const CSeq_align& align, CSeq_align::EScoreType scoreType);
  static long SafeGetNamedScoreInt(const CSeq_align& align, CSeq_align::EScoreType scoreType, bool &ok);
  static void PrintTSeqLocVector(const TSeqLocVector &vec);
  static void PrintTSeqAlignVector(const TSeqAlignVector &alignments, size_t max_align_sets, size_t max_aligns_per_set);
  static void PrintRecordBatch(const std::shared_ptr<arrow::RecordBatch> &rb, int max_rows);
};
