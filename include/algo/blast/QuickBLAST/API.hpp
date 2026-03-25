#pragma once

#include <boost/pool/pool_alloc.hpp>           // boost::pool_allocator
#include <boost/asio/thread_pool.hpp>          // thread_pool
// #include <boost/asio/post.hpp>                 // post
// #include <boost/algorithm/algorithm.hpp>       // algorithm
#include <future>
#include <mutex>
#include <variant>

#include <arrow/array.h>
#include <arrow/buffer.h>
#include <arrow/record_batch.h>
#include <arrow/visit_array_inline.h>
#include <arrow/util/bit_util.h>

// Rcpp header (conversion back to R)
#include <Rcpp.h>

#include <chrono>
#include <iostream>
#include <string_view>
#include <map>
#include <filesystem>
#include <functional>
#include <memory>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <array>
#include <unistd.h>
#include <sys/wait.h>
#include <spawn.h>

#include <filesystem>
#include <system_error>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>

extern char **environ; //posix_spawnp

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#include "omp.h"
#endif

std::map<unsigned int, Rcpp::XPtr<QuickBLAST>> cppObj_list = {};

template <typename T>
using pooled_allocator = boost::pool_allocator<T>;

// common pooled vector types for primitives and strings
template <typename T>
using pooled_vector = std::vector<T, pooled_allocator<T>>;

// convenience factory to create a pooled_vector with reserve
template <typename T>
pooled_vector<T> make_pooled_vector(std::size_t n) {
  pooled_vector<T> v;
  v.reserve(static_cast<typename pooled_vector<T>::size_type>(n));
  return v;
}

enum class ColKind { INT, DOUBLE, STRING, LOGICAL };

struct ColumnData {
  std::string name;
  ColKind kind;
  // use pooled_vector for underlying storage
  pooled_vector<int> int_vals;
  pooled_vector<double> dbl_vals;
  pooled_vector<std::string> str_vals;
  std::vector<bool> str_valid;
  pooled_vector<int> logical_vals; // 0/1 or NA encoded as INT with special NA marker (-2147483648)
  ColumnData() : kind(ColKind::INT) {}
};

// FlatDF: representation of a flattened data.frame in plain C++.
// Each ColumnData must have length = nrows.
struct FlatDF {
  std::vector<ColumnData> cols;
  int64_t nrows = 0;
};

static void add_int_column(FlatDF &out, const std::string &name, int64_t nrows) {
  ColumnData cd;
  cd.name = name;
  cd.kind = ColKind::INT;
  cd.int_vals.resize(static_cast<std::size_t>(nrows), NA_INTEGER); // Pre-fill with NA
  out.cols.push_back(std::move(cd));
}
static void add_double_column(FlatDF &out, const std::string &name, int64_t nrows) {
  ColumnData cd;
  cd.name = name;
  cd.kind = ColKind::DOUBLE;
  cd.dbl_vals.resize(static_cast<std::size_t>(nrows), std::numeric_limits<double>::quiet_NaN());
  out.cols.push_back(std::move(cd));
}
static void add_string_column(FlatDF &out, const std::string &name, int64_t nrows) {
  ColumnData cd;
  cd.name = name;
  cd.kind = ColKind::STRING;
  cd.str_vals.resize(static_cast<std::size_t>(nrows));
  cd.str_valid.resize(static_cast<std::size_t>(nrows), false); // Default to NA
  out.cols.push_back(std::move(cd));
}
static void add_logical_column(FlatDF &out, const std::string &name, int64_t nrows) {
  ColumnData cd;
  cd.name = name;
  cd.kind = ColKind::LOGICAL;
  cd.logical_vals.resize(static_cast<std::size_t>(nrows), INT_MIN);
  out.cols.push_back(std::move(cd));
}


// 1. VISITOR FOR STANDARD 1D COLUMNS
struct PrimitiveVisitor {
  FlatDF& out;
  std::string colname;
  int64_t nrows;
  
  PrimitiveVisitor(FlatDF& o, const std::string& n, int64_t nr) 
    : out(o), colname(n), nrows(nr) {}
  
  // Templated helper for Int8/16/32/64
  template <typename ArrayType>
  arrow::Status ProcessInt(const ArrayType& arr) {
    add_int_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const auto* vals = arr.raw_values();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    // Loop split for max performance (no branching in inner loop)
    if (nulls) {
      for (int64_t i = 0; i < nrows; ++i) {
        if (arrow::bit_util::GetBit(nulls, arr.offset() + i)) {
          cd.int_vals[i] = static_cast<int>(vals[i]);
        }
      }
    } else {
      for (int64_t i = 0; i < nrows; ++i) {
        cd.int_vals[i] = static_cast<int>(vals[i]);
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::Int8Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int16Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int32Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int64Array& arr) { return ProcessInt(arr); }
  
  // Templated helper for Float/Double
  template <typename ArrayType>
  arrow::Status ProcessDouble(const ArrayType& arr) {
    add_double_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const auto* vals = arr.raw_values();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    if (nulls) {
      for (int64_t i = 0; i < nrows; ++i) {
        if (arrow::bit_util::GetBit(nulls, arr.offset() + i)) {
          cd.dbl_vals[i] = static_cast<double>(vals[i]);
        }
      }
    } else {
      for (int64_t i = 0; i < nrows; ++i) {
        cd.dbl_vals[i] = static_cast<double>(vals[i]);
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::FloatArray& arr) { return ProcessDouble(arr); }
  arrow::Status Visit(const arrow::DoubleArray& arr) { return ProcessDouble(arr); }
  
  // Templated helper for Strings / LargeStrings
  template <typename ArrayType>
  arrow::Status ProcessString(const ArrayType& arr) {
    add_string_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    if (nulls) {
      for (int64_t i = 0; i < nrows; ++i) {
        if (arrow::bit_util::GetBit(nulls, arr.offset() + i)) {
          auto sv = arr.GetView(i);
          cd.str_vals[i].assign(sv.data(), sv.size());
          cd.str_valid[i] = true;
        }
      }
    } else {
      for (int64_t i = 0; i < nrows; ++i) {
        auto sv = arr.GetView(i);
        cd.str_vals[i].assign(sv.data(), sv.size());
        cd.str_valid[i] = true;
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::StringArray& arr) { return ProcessString(arr); }
  arrow::Status Visit(const arrow::LargeStringArray& arr) { return ProcessString(arr); }
  
  arrow::Status Visit(const arrow::BooleanArray& arr) {
    add_logical_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    if (nulls) {
      for (int64_t i = 0; i < nrows; ++i) {
        if (arrow::bit_util::GetBit(nulls, arr.offset() + i)) {
          cd.logical_vals[i] = arr.Value(i) ? 1 : 0;
        }
      }
    } else {
      for (int64_t i = 0; i < nrows; ++i) {
        cd.logical_vals[i] = arr.Value(i) ? 1 : 0;
      }
    }
    return arrow::Status::OK();
  }
  
  // Fallback for unsupported types
  arrow::Status Visit(const arrow::Array& arr) {
    add_string_column(out, colname, nrows); 
    return arrow::Status::OK();
  }
};


// 2. VISITOR FOR FLATTENING LIST ELEMENTS
struct ListElementVisitor {
  FlatDF& out;
  std::string colname;
  int64_t nrows;
  const arrow::ListArray& larr;
  int64_t p; // the index within the list to extract
  
  ListElementVisitor(FlatDF& o, const std::string& n, int64_t nr, const arrow::ListArray& la, int64_t idx)
    : out(o), colname(n), nrows(nr), larr(la), p(idx) {}
  
  template <typename ArrayType>
  arrow::Status ProcessInt(const ArrayType& arr) {
    add_int_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const auto* vals = arr.raw_values();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr.value_length(i);
      if (p < len) {
        int64_t idx = larr.value_offset(i) + p;
        if (!nulls || arrow::bit_util::GetBit(nulls, arr.offset() + idx)) {
          cd.int_vals[i] = static_cast<int>(vals[idx]);
        }
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::Int8Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int16Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int32Array& arr) { return ProcessInt(arr); }
  arrow::Status Visit(const arrow::Int64Array& arr) { return ProcessInt(arr); }
  
  template <typename ArrayType>
  arrow::Status ProcessDouble(const ArrayType& arr) {
    add_double_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const auto* vals = arr.raw_values();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr.value_length(i);
      if (p < len) {
        int64_t idx = larr.value_offset(i) + p;
        if (!nulls || arrow::bit_util::GetBit(nulls, arr.offset() + idx)) {
          cd.dbl_vals[i] = static_cast<double>(vals[idx]);
        }
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::FloatArray& arr) { return ProcessDouble(arr); }
  arrow::Status Visit(const arrow::DoubleArray& arr) { return ProcessDouble(arr); }
  
  template <typename ArrayType>
  arrow::Status ProcessString(const ArrayType& arr) {
    add_string_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr.value_length(i);
      if (p < len) {
        int64_t idx = larr.value_offset(i) + p;
        if (!nulls || arrow::bit_util::GetBit(nulls, arr.offset() + idx)) {
          auto sv = arr.GetView(idx);
          cd.str_vals[i].assign(sv.data(), sv.size());
          cd.str_valid[i] = true;
        }
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::StringArray& arr) { return ProcessString(arr); }
  arrow::Status Visit(const arrow::LargeStringArray& arr) { return ProcessString(arr); }
  
  arrow::Status Visit(const arrow::BooleanArray& arr) {
    add_logical_column(out, colname, nrows);
    auto& cd = out.cols.back();
    const uint8_t* nulls = arr.null_bitmap_data();
    
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr.value_length(i);
      if (p < len) {
        int64_t idx = larr.value_offset(i) + p;
        if (!nulls || arrow::bit_util::GetBit(nulls, arr.offset() + idx)) {
          cd.logical_vals[i] = arr.Value(idx) ? 1 : 0;
        }
      }
    }
    return arrow::Status::OK();
  }
  
  arrow::Status Visit(const arrow::Array& arr) {
    add_string_column(out, colname, nrows); 
    return arrow::Status::OK();
  }
};