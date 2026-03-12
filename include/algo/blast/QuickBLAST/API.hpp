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