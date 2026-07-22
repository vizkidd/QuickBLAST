// [[Rcpp::plugins(openmp,cpp20,BH)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppThread)]]

#include <algo/blast/QuickBLAST/commons.hpp>
// #include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/API.hpp>


// Helper function to check if a file exists
inline bool FileExists(const std::string& name) {
  struct stat buffer;   
  return (stat(name.c_str(), &buffer) == 0); 
}

using namespace Rcpp;

//' @name isQuickBLASTLoaded
//' @title Check R <-> C++ (FFI) connection
//'
//' @description This function does nothing than check the connection between the R package and C++ libraries
//' @examples
//' QuickBLAST::isQuickBLASTLoaded()
//' @return String that successfully confirms when the package is loaded properly
//' @export
// [[Rcpp::export]]
bool isQuickBLASTLoaded()
{
 ArrowWrapper *testwrap = new ArrowWrapper();
 std::shared_ptr<ArrowWrapper> testwrap_ = std::make_shared<ArrowWrapper>();
 return true;
}

void PrintClock(std::chrono::time_point<std::chrono::high_resolution_clock> start)
{
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl << std::flush;
}

List rm_null(List x)
{
  int n = x.size();
  LogicalVector to_keep(n);
  for (int i = 0; i < n; i++)
  {
    to_keep[i] = !Rf_isNull(x[i]);
  }
  return x[to_keep];
}

SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
{
  
  // Dispatch based on the data type of the array
  if (type->id() == arrow::Type::STRUCT)
  {
    std::shared_ptr<arrow::StructArray> struct_array = std::static_pointer_cast<arrow::StructArray>(array);
    int num_fields = struct_array->num_fields();
    // int num_fields = arrow_struct_num_fields(struct_array);
    
    // Create an Rcpp list to hold the data frames representing each field of the struct
    Rcpp::List struct_list(num_fields);
    Rcpp::CharacterVector names(num_fields);
    
    for (int i = 0; i < num_fields; i++)
    {
      
      std::shared_ptr<arrow::Array> field_array = struct_array->field(i); // arrow_struct_field(struct_array, i);
      std::shared_ptr<arrow::DataType> field_type = type->field(i)->type();
      std::string field_name = type->field(i)->name();
      names[i] = field_name;
      struct_list[i] = Hits2RList_internal(field_array, field_type, field_name);
    }
    
    struct_list.names() = names;
    
    return struct_list;
  }
  else if (type->id() == arrow::Type::LIST)
  {
    std::shared_ptr<arrow::ListArray> list_array = std::static_pointer_cast<arrow::ListArray>(array);
    std::shared_ptr<arrow::DataType> value_type = type->field(0)->type();
    
    // Convert the list array to an Rcpp list
    Rcpp::List list_values(list_array->length());
    Rcpp::CharacterVector names(list_array->length());
    
    for (int i = 0; i < list_array->length(); i++)
    {
      
      auto sublist_array = list_array->values()->Slice(list_array->value_offset(i), list_array->value_length(i));
      
      names[i] = field_name + "[" + std::to_string(i) + "]";
      auto sublist_name = field_name + "[" + std::to_string(i) + "]";
      list_values[i] = Hits2RList_internal(sublist_array, value_type, sublist_name);
    }
    
    list_values.names() = names;
    
    return list_values;
  }
  else if (type->id() == arrow::Type::STRING || type->id() == arrow::Type::LARGE_STRING)
  {
    auto string_array = std::static_pointer_cast<arrow::StringArray>(array);
    
    Rcpp::StringVector strings(string_array->length());
    
    for (int i = 0; i < string_array->length(); ++i)
    {
      if (string_array->IsValid(i))
      {
        strings[i] = Rcpp::String(string_array->GetString(i));
      }
    }
    
    return strings;
  }
  else if (type->id() == arrow::Type::INT32)
  {
    auto int_array = std::static_pointer_cast<arrow::Int32Array>(array);
    
    Rcpp::IntegerVector ints(int_array->length());
    
    for (int i = 0; i < int_array->length(); ++i)
    {
      if (int_array->IsValid(i))
      {
        ints[i] = int_array->Value(i);
      }
    }
    
    return ints;
  }
  else if (type->id() == arrow::Type::INT64)
  {
    auto int_array = std::static_pointer_cast<arrow::Int64Array>(array);
    
    Rcpp::IntegerVector ints(int_array->length());
    
    for (int i = 0; i < int_array->length(); ++i)
    {
      if (int_array->IsValid(i))
      {
        ints[i] = int_array->Value(i);
      }
    }
    
    return ints;
  }
  else if (type->id() == arrow::Type::DOUBLE)
  { // Use arrow::Type::DOUBLE instead of FLOAT64
    
    auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
    Rcpp::NumericVector doubles(double_array->length());
    
    for (int i = 0; i < double_array->length(); ++i)
    {
      if (double_array->IsValid(i))
      {
        doubles[field_name] = double_array->Value(i);
      }
    }
    
    return doubles;
  }
  else
  {
    // For other data types that don't have a direct conversion, return R_NilValue (NA)
    return R_NilValue;
  }
}

SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
{
  
  // Assuming the schema of the RecordBatch is accessible here
  auto rb_schema = rb->schema();
  
  // Convert each column of the RecordBatch to R objects and store in a list
  Rcpp::List result_list(rb_schema->num_fields()); // arrow_schema_num_fields(rb_schema)
  
  for (int i = 0; i < rb_schema->num_fields(); ++i) // arrow_schema_num_fields(rb_schema)
  {
    RcppThread::checkUserInterrupt();
    auto array = rb->column(i);
    auto field_type = rb_schema->field(i)->type(); // arrow_schema_field_type(rb_schema, i);
    auto field_name = rb_schema->field(i)->name(); // arrow_schema_field_name(rb_schema, i);
    result_list[i] = Hits2RList_internal(array, field_type, field_name);
  }
  
  return result_list;
}

SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector)
{
  Rcpp::List result_list(rb_vector.size());
  
  // Traverse the vector of RecordBatches and convert each RecordBatch
  for (size_t i = 0; i < rb_vector.size(); ++i)
  {
    RcppThread::checkUserInterrupt();
    std::shared_ptr<arrow::RecordBatch> rb = rb_vector[i];
    result_list[i] = Hits2RList(rb);
  }
  
  return result_list;
}

// Helpers: string join for column name components
static std::string join_col(const std::string &a, const std::string &b) {
  if (a.empty()) return b;
  return a + "_" + b;
}

// Forward declarations
static void CollectFlattenedColumns(const std::shared_ptr<arrow::Array>& array,
                                    const std::shared_ptr<arrow::DataType>& type,
                                    const std::string& prefix,
                                    std::vector<std::string>& col_order,
                                    std::unordered_map<std::string, Rcpp::RObject>& cols);

static Rcpp::RObject ConvertPrimitiveValueToR(const std::shared_ptr<arrow::Array>& arr,
                                              const std::shared_ptr<arrow::DataType>& dtype,
                                              int64_t idx);

static Rcpp::List RecordBatchVectorToFlattenedDFList_parallel(const arrow::RecordBatchVector &rbv);

Rcpp::IntegerVector Int8ArrayToR(const std::shared_ptr<arrow::Int8Array>& iarr) {
  int64_t n = iarr->length();
  Rcpp::IntegerVector out(static_cast<int>(n));
  const int8_t* raw = iarr->raw_values(); // zero-overhead access
  for (int64_t i = 0; i < n; ++i) {
    out[i] = iarr->IsValid(i) ? static_cast<int>(raw[i]) : NA_INTEGER;
  }
  return out;
}

Rcpp::IntegerVector Int16ArrayToR(const std::shared_ptr<arrow::Int16Array>& iarr) {
  int64_t n = iarr->length();
  Rcpp::IntegerVector out(static_cast<int>(n));
  const int16_t* raw = iarr->raw_values(); // zero-overhead access
  for (int64_t i = 0; i < n; ++i) {
    out[i] = iarr->IsValid(i) ? static_cast<int>(raw[i]) : NA_INTEGER;
  }
  return out;
}

Rcpp::IntegerVector Int32ArrayToR(const std::shared_ptr<arrow::Int32Array>& iarr) {
  int64_t n = iarr->length();
  Rcpp::IntegerVector out(static_cast<int>(n));
  const int32_t* raw = iarr->raw_values(); // zero-overhead access
  for (int64_t i = 0; i < n; ++i) {
    out[i] = iarr->IsValid(i) ? static_cast<int>(raw[i]) : NA_INTEGER;
  }
  return out;
}

Rcpp::NumericVector Int64ArrayToR(const std::shared_ptr<arrow::Int64Array>& iarr) {
  int64_t n = iarr->length();
  Rcpp::NumericVector out(static_cast<int>(n));
  const int64_t* raw = iarr->raw_values(); // zero-overhead access
  for (int64_t i = 0; i < n; ++i) {
    out[i] = iarr->IsValid(i) ? static_cast<double>(raw[i]) : NA_INTEGER;
  }
  return out;
}

Rcpp::NumericVector DoubleArrayToR(const std::shared_ptr<arrow::DoubleArray>& darr) {
  int64_t n = darr->length();
  Rcpp::NumericVector out(static_cast<int>(n));
  const double* raw = darr->raw_values();
  for (int64_t i = 0; i < n; ++i) {
    out[i] = darr->IsValid(i) ? raw[i] : NA_REAL;
  }
  return out;
}

Rcpp::NumericVector FloatArrayToR(const std::shared_ptr<arrow::FloatArray>& darr) {
  int64_t n = darr->length();
  Rcpp::NumericVector out(static_cast<float>(n));
  const float* raw = darr->raw_values();
  for (int64_t i = 0; i < n; ++i) {
    out[i] = darr->IsValid(i) ? static_cast<double>(raw[i]) : NA_REAL;
  }
  return out;
}

// Convert an Arrow array of primitive/scalar type (length == nrows) to an R vector
static Rcpp::RObject ArrayPrimitivesToR(const std::shared_ptr<arrow::Array>& array,
                                        const std::shared_ptr<arrow::DataType>& type)
{
  using arrow::Int8Array; using arrow::Int16Array; using arrow::Int32Array;
  using arrow::Int64Array; using arrow::DoubleArray; using arrow::FloatArray;
  using arrow::StringArray; using arrow::BooleanArray;
  
  int64_t n = array->length();
  switch (type->id()) {
  case arrow::Type::STRING:{
    auto sarr = std::static_pointer_cast<arrow::StringArray>(array);
    // get actual length from the array (int64_t)
    int64_t n_local = sarr->length();
    // Safety: Rcpp::StringVector constructor expects an 'int' (32-bit) on many builds.
    if (n_local > static_cast<int64_t>(std::numeric_limits<int>::max())) {
      Rcpp::Rcerr << "Arrow string array has too many elements to convert to an R character vector." << std::endl <<std::flush;
    }
    int n_int = static_cast<int>(n_local);
    // std::cout << "converting Arrow StringArray of length " << n_local << std::endl << std::flush;
    Rcpp::StringVector out(n_int);
    // If you expect extremely large individual strings, you might check sarr->value(i).size() here
    for (int64_t i = 0; i < n_local; ++i) {
      int idx = static_cast<int>(i);
      if (sarr->IsValid(i)) {
        // Prefer GetView() to avoid an extra allocation if possible.
        std::string_view sv = sarr->GetView(i);
        // std::cout << "Row : " << idx << std::endl << std::flush; //DEBUG
        // std::cout << "Row Size: " << sv.size() << std::endl << std::flush; //DEBUG
        // std::cout << "Row Data: " << sv.data() << std::endl << std::flush; //DEBUG
        // Construct Rcpp::String from the view (this makes a copy once)
        out[idx] = Rcpp::String(std::string(sv.data(), sv.size()));
      } else {
        out[idx] = NA_STRING;
      }
    }
    // for( int i = 0; i < out.size(); i++ ){
    //   std::cout << "Element " << i << ": " << out(i) << std::endl << std::flush; //DEBUG
    // }
    // std::cout << "here 3.1.4" << std::endl << std::flush; //DEBUG
    return out;
  }
  case arrow::Type::LARGE_STRING: {
    auto sarr = std::static_pointer_cast<arrow::LargeStringArray>(array);
    // get actual length from the array (int64_t)
    int64_t n_local = sarr->length();
    // std::cout << "here 3.0.1" << std::endl << std::flush; //DEBUG
    // Safety: Rcpp::StringVector constructor expects an 'int' (32-bit) on many builds.
    if (n_local > static_cast<int64_t>(std::numeric_limits<int>::max())) {
      Rcpp::Rcerr << "Arrow string array has too many elements to convert to an R character vector."<< std::endl <<std::flush;
    }
    int n_int = static_cast<int>(n_local);
    // std::cout << "converting Arrow StringArray of length " << n_local << std::endl << std::flush;
    Rcpp::StringVector out(n_int);
    // If you expect extremely large individual strings, you might check sarr->value(i).size() here
    for (int64_t i = 0; i < n_local; ++i) {
      int idx = static_cast<int>(i);
      if (sarr->IsValid(i)) {
        // Prefer GetView() to avoid an extra allocation if possible.
        std::string_view sv = sarr->GetView(i);
        // std::cout << "Row : " << idx << std::endl << std::flush; //DEBUG
        // std::cout << "Row Size: " << sv.size() << std::endl << std::flush; //DEBUG
        // std::cout << "Row Data: " << sv.data() << std::endl << std::flush; //DEBUG
        // Construct Rcpp::String from the view (this makes a copy once)
        out[idx] = Rcpp::String(std::string(sv.data(), sv.size()));
      } else {
        out[idx] = NA_STRING;
      }
    }
    // for( int i = 0; i < out.size(); i++ ){
    //   std::cout << "Element " << i << ": " << out(i) << std::endl << std::flush; //DEBUG
    // }
    // std::cout << "here 3.0.4" << std::endl << std::flush; //DEBUG
    return out;
  }
  case arrow::Type::BOOL: {
    auto barr = std::static_pointer_cast<BooleanArray>(array);
    Rcpp::LogicalVector out(n);
    for (int64_t i = 0; i < n; ++i) out[i] = barr->IsValid(i) ? static_cast<int>(barr->Value(i)) : NA_LOGICAL;
    return out;
  }
  case arrow::Type::INT8: {
    auto iarr = std::static_pointer_cast<Int8Array>(array);
    // Rcpp::IntegerVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = iarr->IsValid(i) ? static_cast<int>(iarr->Value(i)) : NA_INTEGER;
    // return out;
    return RObject(Int8ArrayToR(iarr));
  }
  case arrow::Type::INT16: {
    auto iarr = std::static_pointer_cast<arrow::Int16Array>(array);
    // Rcpp::IntegerVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = iarr->IsValid(i) ? static_cast<int>(iarr->Value(i)) : NA_INTEGER;
    // return out;
    return RObject(Int16ArrayToR(iarr));
  }
  case arrow::Type::INT32: {
    auto iarr = std::static_pointer_cast<Int32Array>(array);
    return Rcpp::RObject(Int32ArrayToR(iarr));
    // Rcpp::IntegerVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = iarr->IsValid(i) ? iarr->Value(i) : NA_INTEGER;
    // return out;
  }
  case arrow::Type::INT64: {
    auto iarr = std::static_pointer_cast<Int64Array>(array);
    // // R doesn't have 64-bit ints by default; use double
    // Rcpp::NumericVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = iarr->IsValid(i) ? static_cast<double>(iarr->Value(i)) : NA_REAL;
    // return out;
    return Rcpp::RObject(Int64ArrayToR(iarr));
  }
  case arrow::Type::FLOAT: {
    auto farr = std::static_pointer_cast<FloatArray>(array);
    // Rcpp::NumericVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = farr->IsValid(i) ? static_cast<double>(farr->Value(i)) : NA_REAL;
    // return out;
    return Rcpp::RObject(FloatArrayToR(farr));
  }
  case arrow::Type::DOUBLE: {
    auto darr = std::static_pointer_cast<DoubleArray>(array);
    return Rcpp::RObject(DoubleArrayToR(darr));
    // Rcpp::NumericVector out(n);
    // for (int64_t i = 0; i < n; ++i) out[i] = darr->IsValid(i) ? darr->Value(i) : NA_REAL;
    // return out;
  }
  default: {
    // Unknown primitive: hand back list of NAs
    Rcpp::List out(n);
    for (int64_t i = 0; i < n; ++i) out[i] = R_NilValue;
    return out;
  }
  }
}

// CollectFlattenedColumns: recursively populates cols map and col_order
static void CollectFlattenedColumns(const std::shared_ptr<arrow::Array>& array,
                                    const std::shared_ptr<arrow::DataType>& type,
                                    const std::string& prefix,
                                    std::vector<std::string>& col_order,
                                    std::unordered_map<std::string, Rcpp::RObject>& cols)
{
  if (!array) return;
  int64_t nrows = array->length();
  
  if (type->id() == arrow::Type::STRUCT) {
    auto sarr = std::static_pointer_cast<arrow::StructArray>(array);
    int nf = sarr->num_fields();
    for (int f = 0; f < nf; ++f) {
      std::string child_name = type->field(f)->name();
      std::string colname = child_name; //join_col(prefix, child_name);
      // std::cout << "Colname: " << colname << std::endl << std::flush; //DEBUG
      auto child_array = sarr->field(f);
      auto child_type = type->field(f)->type();
      CollectFlattenedColumns(child_array, child_type, colname, col_order, cols);
    }
    return;
  }
  
  if (type->id() == arrow::Type::LIST) {
    auto larr = std::static_pointer_cast<arrow::ListArray>(array);
    auto val_type = type->field(0)->type();
    auto values = larr->values();           // array of all elements flattened
    // compute max length across rows
    int64_t maxlen = 0;
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (len > maxlen) maxlen = len;
    }
    // For each position p build a column
    for (int64_t p = 0; p < maxlen; ++p) {
      // if element type is primitive -> create atomic column
      if (val_type->id() == arrow::Type::STRUCT) {
        // We need to expand each subfield of the struct element at position p:
        auto vals_struct = std::static_pointer_cast<arrow::StructArray>(values);
        int sub_nf = vals_struct->num_fields();
        for (int sf = 0; sf < sub_nf; ++sf) {
          std::string sf_name = val_type->field(sf)->name();
          std::string colname = prefix + "_" + std::to_string(p) + "_" + sf_name;
          // build atomic vector for this subfield at index offset+p per row
          // find subfield array
          auto subfield_array = vals_struct->field(sf);
          // create atomic vector depending on subfield type
          auto sub_type = val_type->field(sf)->type();
          // create appropriate R vector and fill
          int64_t n = nrows;
          // Use ArrayPrimitivesToR on a per-index basis: iterate rows and read value at values_index = offset + p
          // For performance, handle common scalar types manually:
          if (sub_type->id() == arrow::Type::STRING || sub_type->id() == arrow::Type::LARGE_STRING) {
            Rcpp::StringVector col(n);
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                if (subfield_array->IsValid(idx)) {
                  col[i] = Rcpp::String(std::static_pointer_cast<arrow::StringArray>(subfield_array)->GetString(idx));
                } else col[i] = NA_STRING;
              } else col[i] = NA_STRING;
            }
            cols[colname] = col;
            col_order.push_back(colname);
          } else if (sub_type->id() == arrow::Type::INT8 || sub_type->id() == arrow::Type::INT16 || sub_type->id() == arrow::Type::INT32 || sub_type->id() == arrow::Type::INT64) {
            Rcpp::IntegerVector col(n);
            auto sub_i32 = std::static_pointer_cast<arrow::Int64Array>(subfield_array); // safe if int32; else use Value and cast
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                if (subfield_array->IsValid(idx)) col[i] = static_cast<int>(subfield_array->type()->id() == arrow::Type::INT32 ? std::static_pointer_cast<arrow::Int32Array>(subfield_array)->Value(idx) : static_cast<int>(std::static_pointer_cast<arrow::Int64Array>(subfield_array)->Value(idx)));
                else col[i] = NA_INTEGER;
              } else col[i] = NA_INTEGER;
            }
            cols[colname] = col;
            col_order.push_back(colname);
          } else if (sub_type->id() == arrow::Type::DOUBLE || sub_type->id() == arrow::Type::FLOAT) {
            Rcpp::NumericVector col(n);
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                if (subfield_array->IsValid(idx)) {
                  if (sub_type->id() == arrow::Type::DOUBLE) col[i] = std::static_pointer_cast<arrow::DoubleArray>(subfield_array)->Value(idx);
                  else col[i] = static_cast<double>(std::static_pointer_cast<arrow::FloatArray>(subfield_array)->Value(idx));
                } else col[i] = NA_REAL;
              } else col[i] = NA_REAL;
            }
            cols[colname] = col;
            col_order.push_back(colname);
          } else {
            // Fallback: build list column per row (nil or scalar)
            Rcpp::List col(n);
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                // convert single-element slice
                auto single = subfield_array->Slice(idx, 1);
                col[i] = ArrayPrimitivesToR(single, sub_type);
              } else col[i] = R_NilValue;
            }
            cols[colname] = col;
            col_order.push_back(colname);
          }
        } // end subfields
      } else {
        // val_type is primitive: create column prefix_p of appropriate atomic type
        std::string colname = prefix + "_" + std::to_string(p);
        int64_t n = nrows;
        switch (val_type->id()) {
        case arrow::Type::STRING:
        case arrow::Type::LARGE_STRING: {
          // std::cout << "here 4.0.1" << std::endl << std::flush; //DEBUG
          Rcpp::StringVector col(n);
          // std::cout << "here 4.0.2" << std::endl << std::flush; //DEBUG
          auto sval = std::static_pointer_cast<arrow::StringArray>(values);
          // std::cout << "here 4.0.3" << std::endl << std::flush; //DEBUG
          for (int64_t i = 0; i < n; ++i) {
            int64_t len = larr->value_length(i);
            if (p < len) {
              int64_t idx = larr->value_offset(i) + p;
              // std::cout << "here 4.0.4.1" << std::endl << std::flush; //DEBUG
              col[i] = sval->IsValid(idx) ? Rcpp::String(sval->GetString(idx)) : NA_STRING;
              // std::cout << "here 4.0.4.2" << std::endl << std::flush; //DEBUG
            } else col[i] = NA_STRING;
          }
          // std::cout << "here 4.0.5" << std::endl << std::flush; //DEBUG
          cols[colname] = col;
          // std::cout << "here 4.0.6" << std::endl << std::flush; //DEBUG
          col_order.push_back(colname);
          break;
        }
        case arrow::Type::BOOL: {
          Rcpp::LogicalVector col(n);
          auto barr = std::static_pointer_cast<arrow::BooleanArray>(values);
          for (int64_t i = 0; i < n; ++i) {
            int64_t len = larr->value_length(i);
            if (p < len) {
              int64_t idx = larr->value_offset(i) + p;
              col[i] = barr->IsValid(idx) ? static_cast<int>(barr->Value(idx)) : NA_LOGICAL;
            } else col[i] = NA_LOGICAL;
          }
          cols[colname] = col;
          col_order.push_back(colname);
          break;
        }
        case arrow::Type::INT8:
        case arrow::Type::INT16:
        case arrow::Type::INT32: {
          Rcpp::IntegerVector col(n);
          // handle as Int64Array then cast, as Arrow may be mixed
          for (int64_t i = 0; i < n; ++i) {
            int64_t len = larr->value_length(i);
            if (p < len) {
              int64_t idx = larr->value_offset(i) + p;
              if (values->IsValid(idx)) {
                // try common int arrays
                if (values->type_id() == arrow::Type::INT32) {
                  col[i] = std::static_pointer_cast<arrow::Int32Array>(values)->Value(idx);
                } else if (values->type_id() == arrow::Type::INT64) {
                  col[i] = static_cast<int>(std::static_pointer_cast<arrow::Int64Array>(values)->Value(idx));
                } else if (values->type_id() == arrow::Type::INT16) {
                  col[i] = static_cast<int>(std::static_pointer_cast<arrow::Int16Array>(values)->Value(idx));
                } else {
                  col[i] = NA_INTEGER;
                }
              } else col[i] = NA_INTEGER;
            } else col[i] = NA_INTEGER;
          }
          cols[colname] = col;
          col_order.push_back(colname);
          break;
        }
        case arrow::Type::INT64: {
          Rcpp::NumericVector col(n);
          for (int64_t i = 0; i < n; ++i) {
            int64_t len = larr->value_length(i);
            if (p < len) {
              int64_t idx = larr->value_offset(i) + p;
              if (values->IsValid(idx)) col[i] = static_cast<double>(std::static_pointer_cast<arrow::Int64Array>(values)->Value(idx));
              else col[i] = NA_REAL;
            } else col[i] = NA_REAL;
          }
          cols[colname] = col;
          col_order.push_back(colname);
          break;
        }
        case arrow::Type::FLOAT:
        case arrow::Type::DOUBLE: {
          Rcpp::NumericVector col(n);
          if (values->type_id() == arrow::Type::DOUBLE) {
            auto darr = std::static_pointer_cast<arrow::DoubleArray>(values);
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                col[i] = darr->IsValid(idx) ? darr->Value(idx) : NA_REAL;
              } else col[i] = NA_REAL;
            }
          } else {
            auto farr = std::static_pointer_cast<arrow::FloatArray>(values);
            for (int64_t i = 0; i < n; ++i) {
              int64_t len = larr->value_length(i);
              if (p < len) {
                int64_t idx = larr->value_offset(i) + p;
                col[i] = farr->IsValid(idx) ? static_cast<double>(farr->Value(idx)) : NA_REAL;
              } else col[i] = NA_REAL;
            }
          }
          cols[colname] = col;
          col_order.push_back(colname);
          break;
        }
        default: {
          // fallback: fill list column with single-element conversions
          Rcpp::List col(n);
          for (int64_t i = 0; i < n; ++i) {
            int64_t len = larr->value_length(i);
            if (p < len) {
              int64_t idx = larr->value_offset(i) + p;
              auto single = values->Slice(idx, 1);
              col[i] = ArrayPrimitivesToR(single, val_type);
            } else col[i] = R_NilValue;
          }
          cols[colname] = col;
          col_order.push_back(colname);
        }
        } // switch val_type
      } // else primitive val_type
    } // for p
    return;
  } // LIST
  
  // Primitive scalar at this level: just convert whole array to R atomic vector
  {
    Rcpp::RObject v = ArrayPrimitivesToR(array, type);
    // add to map if not already present
    if (cols.find(prefix) == cols.end()) {
      cols[prefix] = v;
      col_order.push_back(prefix);
    }
    return;
  }
}

// Convert a single RecordBatch into a flattened data.frame (Rcpp::List with df attrs)
static Rcpp::List RecordBatchToFlattenedDF(const std::shared_ptr<arrow::RecordBatch>& rb)
{
  if (!rb) {
    // allocate on the heap so XPtr owns a valid pointer
    Rcpp::List sentinel = Rcpp::List::create();
    sentinel["mangled_rb"] = 0;
    // return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
    return sentinel;
  }
  
  if (rb->num_rows() == 0) {
    // allocate on the heap so XPtr owns a valid pointer
    Rcpp::List sentinel =  Rcpp::List::create();
    sentinel["number_of_hits"] = 0;
    // return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
    return sentinel;
  }
  
  auto schema = rb->schema();
  int nfields = schema->num_fields();
  int64_t nrows = rb->num_rows();
  
  // if (nrows == 0) {
  //   // allocate on the heap so XPtr owns a valid pointer
  //   Rcpp::List* sentinel = new Rcpp::List();
  //   (*sentinel)["number_of_hits"] = 0;
  //   return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
  // }
  
  // map of column name -> R vector
  std::unordered_map<std::string, Rcpp::RObject> cols;
  std::vector<std::string> col_order;
  
  for (int i = 0; i < nfields; ++i) {
    RcppThread::checkUserInterrupt();
    std::string fname = schema->field(i)->name();
    auto ftype = schema->field(i)->type();
    auto arr = rb->column(i);
    CollectFlattenedColumns(arr, ftype, fname, col_order, cols);
  }
  
  if (col_order.empty()) return Rcpp::List::create(); 
  
  // Build Rcpp::List in col_order
  int ncols = static_cast<int>(col_order.size());
  Rcpp::List out(ncols);
  Rcpp::CharacterVector names(ncols);
  for (int j = 0; j < ncols; ++j) {
    const std::string &cn = col_order[j];
    auto it = cols.find(cn);
    if (it != cols.end()) {
      out[j] = std::move(it->second); // moves RObject out of the map (leaves map entry empty)
    } else {
      out[j] = R_NilValue;
    }
    names[j] = cn;
  }
  // If you moved values out of cols and don't need it anymore:
  cols.clear();
  
  out.attr("names") = names;
  out.attr("class") = Rcpp::CharacterVector::create("data.frame");
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(nrows));
  // return Rcpp::XPtr<Rcpp::List>(out, true);
  return out;
}

// Entrypoint: RecordBatchVector -> list of flattened data.frames
// [[Rcpp::export]]
Rcpp::List RecordBatchVectorToFlattenedDFList(SEXP rbv_sexp) {
  // The package likely has direct access to the C++ RecordBatchVector pointer in code.
  // Here we assume rbv_sexp wraps a pointer to shared_ptr<arrow::RecordBatchVector> or you
  // will call this function from C++ with the pointer directly. For convenience we implement
  // a wrapper that does not try to parse SEXP: in your code call RecordBatchVectorToFlattenedDFList directly.
  Rcpp::Rcerr << "This wrapper is intended to be called from C++ with a std::shared_ptr<arrow::RecordBatchVector>. Use the C++ entrypoint instead." << std::endl <<std::flush;
  Rcpp::List sentinel = Rcpp::List::create();
  sentinel["number_of_hits"] = 0;
  // return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
  return sentinel;
}

// Alternate C++ entrypoint (call this from your C++ code where you have rbv):
static Rcpp::List RecordBatchVectorToFlattenedDFList_cpp(const arrow::RecordBatchVector& rbv)
{
  // if (!rbv) {
  //   // allocate on the heap so XPtr owns a valid pointer
  //   Rcpp::List sentinel = Rcpp::List::create();
  //   sentinel["rbv_hits"] = 0;
  //   // return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
  //   return sentinel;
  // }
  
  if (rbv.size() == 0) {
    // allocate on the heap so XPtr owns a valid pointer
    Rcpp::List sentinel = Rcpp::List::create();
    sentinel["rbv_hits"] = 0;
    // return Rcpp::XPtr<Rcpp::List>(sentinel, true); // 'true' => XPtr will delete on GC
    return sentinel;
  }
  Rcpp::List out(rbv.size());
  for (size_t i = 0; i < rbv.size(); ++i) {
    RcppThread::checkUserInterrupt();
    // std::cout << "Processing RecordBatch: " << i << std::endl << std::flush; //DEBUG
    // std::cout << (*rbv)[i]->ToString() << std::endl << std::flush; //DEBUG
    // Rcpp::XPtr rb((*rbv)[i].get(), true);
    out[i] = RecordBatchToFlattenedDF(rbv[i]);
  }
  arrow::MemoryPool* mp = arrow::default_memory_pool();
  // #ifdef ARROW_HAVE_MEMORY_POOL_RELEASE
  mp->ReleaseUnused();
  // #endif
  // return Rcpp::XPtr<Rcpp::List>(out, true);
  return out;
}

// std::vector<std::string> getFilesInDir_boost(const std::string &folderPath, const std::string &extension) {
//   std::vector<std::string> outFiles;
//   boost::system::error_code ec;
//   boost::filesystem::path p(folderPath);
//   
//   if (!boost::filesystem::exists(p, ec) || ec) return outFiles;
//   boost::filesystem::directory_iterator it(p, ec), end;
//   for (; it != end && !ec; ++it) {
//     try {
//       if (!boost::filesystem::is_regular_file(*it, ec) || ec) { ec.clear(); continue; }
//       if (!extension.empty()) {
//         if (it->path().extension() == extension)
//           outFiles.emplace_back(it->path().filename().string());
//       } else {
//         outFiles.emplace_back(it->path().filename().string());
//       }
//     } catch (const std::exception &e) {
//       // continue on errors; do not let exceptions escape to R
//       ec.clear();
//       continue;
//     }
//   }
//   return outFiles;
// }

std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension)
{
  std::vector<std::string> outFiles;
  std::error_code ec;
  
  if (!extension.empty())
  {
    for (const auto &entry : std::filesystem::directory_iterator(folderPath, ec))
    {
      if (ec) {
        Rcpp::stop("[getFilesInDir()] Cannot access directory: " + folderPath);
      }
      if (entry.is_regular_file() && entry.path().extension() == extension)
      {
        outFiles.emplace_back(entry.path().filename().string());
      }
    }
  }
  else
  {
    for (const auto &entry : std::filesystem::directory_iterator(folderPath,ec))
    {
      if (ec) {
        Rcpp::stop("[getFilesInDir()] Cannot access directory: " + folderPath);
      }
      if (entry.is_regular_file())
      {
        outFiles.emplace_back(entry.path().filename().string());
      }
    }
  }
  return outFiles;
}

std::string getFilenameWithoutExtension(const std::string &filename)
{
  // Find the position of the last dot (.)
  size_t dotPos = filename.find_last_of('.');
  if (dotPos != std::string::npos)
  {
    // Return the substring from the beginning up to the last dot position
    return filename.substr(0, dotPos);
  }
  else
  {
    // If no dot is found, return the entire filename as it is
    return filename;
  }
}

std::vector<std::string> getFilesInDir_std(const std::string &folderPath, const std::string &extension) {
  std::vector<std::string> outFiles;
  std::error_code ec;
  std::filesystem::path p(folderPath);
  
  if (!std::filesystem::exists(p, ec) || ec) return outFiles;
  
  std::filesystem::directory_options options = std::filesystem::directory_options::skip_permission_denied;
  std::filesystem::directory_iterator it(p, options, ec);
  if (ec) return outFiles;
  
  std::filesystem::directory_iterator end;
  while (it != end) {
    std::error_code ec_entry;
    if (it->is_regular_file(ec_entry) && !ec_entry) {
      if (extension.empty() || it->path().extension() == extension) {
        outFiles.push_back(it->path().filename().string());
      }
    }
    std::error_code ec_inc;
    it.increment(ec_inc);
    if (ec_inc) break;
  }
  return outFiles;
}

// 2. C++17 Standard String Trimmer (Replaces boost::algorithm::trim)
auto trim_string = [](std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
};

std::string ConvertBLASTOptions2String_std(SEXP options) {
  if (TYPEOF(options) == STRSXP) return Rcpp::as<std::string>(options);
  
  Rcpp::List lst = Rcpp::as<Rcpp::List>(options);
  Rcpp::CharacterVector names = lst.names();
  std::string out;
  
  for (int i = 0; i < lst.size(); ++i) {
    std::string name = Rcpp::as<std::string>(names[i]);
    std::string val = Rcpp::as<std::string>(lst[i]);
    
    trim_string(name);
    trim_string(val);
    
    out += "-" + name + " " + val + " ";
  }
  
  trim_string(out);
  return out;
}

// std::string ConvertBLASTOptions2String_boost(SEXP options) {
//   if (TYPEOF(options) == STRSXP) return Rcpp::as<std::string>(options);
//   
//   Rcpp::List lst = Rcpp::as<Rcpp::List>(options);
//   Rcpp::CharacterVector names = lst.names();
//   std::string out;
//   for (int i = 0; i < lst.size(); ++i) {
//     std::string name = Rcpp::as<std::string>(names[i]);
//     std::string val = Rcpp::as<std::string>(lst[i]);
//     // trim whitespace
//     boost::algorithm::trim(name);
//     boost::algorithm::trim(val);
//     out += "-" + name + " " + val + " ";
//   }
//   boost::algorithm::trim(out);
//   return out;
// }

std::string ConvertBLASTOptions2String(SEXP options)
{
  std::string options_;
  int typ = TYPEOF(options);
  switch (typ)
  {
  case LISTSXP:
  case VECSXP:
  {
    Rcpp::List options_list = Rcpp::as<Rcpp::List>(options);
    Rcpp::StringVector options_names = options_list.names();
    
    for (int i = 0; i < options_list.size(); ++i)
    {
      std::string name = Rcpp::as<std::string>(options_names[i]);
      SEXP value = options_list[i];
      
      options_ += "-" + name + " " + as<std::string>(value) + " ";
    }
    break;
  }
  case STRSXP:
  {
    
    options_ = Rcpp::as<std::string>(options);
    break;
  }
  default:
    Rcpp::Rcerr << "Only named list or string allowed for BLAST options : Check QuickBLAST::GetAvailableBLASTOptions()" << std::endl << std::flush;
    break;
  }
  return options_;
}

//' @name CreateQuickBLASTInstance
//'
//' @title Create new QuickBLAST instance with seq_type, strand, program and BLAST options.
//'
//' @description Create a new QuickBLAST C++ object with seq_type, strand, program and BLAST options, which can be used in QuickBLAST::BLAST2Files() and QuickBLAST::BLAST2Seqs()
//'
//' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param seq_type (int) 0 - (eNucleotide) (OR) 1 - (eProtein)
//' @param strand (int) 0 - (ePlus) (OR) 1 - (eMinus)
//' @param program (string) Name of the BLAST program
//' @param options (string (or) Named List) List of BLAST options - check QuickBLAST::GetAvailableBLASTOptions(). String should be of the format "-option1 value1 -option2 value2". If empty, default values (per program) are used.
//' @param save_sequences (bool) Save Full Sequences to output?. (Default: FALSE)
//' @param save_hsp_sequences (bool) Save HSP Sequences to output?. (Default: FALSE)
//' @return (\code{Rcpp::XPtr<QuickBLAST>}) Pointer to a QuickBLAST Instance (Cannot be used in R)
//' 
//' @note Set save_sequences AND/OR save_hsp_sequences when using Genomes
//' 
//' @examples
//' \donttest{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = FALSE,
//'   save_hsp_sequences = FALSE
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP CreateQuickBLASTInstance(const int seq_type, const int strand, SEXP program, SEXP options = R_NilValue, const bool save_sequences = false, const bool save_hsp_sequences = false) 
{
 try
 {
   
   QuickBLAST::ESeqType seq_type_ = static_cast<QuickBLAST::ESeqType>(seq_type);
   QuickBLAST::EStrand strand_ = static_cast<QuickBLAST::EStrand>(strand);
   std::string program_ = Rcpp::as<std::string>(program);
   
   std::string options_ = options == R_NilValue ? "" : ConvertBLASTOptions2String_std(options);
   bool save_sequences_ = save_sequences; //static_cast<bool>(save_sequences);
   bool save_hsp_sequences_ = save_hsp_sequences;
   QuickBLAST *objPtr = new QuickBLAST(seq_type_, strand_, program_, options_, save_sequences_, save_hsp_sequences_);
   Rcpp::XPtr<QuickBLAST> ptr(objPtr, true);
   ptr.attr("seq_type") = seq_type;
   ptr.attr("strand") = strand;
   ptr.attr("program") = program_;
   ptr.attr("options") = options_;
   ptr.attr("save_sequences") = save_sequences_;
   ptr.attr("save_hsp_sequences") = save_hsp_sequences_;
   // ptr.attr("num_threads") = threads;
   // ptr.attr("class") = "QuickBLAST_XPtr";
   Rf_classgets(ptr, Rf_mkString("QuickBLAST_XPtr"));
   unsigned int list_size = cppObj_list.size();
   
   cppObj_list.insert(std::make_pair(list_size, ptr));
   
   return ptr;
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("CreateQuickBLASTInstance() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("CreateQuickBLASTInstance(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 } catch(const std::exception &e){
   Rcpp::Rcerr << std::string("CreateQuickBLASTInstance() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << std::string("CreateQuickBLASTInstance() - Unknown Exception") << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}

//' @name GetInstanceCount
//' @title Get count of QuickBLAST instances stored in C++ side
//'
//' @description This function gives the size of the list of QuickBLAST C++ object list
//'
//' @return Count of QuickBLAST instances
//' @examples
//' \donttest{
//' QuickBLAST::GetInstanceCount()
//' }
//' @export
// [[Rcpp::export]]
unsigned int GetInstanceCount()
{
 return (unsigned int)cppObj_list.size();
}

Rcpp::XPtr<QuickBLAST> ResolveQuickBLASTInstance(SEXP inst)
{
  switch (TYPEOF(inst)) {
  case INTSXP: 
  case REALSXP: {
    // numeric 
    unsigned int ptr_id = Rcpp::as<unsigned int>(inst);
    try {
      return cppObj_list.at(ptr_id);
    } catch (const std::out_of_range& e) {
      Rcpp::stop(std::string("ResolveQuickBLASTInstance() - Error fetching QuickBLAST instance (") + std::to_string(ptr_id) + std::string("): ") + e.what());
    }
    break;
  }
  case EXTPTRSXP: {
    // 1. Check if it has the correct class
    if (!Rf_inherits(inst, "QuickBLAST_XPtr")) {
    Rcpp::stop("ResolveQuickBLASTInstance(): Expected external pointer of class 'QuickBLAST_XPtr'");
  }
    
    // 2. Check if the internal C++ pointer is NULL (dead pointer from loaded workspace)
    if (R_ExternalPtrAddr(inst) == NULL) {
      Rcpp::Rcout << "Dead pointer detected. Attempting to reload..." << std::endl;
      
      RObject rob(inst);
      
      // Safely check if required attributes exist before trying to read them
      if (!rob.hasAttribute("seq_type") || !rob.hasAttribute("program")) {
        Rcpp::stop("Cannot reload QuickBLAST pointer: Missing required attributes in saved R object.");
      }
      
      int seq_type           = Rcpp::as<int>(rob.attr("seq_type"));
      int strand             = Rcpp::as<int>(rob.attr("strand"));
      bool save_sequences    = Rcpp::as<bool>(rob.attr("save_sequences"));
      bool save_hsp_sequences= Rcpp::as<bool>(rob.attr("save_hsp_sequences"));
      SEXP s_prog            = rob.attr("program");
      
      SEXP s_opts = R_NilValue;
      if (rob.hasAttribute("options")) {
        s_opts = rob.attr("options");
      }
      
      // Re-instantiate
      SEXP newinst = CreateQuickBLASTInstance(seq_type, strand, s_prog, s_opts,
                                              save_sequences, save_hsp_sequences);
      
      // Update the original SEXP so R knows about the new memory address
      R_SetExternalPtrAddr(inst, R_ExternalPtrAddr(newinst));
    }
    
    // 3. Return the guaranteed-valid pointer
    return Rcpp::XPtr<QuickBLAST>(inst);
  }
  case NILSXP: {
    Rcpp::stop("ResolveQuickBLASTInstance(): pointer is NULL");
    break;
  }
  default:
    Rcpp::stop(std::string("Unsupported SEXP type: ") + Rf_type2char(TYPEOF(inst)));
  }
}

//' @name GetInstanceID
//' @title Get ID/Index of a QuickBLAST instance stored in C++ side
//'
//' @description This function fetches the ID/Index of a QuickBLAST instance of a \code{Rcpp::XPtr<QuickBLAST>} stored in C++ side.
//'
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//'
//' @return (unsigned int) ID/Index of the QuickBLAST instance pointer, FALSE otherwise
//' @examples
//' \donttest{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F,
//'   num_threads=24
//' )
//' QuickBLAST::GetInstanceID(
//'   blastp_inst
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP GetInstanceID(SEXP ptr)
{
 try
 {
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   // Find the first key associated with the target value
   auto it = std::find_if(cppObj_list.begin(), cppObj_list.end(),
                          [&](const std::pair<unsigned int, Rcpp::XPtr<QuickBLAST>>& pair) {
                            return pair.second.get() == ptr_.get();
                          });
   
   if (it != cppObj_list.end()) {
     return Rcpp::wrap(it->first);
   }
   return Rcpp::wrap(false);
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("GetInstanceID() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("GetInstanceID(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("GetInstanceID() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "GetInstanceID() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}

//' @name GetQuickBLASTInstance
//' @title Get QuickBLAST instance stored in C++ side at ID/Index
//'
//' @description This function fetches the QuickBLAST instance of a \code{Rcpp::XPtr<QuickBLAST>} at ID/Index stored in C++ side.
//'
//' @param ptr_id (unsigned int) ID/Index of Pointer to a QuickBLAST instance (in C++ side).
//' @examples
//' \donttest{
//' # QuickBLAST::GetQuickBLASTInstance(0)
//' }
//' @return (\code{Rcpp::XPtr<QuickBLAST>}) Pointer to a QuickBLAST instance, FALSE otherwise
//' @export
// [[Rcpp::export]]
SEXP GetQuickBLASTInstance(unsigned int ptr_id)
{
 try {
   return cppObj_list.at(ptr_id);
 } catch (const std::out_of_range& e) {
   Rcpp::Rcerr << "Error fetching QuickBLAST instance (" << ptr_id << "): " << e.what() << std::endl;
 }
 return Rcpp::wrap(false);
}

unsigned int DetectThreadLimit(unsigned int num_threads){
  unsigned int hw = std::thread::hardware_concurrency();
  if (hw == 0) hw = 1;                      // fallback if implementation returns 0
  unsigned int threads = num_threads > 0 ? num_threads > hw ? hw : num_threads : hw;
  
  // optional: respect OMP_NUM_THREADS if user set it and num_threads==0
  if (num_threads == 0) {
    const char* omp_env = std::getenv("OMP_NUM_THREADS");
    if (omp_env) {
      unsigned long envv = std::stoul(omp_env);
      if (envv > 0) threads = static_cast<unsigned int>(envv);
    }
  }
  
  return threads;
}

//' @name DeleteQuickBLASTInstance
//' @title Delete a QuickBLAST instance stored in C++ side
//'
//' @description This function deletes a QuickBLAST instance based on the instance ID
//'
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//'
//' @return TRUE - if the instance is deleted successfully, throws error otherwise
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::DeleteQuickBLASTInstance(
//'   QuickBLAST::GetInstanceID(
//'     blastp_inst
//'   )
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP DeleteQuickBLASTInstance(SEXP ptr)
{
 try
 {
   if(ptr == R_NilValue){ 
     Rcpp::Rcerr << "DeleteQuickBLASTInstance(): Input pointer/ID cannot be NULL." << std::endl << std::flush;
     return Rcpp::wrap(false); 
   }
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   auto ptr_id = Rcpp::as<unsigned int>(GetInstanceID(ptr));
   ptr_.release();
   
   cppObj_list.erase(ptr_id);
   
   return Rcpp::wrap(true);
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("DeleteQuickBLASTInstance() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("DeleteQuickBLASTInstance(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("DeleteQuickBLASTInstance() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "DeleteQuickBLASTInstance() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}

//' @name GetQuickBLASTOptions
//'
//' @title Get BLAST options for a QuickBLAST instance as a string.
//'
//' @description Get the BLAST options for a QuickBLAST instance.
//'
//' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()], [QuickBLAST::SetQuickBLASTOptions()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @return (string) BLAST options as std::string
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "tblastn",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::GetQuickBLASTOptions(
//'   blastp_inst
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP GetQuickBLASTOptions(SEXP ptr)
{
 try
 {
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   return Rcpp::wrap(ptr_->GetQuickBLASTOptionString());
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("GetQuickBLASTOptions() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("GetQuickBLASTOptions(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("GetQuickBLASTOptions() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "GetQuickBLASTOptions() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}

//' @name SetQuickBLASTOptions
//'
//' @title Set BLAST options for a QuickBLAST instance.
//'
//' @description Set/Modify the BLAST options for a QuickBLAST instance.
//'
//' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param program_name (string) Name of the BLAST program
//' @param options (string (or) Named List) List of BLAST options - check QuickBLAST::GetAvailableBLASTOptions(). String should be of the format "-option1 value1 -option2 value2"
//' @param verbose (bool) Verbose?
//' @return (bool) TRUE - if options set for the QuickBLAST instance, FALSE otherwise
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "tblastn",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::SetQuickBLASTOptions(
//'   blastp_inst,
//'   "blastp",
//'   "-evalue 1"
//' )
//' }
//' @export
// [[Rcpp::export]]
bool SetQuickBLASTOptions(SEXP ptr, SEXP program_name, SEXP options, bool verbose = true)
{
 try
 {
   std::string program_name_ = Rcpp::as<std::string>(program_name);
   std::string options_ = ConvertBLASTOptions2String_std(options);
   
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   ptr_->GetQuickBLASTOptions() = *ptr_->SetQuickBLASTOptions(program_name_, options_, CBlastOptions::EAPILocality::eLocal, verbose);
   return true;
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("SetQuickBLASTOptions() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("SetQuickBLASTOptions(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("SetQuickBLASTOptions() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "SetQuickBLASTOptions() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}

//TODO: Print Alignments and write an overload the infix %BLAST% operator to call BLAST2Seqs()

//' @name BLAST2Seqs
//'
//' @title BLAST 2 Sequence strings
//'
//' @description BLAST 2 nucleotide or protein strings with a QuickBLAST instance.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query sequence
//' @param subject (string) Subject sequence.
//' @param verbose (bool) Verbosity (Default: TRUE).
//' @return (\code{Rcpp::XPtr<QuickBLAST>}) Pointer to a QuickBLAST Instance (Cannot be used in R)
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::BLAST2Seqs(
//'   blastp_inst,
//'   "MQILLVEDDNTLFQELKKELEQWDFNV
//'   AGIEDFGKVMDTFESFNPEIVILDVQLP
//'   KYDGFYWCRKMREVSNVPILFLSSRDNP
//'   MDQVMSMELGADDYMQKPFYTNVLIAKL
//'   QAIYRRVYEFTAEEKRTLTWQDAVVDLS
//'   KDSIQKGDDTIFLSKTEMIILEILITKK
//'   NQIVSRDTIITALWDDEAFVSDNTLTVN
//'   VNRLRKKLSEISMDSAIETKVGKGYMAHE",
//'   "MQILLVEDDNTLFQELKKELEQWDFNV
//'   AGIEDFGKVMDTFESFNPEIVILDVQLP
//'   KYDGFYWCRKMREVSNVPILFLSSRDNP
//'   MDQVMSMELGADDYMQKPFYTNVLIAKL
//'   QAIYRRVYEFTAEEKRTLTWQDAVVDLS
//'   KDSIQKGDDTIFLSKTEMIILEILITKK
//'   NQIVSRDTIITALWDDEAFVSDNTLTVN
//'   VNRLRKKLSEISMDSAIETKVGKGYMAHE"
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP BLAST2Seqs(SEXP ptr, SEXP query, SEXP subject, bool verbose = true)
{
 try{
   
   if (query == R_NilValue || TYPEOF(query) != STRSXP || Rf_length(query) != 1) {
     Rcpp::Rcerr << "query must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject == R_NilValue || TYPEOF(subject) != STRSXP || Rf_length(subject) != 1) {
     Rcpp::Rcerr << "subject must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   auto start = std::chrono::high_resolution_clock::now();
   
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   std::string query_ = Rcpp::as<std::string>(query);
   std::string subject_ = Rcpp::as<std::string>(subject);
   
   std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
   
   std::shared_ptr<arrow::RecordBatch> ret_rb = ptr_->BLAST_seqs(query_, subject_, verbose); //
   
   if (ret_rb)
   {
     try{
       arrow::Status rb_sts = ret_rb->ValidateFull();
       
       if (!rb_sts.ok())
       {
         Rcpp::Rcerr << "[BLAST2Seqs()] ERR : Invalid RB :" << rb_sts.message().c_str() << std::endl << rb_sts.detail()->ToString().c_str() << std::endl <<std::flush;
       }
       else
       {
         if(ret_rb->num_rows() <= 0){
           Rcpp::Rcerr << "[BLAST2Seqs()] No alignments could be computed" << std::endl <<std::flush;
         }
         ret_val->emplace_back(ret_rb); 
       }
     }catch(...){
       // ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
       Rcpp::Rcerr << "[BLAST2Seqs()] Unknown Error" << std::endl <<std::flush;
       return Rcpp::wrap(false); 
     }
   }
   
   if(verbose)
     PrintClock(start);
   
   // Rcpp::List ret_vals_ = as<List>(Hits2RList(*ret_val));
   Rcpp::Rcout << ret_val->size() << std::endl << std::flush; //DEBUG
   if(ret_val->size() > 0){
     auto ret_vals_ = RecordBatchVectorToFlattenedDFList_parallel(*ret_val);
     ret_val->clear();
     ret_val->shrink_to_fit();
     return ret_vals_; //rm_null(ret_vals_);
   }else{
     // Rcpp::Function invisible("invisible");
     // return invisible(R_NilValue);
     return Rcpp::wrap(false);
   }
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Seqs() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLAST2Seqs(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Seqs() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLAST2Seqs() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}



//' @name BLAST2Folders
//'
//' @title BLAST 2 Folders
//'
//' @description BLAST 2 Folders with FASTA files containing nucleotide or protein sequences with a QuickBLAST instance. The files from query and subject folders are selected with the extension parameter
//'
//' @note Only FASTA files are supported by this function, use [QuickBLAST::BLAST2DBs()] if inputs are BLAST DBs.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query folder
//' @param subject (string) Subject folder.
//' @param extension (string) Extension of files in folder.
//' @param out_folder (string) Ouput Folder (Required).
//' @param out_format (string) Ouput Format. 'ipc'/'csv'/'parquet' (Optional) (Default: 'parquet').
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param reciprocal_hits (bool) Perform Bi-directional (Reciprocal => query <-> subject) BLAST? (Default: FALSE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size - Size of file write buffer (Optional).
//' @param verbose (bool) Verbosity (Defaut: TRUE).
//' @return (bool) TRUE - on success, FALSE - Otherwise. (Results are not returned as R Lists to reduce overhead)
//' @export
// [[Rcpp::export]]
bool BLAST2Folders(SEXP ptr, SEXP query, SEXP subject, SEXP extension, SEXP out_folder, SEXP out_format = R_NilValue, unsigned int num_threads = 0, bool reciprocal_hits = false, unsigned int min_batch_size = 0, bool verbose = true)
{
 try{
   auto start = std::chrono::high_resolution_clock::now();
   unsigned int threads = DetectThreadLimit(num_threads);
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   if (query == R_NilValue || TYPEOF(query) != STRSXP || Rf_length(query) != 1) {
     Rcpp::Rcerr << "query must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject == R_NilValue || TYPEOF(subject) != STRSXP || Rf_length(subject) != 1) {
     Rcpp::Rcerr << "subject must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (extension == R_NilValue || TYPEOF(extension) != STRSXP || Rf_length(extension) != 1) {
     Rcpp::Rcerr << "extension must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::string query_ = Rcpp::as<std::string>(query);
   std::string subject_ = Rcpp::as<std::string>(subject);
   std::string out_folder_ = Rcpp::as<std::string>(out_folder);
   std::string out_format_ = out_format == R_NilValue ? "parquet" : Rcpp::as<std::string>(out_format);
   std::string extension_ = Rcpp::as<std::string>(extension);
   
   std::filesystem::path outPath(out_folder_);
   std::filesystem::create_directory(outPath);
   if (!std::filesystem::is_empty(outPath))
   {
     Rcpp::Rcerr << std::string("[BLAST2Folders()] out_folder : Folder must be empty -") + out_folder_ << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::filesystem::path query_path(query_);
   std::filesystem::path subject_path(subject_);
   std::vector<std::string> queryFiles = getFilesInDir_std(query_, extension_);
   std::vector<std::string> subjectFiles = getFilesInDir_std(subject_, extension_);
   
   if(queryFiles.empty()){
     Rcpp::Rcerr << "query directory is empty"<< std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if(subjectFiles.empty()){
     Rcpp::Rcerr << "subject directory is empty" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   int iterations = queryFiles.size() * subjectFiles.size();
   
   Progress progress_bar(iterations, verbose);
   
   for (int i = 0; i < (int)queryFiles.size(); i++)
   {
     for (int j = 0; j < (int)subjectFiles.size(); j++)
     {
       if(Progress::check_abort()){
         Rcpp::Rcerr << "Aborted by Progress::check_abort()." << std::endl <<std::flush;
         return Rcpp::wrap(false); 
       }
       RcppThread::checkUserInterrupt();
       progress_bar.increment();
       if (queryFiles[i] == subjectFiles[j])
       {
         continue; //Same file, skipping
       }
       
       std::string base_name = getFilenameWithoutExtension(queryFiles[i]) + "-" + getFilenameWithoutExtension(subjectFiles[j]);
       std::filesystem::path outFile_ = outPath / (base_name + ".hits");
       std::string outFileStr = outFile_.string();
       std::filesystem::path query_input = query_path / queryFiles[i];
       std::filesystem::path subject_input = subject_path / subjectFiles[j];
       
       if (!std::filesystem::exists(query_input.string()) || !std::filesystem::exists(subject_input.string()))
       {
         Rcpp::Rcerr << "[BLAST2Folders()] Warn : File not found : " << query_input.string() << " or " << subject_input.string() << std::endl;
         continue;
       }
       
       if (!reciprocal_hits)
       {
         std::string base_name_rbh = getFilenameWithoutExtension(subjectFiles[j]) + "-" + getFilenameWithoutExtension(queryFiles[i]);
         std::filesystem::path outFile_rbh = outPath / (base_name_rbh + ".hits");
         
         if (std::filesystem::exists(outFile_rbh.string()))
         {
           continue;
         }
       }
       int one_dim_index = i + queryFiles.size() * j;
       if(verbose)
         Rcpp::Rcout << "BLASTing : " << base_name << std::endl;
       static_cast<void>(ptr_->BLAST_files(query_input.string(), subject_input.string(), outFileStr, out_format_, threads, /*return_values*/ false, min_batch_size, verbose));
       
     }
   }
   
   if(verbose)
     PrintClock(start);
   
   return Rcpp::wrap(true);
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Folders() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLAST2Folders(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Folders() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLAST2Folders() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}

//' @name BLAST1Folder
//'
//' @title BLAST Files within a Folder
//'
//' @description BLAST FASTA files containing nucleotide or protein sequences, within a folder with a QuickBLAST instance. The files from the folder are selected with the extension parameter and BLAST'd against each other.
//'
//' @note Only FASTA files are supported by this function, use [QuickBLAST::BLAST2DBs()] if inputs are BLAST DBs.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param input_folder (string) Input folder
//' @param extension (string) Extension of files in folder.
//' @param out_folder (string) Ouput Folder (Required).
//' @param out_format (string) Ouput Format. 'ipc'/'csv'/'parquet' (Optional) (Default: 'parquet').
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param reciprocal_hits (bool) Perform Bi-directional (Reciprocal => query <-> subject) BLAST? (Default: FALSE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size - Size of file write buffer (Optional).
//' @param verbose (bool) Verbosity (Defulat: TRUE).
//' @return (bool) TRUE - on success, FALSE - Otherwise. (Results are not returned as R Lists to reduce overhead)
//' @export
// [[Rcpp::export]]
bool BLAST1Folder(SEXP ptr, SEXP input_folder, SEXP extension, SEXP out_folder, SEXP out_format = R_NilValue, unsigned int num_threads = 0, bool reciprocal_hits = false, unsigned int min_batch_size = 0, bool verbose = true)
{
 
 try{
   auto start = std::chrono::high_resolution_clock::now();
   unsigned int threads = DetectThreadLimit(num_threads);
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   if (out_folder == R_NilValue || TYPEOF(out_folder) != STRSXP || Rf_length(out_folder) != 1) {
     Rcpp::Rcerr << "out_folder must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (input_folder == R_NilValue || TYPEOF(input_folder) != STRSXP || Rf_length(input_folder) != 1) {
     Rcpp::Rcerr << "input_folder must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (extension == R_NilValue || TYPEOF(extension) != STRSXP || Rf_length(extension) != 1) {
     Rcpp::Rcerr << "extension must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::string input_folder_ = Rcpp::as<std::string>(input_folder);
   std::string out_folder_ = Rcpp::as<std::string>(out_folder);
   std::string extension_ = Rcpp::as<std::string>(extension);
   std::string out_format_ = out_format == R_NilValue ? "parquet" : Rcpp::as<std::string>(out_format);
   
   if(input_folder_.empty() || out_folder_.empty() || out_format_.empty()){
     Rcpp::Rcerr << "Input cannot be empty: input_folder || out_folder || out_format" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::filesystem::path outPath(out_folder_);
   std::filesystem::create_directory(outPath);
   if (!std::filesystem::is_empty(outPath))
   {
     Rcpp::Rcerr << std::string("[BLAST1Folder()] out_folder : Folder must be empty - ") + out_folder_ << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::vector<std::string> inputFiles = getFilesInDir_std(input_folder_, extension_);
   
   if(inputFiles.empty()){
     Rcpp::Rcerr << "input directory is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::filesystem::path folder = input_folder_;
   
   int iterations = inputFiles.size() * inputFiles.size();
   Progress progress_bar(iterations, verbose);
   
   for (int i = 0; i < (int)inputFiles.size(); i++)
   {
     for (int j = 0; j < (int)inputFiles.size(); j++)
     {
       if(Progress::check_abort()){
         Rcpp::Rcerr << "Aborted by Progress::check_abort()." << std::endl <<std::flush;
         return Rcpp::wrap(false); 
       }
       RcppThread::checkUserInterrupt();
       progress_bar.increment();
       if (inputFiles[i] == inputFiles[j])
       {
         continue; //Same file, skipping
       }
       
       std::string base_name = getFilenameWithoutExtension(inputFiles[i]) + "-" + getFilenameWithoutExtension(inputFiles[j]);
       std::filesystem::path outFile_ = outPath / (base_name + ".hits");
       std::string outFileStr = outFile_.string();
       
       if (!reciprocal_hits)
       {
         std::string base_name_rbh = getFilenameWithoutExtension(inputFiles[j]) + "-" + getFilenameWithoutExtension(inputFiles[i]);
         std::filesystem::path outFile_rbh = outPath / (base_name_rbh + ".hits");
         
         if (std::filesystem::exists(outFile_rbh.string()))
         {
           continue;
         }
       }
       int one_dim_index = i + inputFiles.size() * j;
       std::filesystem::path qry_filename = inputFiles[i];
       std::filesystem::path qry_filePath = folder / qry_filename;
       std::filesystem::path subj_filename = inputFiles[j];
       std::filesystem::path subj_filePath = folder / subj_filename;
       if(verbose)
         Rcpp::Rcout << "BLASTing : " << base_name << std::endl;
       static_cast<void>(ptr_->BLAST_files(qry_filePath.string(), subj_filePath.string(), outFileStr, out_format_, threads, /*return_values*/ false, min_batch_size, verbose));
       
     }
   }
   
   if(verbose)
     PrintClock(start);
   
   return Rcpp::wrap(true);
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLAST1Folder() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLAST1Folder(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLAST1Folder() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLAST1Folder() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}


//' @name BLAST2Files
//'
//' @title BLAST 2 Files
//'
//' @description BLAST 2 FASTA files containing nucleotide or protein sequences with a QuickBLAST instance.
//'
//' @note Only FASTA files are supported by this function, use [QuickBLAST::BLAST2DBs()] if inputs are BLAST DBs.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2DBs()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query file
//' @param subject (string) Subject file
//' @param out_file (string) Ouput file (Optional)
//' @param out_format (string) Ouput Format. 'ipc'/'csv'/'parquet' (Optional) (Default: 'parquet').
// ' @param seq_limit (int) Batch Size to BLAST at a time. { -1 = Whole File, 0 - One sequence at a time or > 0 } - Size of BLAST sequences buffer (Optional)
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size - Size of file write buffer (Optional).
//' @param verbose (bool) Verbosity (Default: TRUE).
//' @return (SEXP) Rcpp::List - if return_values == TRUE, out_file - Otherwise.
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::BLAST2Files(
//'   ptr = blastp_inst,
//'   query = system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   subject = system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   out_file = "test.arrow",
//'   out_format = "parquet",
//'   return_values = F,
//'   min_batch_size = 1024
//' )
//' QuickBLAST::BLAST2Files(
//'   ptr = blastp_inst,
//'   query = system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   subject = system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   out_file = "test.arrow",
//'   return_values = T,
//'   min_batch_size = 0,
//'   seq_limit = 0
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP BLAST2Files(SEXP ptr, SEXP query, SEXP subject, SEXP out_file = R_NilValue, SEXP out_format = R_NilValue, unsigned int num_threads = 0, bool return_values = true, unsigned int min_batch_size = 0, bool verbose = true)
{
 try{
   auto start = std::chrono::high_resolution_clock::now();
   unsigned int threads = DetectThreadLimit(num_threads);
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   if(out_file == R_NilValue && return_values == false){
     Rcpp::Rcerr << "Error: Both out_file cannot be empty and return_values == FALSE" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (query == R_NilValue || TYPEOF(query) != STRSXP || Rf_length(query) != 1) {
     Rcpp::Rcerr << "Query must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject == R_NilValue || TYPEOF(subject) != STRSXP || Rf_length(subject) != 1) {
     Rcpp::Rcerr << "subject must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::string query_ = Rcpp::as<std::string>(query);
   std::string subject_ = Rcpp::as<std::string>(subject);
   std::string out_file_ = out_file == R_NilValue ? "" : Rcpp::as<std::string>(out_file);
   std::string out_format_ = out_format == R_NilValue ? "parquet" : Rcpp::as<std::string>(out_format);   
   if (query_.empty()) {
     Rcpp::Rcerr << "query is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject_.empty()) {
     Rcpp::Rcerr << "subject is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (out_format_.empty()) {
     Rcpp::Rcerr << "out_format is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (!std::filesystem::exists(query_)) {
     Rcpp::Rcerr << "[BLAST2Files()] Error: Query file does not exist or is not accessible: " << query_ << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   if (!std::filesystem::exists(subject_)) {
     Rcpp::Rcerr << "[BLAST2Files()] Error: Subject file does not exist or is not accessible: " << subject_ << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   
   std::shared_ptr<arrow::RecordBatchVector> ret_vals = std::make_shared<arrow::RecordBatchVector>();
   
   if (return_values)
   {
     ret_vals = ptr_->BLAST_files(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size, verbose); 
   }
   else
   {
     static_cast<void>(ptr_->BLAST_files(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size, verbose));
   }
   
   if(verbose)
     PrintClock(start);
   
   if(return_values){
     if(ret_vals->size() > 0){
       auto ret_vals_ = RecordBatchVectorToFlattenedDFList_parallel(*ret_vals);
       ret_vals->clear();
       ret_vals->shrink_to_fit();
       return ret_vals_; //rm_null(ret_vals_);
     }else{
       return Rcpp::wrap(false);
     }
   }else{
     return Rcpp::wrap(true);
   }
   
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Files() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLAST2Files(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLAST2Files() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLAST2Files() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}


//' @name RemoteBLAST
//'
//' @title BLAST query against remote NCBI DBs
//'
//' @description BLAST the input query against remote NCBI DBs (one sequence at a time - to respect rate limits)
//'
//' @note Check BLAST Guide (\url{https://blast.ncbi.nlm.nih.gov/BLAST_guide.pdf}) and NCBI BLAST (\url{https://blast.ncbi.nlm.nih.gov/Blast.cgi}) (Program -> Choose DB/Search set) for database names.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param database (string) Name of the remote NCBI DB - Check note for reference and supported values.
//' @param query_input (Rcpp::List) (Named) List of input queries (Sequences, Files, Folders - type is determined by input_type parameter)
//' @param input_type (QuickBLAST::EInputType) Input type (Check [QuickBLAST::GetQuickBLASTEnums()])
//' @param outFile (string) Output file name (Optional)
//' @param outFormat (string) Format of Output File (Required for outFile) (Default: 'parquet')
//' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
//' @param max_poll_seconds (int) Max seconds to wait for RemoteBLAST (Default: 360) (Optional)
//' @param poll_interval_ms (int) Milliseconds wait-time between polling RemoteBLAST service (Default(4s): 4000) (Optional)
//' @param verbose (bool) Verbosity (Default: TRUE).
//' @return (SEXP) Rcpp::List - if return_values == TRUE, outFile - Otherwise.
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::RemoteBLAST(
//'   blastp_inst,
//'   query_input="MQILLVEDDNTLFQELKKELEQWDFN
//'   VAGIEDFGKVMDTFESFNPEIVILDVQLPKYDGFYWCRK
//'   MREVSNVPILFLSSRDNPMDQVMSMELGADDYMQKPFYT
//'   NVLIAKLQAIYRRVYEFTAEEKRTLTWQDAVVDLSKDSI
//'   QKGDDTIFLSKTEMIILEILITKKNQIVSRDTIITALWD
//'   DEAFVSDNTLTVNVNRLRKKLSEISMDSAIETKVGKGYMAHE",
//'   database= "pdb",
//'   input_type=1,
//'   return_values=T
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP RemoteBLAST(SEXP ptr, SEXP database, SEXP query_input, int input_type, SEXP outFile = R_NilValue, SEXP outFormat = R_NilValue, bool return_values = true, unsigned int max_poll_seconds = 360, unsigned int poll_interval_ms = 4000, bool verbose = true)
{
 try{
   auto start = std::chrono::high_resolution_clock::now();
   
   if (outFile == R_NilValue && return_values == false) {
     Rcpp::Rcerr << "Error: Both outFile cannot be empty and return_values == FALSE" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (query_input == R_NilValue || TYPEOF(query_input) != STRSXP || Rf_length(query_input) != 1) {
     Rcpp::Rcerr << "query_input must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   std::string program_ = ptr_->GetProgram();
   std::string database_ = Rcpp::as<std::string>(database);
   QuickBLAST::EInputType input_type_ = static_cast<QuickBLAST::EInputType>(input_type);
   std::string outFile_ = outFile == R_NilValue ? "" : Rcpp::as<std::string>(outFile); 
   std::string outFormat_ = outFormat == R_NilValue ? "parquet" : Rcpp::as<std::string>(outFormat); 
   std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
   
   switch(input_type_){
   case QuickBLAST::EInputType::eFile: {
     Rcpp::Rcerr << "Only QuickBLAST::EInputType::eSequenceString is implemented." << std::endl <<std::flush;
     break;
   }
   case QuickBLAST::EInputType::eSequenceString:{
     Rcpp::List query_input_ = Rcpp::as<Rcpp::List>(query_input);
     if(return_values){
       ret_val = ptr_->BLAST_remote(program_, database_, query_input_, input_type_, outFile_, outFormat_, return_values, max_poll_seconds, poll_interval_ms, verbose);
     }else{
       static_cast<void>(ptr_->BLAST_remote(program_, database_, query_input_, input_type_, outFile_, outFormat_, return_values, max_poll_seconds, poll_interval_ms, verbose));
       // ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
     }
     break;
   }
   case QuickBLAST::EInputType::eFolder:{
     Rcpp::Rcerr << "Only QuickBLAST::EInputType::eSequenceString is implemented." << std::endl <<std::flush;
     break;
   }
   case QuickBLAST::EInputType::eBLASTDB:{
     Rcpp::Rcerr << "Only QuickBLAST::EInputType::eSequenceString is implemented." << std::endl <<std::flush;
     break;
   }
   }
   if(verbose)
     PrintClock(start);
   if(return_values){
     if(ret_val->size() > 0){
       auto ret_vals_ = RecordBatchVectorToFlattenedDFList_parallel(*ret_val);
       ret_val->clear();
       ret_val->shrink_to_fit();
       return ret_vals_; //rm_null(ret_vals_);
     }else{
       return Rcpp::wrap(false);
     }
   }else{
     return Rcpp::wrap(true);
   }
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("RemoteBLAST() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("RemoteBLAST(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("RemoteBLAST() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "RemoteBLAST() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}

//' @name isBLASTDB
//'
//' @title Check BLAST DB Files
//'
//' @description Check whether a path/name corresponds to a BLAST DB (heuristic). Check if all the files of a BLAST DB exist (.psq .pog .pin .phr .pos .pto .pot .pdb .ptf .pjs). This checks the directory for filenames that start with the provided basename and have extensions commonly produced by makeblastdb:
//'  - protein:  .phr .pin .psq
//'  - nucleotide: .nhr .nin .nsq
//' It returns a list with a boolean `is_db`, a guessed `type` ("Protein"/"Nucleotide"/"Mixed"/"Unknown"),
//' a character vector of matching `files`, and the `dir` and `name` used.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::MakeBLASTDB()], [QuickBLAST::BLAST2DBs()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param input_db character(1) path to db (path + name) or a bare name (current directory assumed)
//' @return list with keys: is_db (logical), type (string), files (character vector), dir (string), name (string), message (string)
//' @examples
//' \dontrun{
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   "protein_query.db"
//' )
//' QuickBLAST::isBLASTDB(
//'   tools::file_path_sans_ext(
//'     system.file(
//'       "extdata",
//'       "protein_query.db.pin",
//'       package = "QuickBLAST",
//'       mustWork = T
//'     )
//'   )
//' )
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List isBLASTDB(SEXP ptr, SEXP input_db) {
 auto return_false_list = [](std::string msg) {
   return Rcpp::List::create(
     Rcpp::_["is_db"] = false,
     Rcpp::_["type"] = "Unknown",
     Rcpp::_["files"] = Rcpp::StringVector::create(),
     Rcpp::_["dir"] = "",
     Rcpp::_["name"] = "",
     Rcpp::_["db_name"] = "",
     Rcpp::_["seq_count"] = 0,
     Rcpp::_["message"] = msg
   );
 };
 
 try {
   
   if (TYPEOF(input_db) != STRSXP) {
     Rcpp::Rcerr << "isBLASTDB(): input must be a string" << std::endl <<std::flush;
     // return Rcpp::wrap(false); 
     return return_false_list("isBLASTDB(): input must be a string");
   }
   std::string input_db_ = Rcpp::as<std::string>(input_db);
   if (input_db_.empty()) {
     Rcpp::Rcerr << "isBLASTDB(): input string is empty" << std::endl <<std::flush;
     // return Rcpp::wrap(false); 
     return return_false_list("isBLASTDB(): input string is empty");
   }
   
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   std::filesystem::path dbPath(input_db_);
   std::filesystem::path dbDir;
   std::string dbName;
   std::string dbActualName = dbPath.string();
   
   std::error_code ec_cwd; // For safely getting current path
   
   if (dbPath.has_filename()) {
     dbDir = dbPath.parent_path();
     dbName = dbPath.stem().string();
     if (dbDir.empty()) {
       dbDir = std::filesystem::current_path(ec_cwd);
       if (ec_cwd) return return_false_list("Failed to get current_path");
     }
   } else {
     dbDir = dbPath;
     if (dbDir.empty()) {
       dbDir = std::filesystem::current_path(ec_cwd);
       if (ec_cwd) return return_false_list("Failed to get current_path");
     }
     Rcpp::Rcerr << "isBLASTDB(): input appears to be a directory; supply path+name or name (basename)\n";
     return return_false_list("isBLASTDB(): input appears to be a directory");
   }
   
   // SAFELY check if directory exists
   std::error_code ec_dir;
   if (!std::filesystem::exists(dbDir, ec_dir) || ec_dir || 
       !std::filesystem::is_directory(dbDir, ec_dir) || ec_dir) {
       std::string msg = "isBLASTDB(): directory does not exist or cannot be accessed: " + dbDir.string();
     Rcpp::Rcerr << msg << "\n";
     return return_false_list(msg);
   }
   
   // Known BLAST+ index extensions
   const std::vector<std::string> protein_exts = {".phr", ".pin", ".psq"};
   const std::vector<std::string> nucleotide_exts = {".nhr", ".nin", ".nsq"};
   const std::vector<std::string> extra_exts = {".pal", ".pal.idx", ".pal.rev"}; 
   
   std::vector<std::string> matched_files;
   bool matched_protein = false;
   bool matched_nucleotide = false;
   
   std::error_code ec_iter;
   std::filesystem::directory_options options = std::filesystem::directory_options::skip_permission_denied;
   
   // SAFELY initialize iterator
   std::filesystem::directory_iterator it(dbDir, options, ec_iter);
   if (ec_iter) {
     Rcpp::Rcerr << "[isBLASTDB()] Filesystem error creating iterator - " << ec_iter.message() << "\n";
     return return_false_list(std::string("[isBLASTDB()] Filesystem error - ") + ec_iter.message());
   }
   
   // SAFELY iterate using a while loop (DO NOT use range-based 'for')
   std::filesystem::directory_iterator end;
   while (it != end) {
     std::error_code ec_entry;
     
     // Safely check if it's a regular file
     if (it->is_regular_file(ec_entry) && !ec_entry) {
       std::string fname = it->path().filename().string();
       
       // check starts with basename
       if (fname.rfind(dbName, 0) == 0) {
         // record file
         matched_files.push_back(it->path().string());
         
         std::string ext = it->path().extension().string();
         std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
         
         auto protein_exts_check = std::find(protein_exts.begin(), protein_exts.end(), ext);
         auto nuc_exts_check = std::find(nucleotide_exts.begin(), nucleotide_exts.end(), ext);
         
         if (protein_exts_check != protein_exts.end()) matched_protein = true;
         if (nuc_exts_check != nucleotide_exts.end()) matched_nucleotide = true;
         
         // 1. Get the FULL absolute path, not just the filename
         std::string fullFilePath = it->path().string();
         
         // 2. Find the extension from the end of the string (remove the ', 0')
         auto extPos = fullFilePath.rfind(ext); 
         
         switch(ptr_->GetSeqType()){
         case QuickBLAST::ESeqType::eNucleotide: {
           if(matched_nucleotide)  
             if(extPos != std::string::npos) dbActualName = fullFilePath.substr(0, extPos);
             break;
         }
         case QuickBLAST::ESeqType::eProtein: {
           if(matched_protein)  
             if(extPos != std::string::npos) dbActualName = fullFilePath.substr(0, extPos);
             break;
         } 
         }
         
         // auto dbFile = it->path().filename().string();
         // auto extPos = dbFile.rfind(ext, 0);
         // 
         // switch(ptr_->GetSeqType()){
         // case QuickBLAST::ESeqType::eNucleotide: {
         //   if(matched_nucleotide)  
         //     if(extPos != std::string::npos) dbActualName = dbFile.substr(0, extPos);
         //     break;
         // }
         // case QuickBLAST::ESeqType::eProtein: {
         //   if(matched_protein)  
         //     if(extPos != std::string::npos) dbActualName = dbFile.substr(0, extPos);
         //     break;
         // } 
         // }
         
       }
     }
     
     // SAFELY advance the iterator
     std::error_code ec_inc;
     it.increment(ec_inc);
     if (ec_inc) {
       // If we hit an unrecoverable permission error moving to the next file, stop iterating
       break; 
     }
   }
   
   // std::filesystem::path dbPath(input_db_);
   // std::filesystem::path dbDir;
   // std::string dbName;
   // std::string dbActualName = dbPath.string();
   // 
   // if (dbPath.has_filename()) {
   //   dbDir = dbPath.parent_path();
   //   dbName = dbPath.stem().string();
   //   if (dbDir.empty()) dbDir = std::filesystem::current_path();
   // } else {
   //   // If no filename (rare), treat input as directory
   //   dbDir = dbPath;
   //   if (dbDir.empty()) dbDir = std::filesystem::current_path();
   //   Rcpp::Rcerr << "isBLASTDB(): input appears to be a directory; supply path+name or name (basename)" << std::endl <<std::flush;
   //   // return Rcpp::wrap(false);
   //   return return_false_list("isBLASTDB(): input appears to be a directory; supply path+name or name (basename)");
   // }
   // 
   // std::error_code ec_dir;
   // if (!std::filesystem::exists(dbDir, ec_dir) || ec_dir || 
   //     !std::filesystem::is_directory(dbDir, ec_dir) || ec_dir) {
   //   Rcpp::Rcerr << std::string("isBLASTDB(): directory does not exist: ") + dbDir.string() << std::endl <<std::flush;
   //   // return Rcpp::wrap(false); 
   //   return return_false_list(std::string("isBLASTDB(): directory does not exist: ") + dbDir.string());
   // }
   // 
   // // Known BLAST+ index extensions
   // const std::vector<std::string> protein_exts = {".phr", ".pin", ".psq"};
   // const std::vector<std::string> nucleotide_exts = {".nhr", ".nin", ".nsq"};
   // // Some older / alternate names sometimes seen:
   // const std::vector<std::string> extra_exts = {".pal", ".pal.idx", ".pal.rev"}; // optional
   // 
   // std::vector<std::string> matched_files;
   // bool matched_protein = false;
   // bool matched_nucleotide = false;
   // 
   // std::error_code ec;
   // // Tell the iterator to gracefully skip folders without permissions
   // std::filesystem::directory_options options = std::filesystem::directory_options::skip_permission_denied;
   // std::filesystem::directory_iterator it(dbDir, options, ec);
   // 
   // // Check if the base directory itself had an error (like permission denied)
   // if (ec) {
   //   Rcpp::Rcerr << "[isBLASTDB()] Filesystem error - " << ec.message() << std::endl << std::flush;
   //   return return_false_list(std::string("[isBLASTDB()] Filesystem error - ") + ec.message());
   // }
   // 
   // // Iterate directory and look for files that start with dbName
   // for (const auto &entry : it) {
   //   std::error_code ec_entry;
   //   // Skip if it throws a permission error OR if it's not a regular file
   //   if (!entry.is_regular_file(ec_entry) || ec_entry) continue;
   //   std::string fname = entry.path().filename().string();
   //   // check starts with basename (this handles name.ext and name.1.ext etc)
   //   if (fname.rfind(dbName, 0) != 0) continue;
   //   
   //   Rcpp::Rcout << "Matched File: " << entry.path().string() << std::endl << std::flush; //DEBUG
   //   // record file
   //   matched_files.push_back(entry.path().string());
   //   
   //   std::string ext = entry.path().extension().string();
   //   // standardize ext to lower-case
   //   std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
   //   
   //   auto protein_exts_check = std::find(protein_exts.begin(), protein_exts.end(), ext);
   //   auto nuc_exts_check = std::find(nucleotide_exts.begin(), nucleotide_exts.end(), ext);
   //   
   //   if (protein_exts_check != protein_exts.end()) {
   //     matched_protein = true;
   //   }
   //   if (nuc_exts_check != nucleotide_exts.end()) {
   //     matched_nucleotide = true;
   //   }
   //   // also consider extra_exts if needed
   //   if (std::find(extra_exts.begin(), extra_exts.end(), ext) != extra_exts.end()) {
   //     // don't toggle type here, just record file
   //   }
   //   
   //   auto dbFile = entry.path().filename().string();
   //   auto extPos = dbFile.rfind(ext, 0);
   //   
   //   switch(ptr_->GetSeqType()){
   //   case QuickBLAST::ESeqType::eNucleotide: {
   //     if(matched_nucleotide)  
   //       if(extPos != std::string_view::npos)
   //         dbActualName = dbFile.substr(0, extPos);
   //       break;
   //   }
   //   case QuickBLAST::ESeqType::eProtein: {
   //     if(matched_protein)  
   //       if(extPos != std::string_view::npos)
   //         dbActualName = dbFile.substr(0, extPos);
   //       break;
   //   } 
   //   }
   //   
   // }
   
   // Decide result
   // bool is_db = !matched_files.empty();
   bool is_db = (matched_protein || matched_nucleotide);
   std::string type = "Unknown";
   if (matched_protein && !matched_nucleotide) type = "Protein";
   else if (!matched_protein && matched_nucleotide) type = "Nucleotide";
   else if (matched_protein && matched_nucleotide) type = "Mixed";
   else if (is_db) type = "Unknown"; // found matching files but not the canonical ones
   
   std::string warnings = "None";
   switch(ptr_->GetSeqType()){
   case QuickBLAST::ESeqType::eNucleotide: {
     if (matched_protein && !matched_nucleotide) warnings = "isBLASTDB(): QuickBLAST Pointer is configured for eNucleotide but input DB is eProtein";
     else if (!matched_protein && matched_nucleotide) warnings = "None";
     else if (matched_protein && matched_nucleotide) warnings = "isBLASTDB(): QuickBLAST Pointer is configured for eNucleotide but input DB is matches eProtein & eNucleotide";
     break;
   }
   case QuickBLAST::ESeqType::eProtein: {
     if (matched_protein && !matched_nucleotide) warnings = "None";
     else if (!matched_protein && matched_nucleotide) warnings = "isBLASTDB(): QuickBLAST Pointer is configured for eProtein but input DB is eNucleotide";
     else if (matched_protein && matched_nucleotide) warnings = "isBLASTDB(): QuickBLAST Pointer is configured for eProtein but input DB is matches eProtein & eNucleotide";
     break;
   }
   }
   
   std::string message;
   if (!is_db) {
     message = "No matching BLAST DB files found for basename '" + dbName + "' in directory '" + dbDir.string() + "'";
   } else {
     message = "Found " + std::to_string(matched_files.size()) + " file(s) that match basename '" + dbName + "'";
   }
   
   // Build return list
   Rcpp::CharacterVector files(matched_files.size());
   for (size_t i = 0; i < matched_files.size(); ++i) files[i] = matched_files[i];
   
   return Rcpp::List::create(
     Rcpp::_["is_db"] = is_db,
     Rcpp::_["type"]  = type,
     Rcpp::_["files"] = files,
     Rcpp::_["dir"]   = dbDir.string(),
     Rcpp::_["name"]  = dbName,
     Rcpp::_["db_name"]  = dbActualName,
     Rcpp::_["seq_count"] =  matched_protein || matched_nucleotide ? ptr_->SizeOfDB(dbActualName) : 0,
     Rcpp::_["message"] = message,
     Rcpp::_["warnings"] = warnings
   );
 }catch (const Rcpp::exception &e) {
   Rcpp::Rcerr << std::string("isBLASTDB() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch (const std::filesystem::filesystem_error &e) {
   Rcpp::Rcerr << "isBLASTDB() - Filesystem Exception: " << e.what() << "\n";
 }catch (const std::runtime_error &e) {
   Rcpp::Rcerr << std::string("isBLASTDB() - C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch (const std::exception &e) {
   Rcpp::Rcerr << std::string("isBLASTDB() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch (...) {
   Rcpp::Rcerr << "isBLASTDB() - Unknown Exception" << std::endl <<std::flush;
 }
 // return Rcpp::wrap(false); 
 return return_false_list("[isBLASTDB()] Reached end of function.");
}

//@param stdout_opt (Bool/String) - Re-route STDOUT
//@param stderr_opt (Bool/String) - Re-route STDERR

//' @name MakeBLASTDB
//'
//' @title Make on-disk BLAST DB from a FASTA file
//'
//' @description Calls makeblastdb to create a BLAST DB of a FASTA file
//'
//' @note Faitful re-implementation of makeblastdb seemed pointless, hence the system.call() to a the program.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2DBs()], [QuickBLAST::isBLASTDB()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param input_file (string) Path to FASTA file.
//' @param database_name (string) Name of the output DB.
//' @param parse_seqids (bool) TRUE - Checks FASTA headers for malformations (Default: FALSE)
//' @return (bool) DB name on success, FALSE - Otherwise.
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   "protein_query.db"
//' )
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   "protein_subject.db"
//' )
//' QuickBLAST::BLAST2DBs(
//'   ptr=blastp_inst,
//'   query="protein_query.db",
//'   subject="protein_subject.db",
//'   num_threads=24,
//'   out_file="test.db.arrow",
//'   return_values = T
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP MakeBLASTDB(SEXP ptr, SEXP input_file, SEXP database_name, bool parse_seqids = false){
 try{
   
   if (TYPEOF(input_file) != STRSXP || input_file == R_NilValue) {
     Rcpp::Rcerr << "input_file must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false);
   }
   if (TYPEOF(database_name) != STRSXP || database_name == R_NilValue) {
     Rcpp::Rcerr << "database_name must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false);
   }
   
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   std::string input_file_ = Rcpp::as<std::string>(input_file);
   std::ifstream infile(input_file_);
   if (!infile.good()) {
     Rcpp::Rcerr << "MakeBLASTDB Error: " << input_file_ << " does not exist or cannot be read."  << std::endl <<std::flush;
     return Rcpp::wrap(false);
   }
   std::string database_name_ = Rcpp::as<std::string>(database_name);
   std::string parse_seqids_opt;
   if(parse_seqids){
     parse_seqids_opt = "-parse_seqids";
   }else{
     parse_seqids_opt = "";
   }
   
   Environment pkg = Environment::namespace_env("QuickBLAST");
   Function f = pkg[".GetBinPath"];
   std::string mbdb_exe = "makeblastdb";
   // --- 1. Handle Windows .exe Extension ---
#if defined(_WIN32) || defined(__MINGW32__)
   mbdb_exe += ".exe";
#endif
   std::string program_path = Rcpp::as<std::string>( f(Named("file_name")=mbdb_exe) );
   
   //      // --- 1. Handle Windows .exe Extension ---
   // #if defined(_WIN32) || defined(__MINGW32__)
   //      if (program_path.length() < 4 || program_path.substr(program_path.length() - 4) != ".exe") {
   //        program_path += ".exe";
   //      }
   // #endif
   
   // --- 2. Check if the program actually exists ---
   if (!FileExists(program_path)) {
     Rcpp::Rcerr << "[MakeBLASTDB] Error: Executable not found at path: "
                 << program_path << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   
   std::string dbtype;
   switch(ptr_->GetSeqType()){
   case QuickBLAST::ESeqType::eNucleotide: {
     dbtype = "nucl";
     break;
   }
   case QuickBLAST::ESeqType::eProtein: {
     dbtype = "prot";
     break;
   }
   }
   
   // 3. Build the base arguments
   std::vector<std::string> argv = {
     program_path,
     "-in", input_file_,
     "-dbtype", dbtype,
     "-out", database_name_
   };
   
   if (parse_seqids) {
     argv.push_back("-parse_seqids");
   }
   
   // 4. Print the exact command being run for Debugging Context
   Rcpp::Rcout << "[MakeBLASTDB] Executing command: ";
   for (const auto &s : argv) {
     Rcpp::Rcout << s << " ";
   }
   Rcpp::Rcout << std::endl << std::flush;
   
   // Build the C-style argv array
   std::vector<char*> cargv;
   for (const auto &s : argv) {
     cargv.push_back(const_cast<char*>(s.c_str()));
   }
   cargv.push_back(nullptr);
   
   
#if !defined(_WIN32) && !defined(__MINGW32__)
   // --- POSIX WAY (Linux / macOS) ---
   int pid;
   // posix_spawnp returns 0 on success, or an error code directly on failure
   int spawn_status = posix_spawnp(&pid, cargv[0], nullptr, nullptr, cargv.data(), environ);
   
   if (spawn_status == 0) {
     int wait_status;
     if (waitpid(pid, &wait_status, 0) == -1) {
       Rcpp::Rcerr << "[MakeBLASTDB] waitpid failed. errno: " << errno
                   << " (" << std::strerror(errno) << ")" << std::endl << std::flush;
       return Rcpp::wrap(false);
     }
     if (!WIFEXITED(wait_status) || WEXITSTATUS(wait_status) != 0) {
       Rcpp::Rcerr << "[MakeBLASTDB] process failed with non-zero exit status: "
                   << WEXITSTATUS(wait_status) << std::endl << std::flush;
       return Rcpp::wrap(false);
     }
   } else {
     // Use spawn_status for strerror on POSIX
     Rcpp::Rcerr << "[MakeBLASTDB] posix_spawnp failed to start makeblastdb. Error code: "
                 << spawn_status << " (" << std::strerror(spawn_status) << ")"
                 << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   
#else
   // --- WINDOWS WAY (MinGW / Rtools) ---
   int status = -1;
   intptr_t pid = _spawnvp(_P_NOWAIT, cargv[0], (char *const *)cargv.data());
   
   if (pid == -1) {
     // Use errno for strerror on Windows
     Rcpp::Rcerr << "[MakeBLASTDB] _spawnvp failed. Path: " << program_path
                 << " | errno: " << errno << " (" << std::strerror(errno) << ")"
                 << std::endl << std::flush;
     return Rcpp::wrap(false);
   } else {
     if (_cwait(&status, pid, WAIT_CHILD) == -1) {
       Rcpp::Rcerr << "[MakeBLASTDB] _cwait failed to track process. errno: "
                   << errno << " (" << std::strerror(errno) << ")"
                   << std::endl << std::flush;
       return Rcpp::wrap(false);
     }
     if (status != 0) {
       Rcpp::Rcerr << "[MakeBLASTDB] makeblastdb process returned non-zero status: "
                   << status << std::endl << std::flush;
       return Rcpp::wrap(false);
     }
   }
#endif
   
   // Environment pkg = Environment::namespace_env("QuickBLAST"); //https://teuder.github.io/rcpp4everyone_en/230_R_function.html#function
   // Function f = pkg[".GetLibsPath"];
   // std::string program_path = Rcpp::as<std::string>( f(Named("file_name")="makeblastdb") );
   //
   // std::string dbtype;
   // switch(ptr_->GetSeqType()){
   // case QuickBLAST::ESeqType::eNucleotide:{
   //   dbtype = "nucl";
   //   break;
   // }
   // case QuickBLAST::ESeqType::eProtein:{
   //   dbtype = "prot";
   //   break;
   // }
   // }
   //
   // std::vector<std::string> argv = {
   //   program_path, "-in", input_file_, "-dbtype", dbtype, "-out", database_name_, parse_seqids_opt
   // };
   //
   // // build argv char* array
   // std::vector<char*> cargv;
   // for (const auto &s : argv) cargv.push_back(const_cast<char*>(s.c_str()));
   // cargv.push_back(nullptr);
   //
   // pid_t pid;
   // int status = posix_spawnp(&pid, cargv[0], nullptr, nullptr, cargv.data(), environ);
   
   //  if (status != 0) {
   //    throw std::runtime_error("posix_spawnp failed: " + std::to_string(status));
   //  }
   //  if (waitpid(pid, &status, 0) == -1) {
   //    throw std::runtime_error("waitpid failed");
   //  }
   //  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
   //    throw std::runtime_error("makeblastdb failed (exit code " + std::to_string(WEXITSTATUS(status)) + ").");
   //  }
   
   if(!Rcpp::as<bool>(isBLASTDB(ptr, database_name)["is_db"])){
     return Rcpp::wrap(false);
   }
   
   Rcpp::Rcout << "makeblastdb finished successfully" << std::endl<< std::flush;
   
   return Rcpp::wrap(database_name_);
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("MakeBLASTDB() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("MakeBLASTDB(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("MakeBLASTDB() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "MakeBLASTDB() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false);
}

// SEXP MakeBLASTDB(SEXP ptr, SEXP input_file, SEXP database_name, bool parse_seqids = false, 
//                 SEXP stdout_opt = R_NilValue, SEXP stderr_opt = R_NilValue){
//  try {
//    if (TYPEOF(input_file) != STRSXP || input_file == R_NilValue) {
//      Rcpp::stop("input_file must be a single string (character vector of length 1)");
//    }
//    if (TYPEOF(database_name) != STRSXP || database_name == R_NilValue) {
//      Rcpp::stop("database_name must be a single string (character vector of length 1)");
//    }
//    
//    // Evaluate output routing: Default to TRUE if nothing is passed
//    Rcpp::RObject out_target = (stdout_opt == R_NilValue) ? Rcpp::RObject(Rcpp::wrap(true)) : Rcpp::RObject(stdout_opt);
//    Rcpp::RObject err_target = (stderr_opt == R_NilValue) ? Rcpp::RObject(Rcpp::wrap(true)) : Rcpp::RObject(stderr_opt);
//    
//    Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
//    std::string input_file_ = Rcpp::as<std::string>(input_file);
//    std::ifstream infile(input_file_);
//    if (!infile.good()) {
//      Rcpp::stop("MakeBLASTDB Error: " + input_file_ + " does not exist or cannot be read.");
//    }
//    std::string database_name_ = Rcpp::as<std::string>(database_name);
//    
//    Environment pkg = Environment::namespace_env("QuickBLAST"); 
//    Function f = pkg[".GetLibsPath"];
//    std::string mbdb_exe = "makeblastdb";
//    
// #if defined(_WIN32) || defined(__MINGW32__)
//    mbdb_exe += ".exe";
// #endif
//    std::string program_path = Rcpp::as<std::string>( f(Named("file_name")=mbdb_exe) );
//    
//    if (!FileExists(program_path)) {
//      Rcpp::stop("[MakeBLASTDB] Error: Executable not found at path: " + program_path);
//    }
//    
//    std::string dbtype;
//    switch(ptr_->GetSeqType()){
//    case QuickBLAST::ESeqType::eNucleotide: { dbtype = "nucl"; break; }
//    case QuickBLAST::ESeqType::eProtein:    { dbtype = "prot"; break; }
//    }
//    
//    // Build arguments specifically for system2
//    Rcpp::CharacterVector r_args;
//    r_args.push_back("-in");
//    r_args.push_back(input_file_);
//    r_args.push_back("-dbtype");
//    r_args.push_back(dbtype);
//    r_args.push_back("-out");
//    r_args.push_back(database_name_);
//    
//    if (parse_seqids) {
//      r_args.push_back("-parse_seqids");
//    }
//    
//    // Dynamically access R's system2
//    Rcpp::Environment base_env = Rcpp::Environment::base_env();
//    Rcpp::Function system2 = base_env["system2"];
//    
//    // Execute system2 safely via R's C++ API, passing the user-defined routing options
//    SEXP res = system2(
//      Rcpp::Named("command") = program_path,
//      Rcpp::Named("args")    = r_args,
//      Rcpp::Named("stdout")  = out_target, 
//      Rcpp::Named("stderr")  = err_target
//    );
//    
//    // Handle the variable return type from system2
//    int status = 0;
//    if (TYPEOF(res) == INTSXP || TYPEOF(res) == REALSXP) {
//      // Normal execution returning exit status
//      status = Rcpp::as<int>(res);
//    } else if (TYPEOF(res) == STRSXP) {
//      // If the user passed stdout=TRUE, system2 returns a CharacterVector 
//      // of the output and attaches the exit code as a "status" attribute.
//      Rcpp::CharacterVector cv(res);
//      if (cv.hasAttribute("status")) {
//        status = Rcpp::as<int>(cv.attr("status"));
//      }
//    }
//    
//    if (status != 0) {
//      Rcpp::stop("[MakeBLASTDB] process failed with non-zero exit status: " + std::to_string(status));
//    }
//    
//    if(!Rcpp::as<bool>(isBLASTDB(ptr, database_name)["is_db"])){
//      Rcpp::stop("[MakeBLASTDB] Database validation failed after building.");
//    }
//    
//    return Rcpp::wrap(database_name_);
//    
//  } catch(const Rcpp::exception &e) {
//    throw;
//  } catch(const std::exception &e) {
//    Rcpp::stop(std::string("MakeBLASTDB(): C++ Error : ") + e.what());
//  } catch(...) {
//    Rcpp::stop("MakeBLASTDB() - Unknown C++ Exception");
//  }
// }


//' @name BLAST2DBs
//'
//' @title BLAST 2 on-disk BLAST DBs
//'
//' @description Calls makeblastdb to create a BLAST DB of a FASTA file
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::MakeBLASTDB()], [QuickBLAST::isBLASTDB()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query DB
//' @param subject (string) Subject DB
//' @param out_file (string) Ouput file (Optional)
//' @param out_format (string) Ouput Format. 'ipc'/'csv'/'parquet' (Optional) (Default: 'parquet').
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param refresh_db (bool) If TRUE, re-creates the DBs
//' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size - Size of file write buffer (Optional).
//' @param enable_chunking (bool) Chunk large sequences? (Default: FALSE)
//' @param chunk_size (int) Size of chunks (Default: 50000)
//' @param overlap (int) Overlap between chunks (Default: 1000)
//' @param verbose (bool) Verbose? (Default: TRUE)
//' @return (SEXP) Rcpp::List - if return_values == TRUE, out_file - Otherwise.
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::BLAST2DBs(
//'   ptr=blastp_inst,
//'   query=system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   subject=system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   num_threads=24,
//'   out_file="test.db.arrow",
//'   return_values = T
//' )
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ), 
//'   "protein_query.db"
//' )
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   "protein_subject.db"
//' )
//' QuickBLAST::BLAST2DBs(
//'   ptr=blastp_inst,
//'   query="protein_query.db",
//'   subject="protein_subject.db",
//'   num_threads=24,
//'   out_file="test.db.arrow",
//'   return_values = T
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP BLAST2DBs(SEXP ptr, SEXP query, SEXP subject, SEXP out_file = R_NilValue, SEXP out_format = R_NilValue, unsigned int num_threads = 0, bool refresh_db = false, bool return_values = true, unsigned int min_batch_size = 0, const bool& enable_chunking = false, unsigned int chunk_size = 50000, unsigned int overlap = 1000, bool verbose = true){
 try{
   auto start = std::chrono::high_resolution_clock::now();
   unsigned int threads = DetectThreadLimit(num_threads);
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   if(out_file == R_NilValue && return_values == false){
     Rcpp::Rcerr << "Error: Both out_file cannot be empty and return_values == FALSE" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (query == R_NilValue || TYPEOF(query) != STRSXP || Rf_length(query) != 1) {
     Rcpp::Rcerr << "query must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject == R_NilValue || TYPEOF(subject) != STRSXP || Rf_length(subject) != 1) {
     Rcpp::Rcerr << "subject must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::string query_ = Rcpp::as<std::string>(query);
   std::string subject_ = Rcpp::as<std::string>(subject);
   
   // --- QUERY ---
   Rcpp::List q_db_info = isBLASTDB(ptr, query);
   if(refresh_db || !Rcpp::as<bool>(q_db_info["is_db"])) {
     // Not a BLAST DB, assume FASTA and try making BLAST DB
     query_ = query_ + ".db";
     MakeBLASTDB(ptr, query, Rcpp::wrap(query_), false);
   } else {
     // It IS a DB, but we MUST update query_ to the exactly resolved base path!
     // (Replace "db_name" with whatever key you used in isBLASTDB)
     query_ = Rcpp::as<std::string>(q_db_info["db_name"]); 
   }
   Rcpp::Rcout << "Q :" << query_ << std::endl << std::flush; //DEBUG
   
   // --- SUBJECT ---
   Rcpp::List s_db_info = isBLASTDB(ptr, subject);
   if(refresh_db || !Rcpp::as<bool>(s_db_info["is_db"])) {
     // Not a BLAST DB, assume FASTA and try making BLAST DB
     subject_ = subject_ + ".db";
     MakeBLASTDB(ptr, subject, Rcpp::wrap(subject_), false);
   } else {
     // Update subject_ to the exactly resolved base path
     subject_ = Rcpp::as<std::string>(s_db_info["db_name"]);
   }
   Rcpp::Rcout << "S :" << subject_ << std::endl << std::flush; //DEBUG
   
   // // Rcpp::Rcout << "Q isBLASTDB():" << Rcpp::as<bool>(isBLASTDB(ptr, query_)["is_db"]) << std::endl << std::flush; //DEBUG
   // if(!Rcpp::as<bool>(isBLASTDB(ptr, query_)["is_db"])){
   //   //Not a BLAST DB, assume FASTA and try making BLAST DB
   //   query_ = query_ + ".db";
   //   MakeBLASTDB(ptr, query, Rcpp::wrap(query_));
   // }
   // Rcpp::Rcout << "Q :" << query_ << std::endl << std::flush; //DEBUG
   // // Rcpp::Rcout << "S isBLASTDB():" << Rcpp::as<bool>(isBLASTDB(ptr, subject_)["is_db"]) << std::endl << std::flush; //DEBUG
   // if(!Rcpp::as<bool>(isBLASTDB(ptr, subject_)["is_db"])){
   //   //Not a BLAST DB, assume FASTA and try making BLAST DB
   //   subject_ = subject_ + ".db";
   //   MakeBLASTDB(ptr, subject, Rcpp::wrap(subject_));
   // }
   // Rcpp::Rcout << "S :" << subject_ << std::endl << std::flush; //DEBUG
   
   std::string out_file_ = out_file == R_NilValue ? "" : Rcpp::as<std::string>(out_file);
   std::string out_format_ = out_format == R_NilValue ? "parquet" : Rcpp::as<std::string>(out_format);   
   
   if (query_.empty()) {
     Rcpp::Rcerr << "query is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject_.empty()) {
     Rcpp::Rcerr << "subject is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (out_format_.empty()) {
     Rcpp::Rcerr << "out_format is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::shared_ptr<arrow::RecordBatchVector> ret_vals = std::make_shared<arrow::RecordBatchVector>();
   
   if (return_values)
   {
     ret_vals = ptr_->BLAST_dbs(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size, enable_chunking, chunk_size, overlap, verbose);
   }
   else
   {
     static_cast<void>(ptr_->BLAST_dbs(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size, enable_chunking, chunk_size, overlap, verbose));
   }
   
   if(verbose)
     PrintClock(start);
   
   if(return_values){
     if(ret_vals->size() > 0){
       auto ret_vals_ = RecordBatchVectorToFlattenedDFList_parallel(*ret_vals);
       ret_vals->clear();
       ret_vals->shrink_to_fit();
       return ret_vals_; //rm_null(ret_vals_);
     }else{
       return Rcpp::wrap(false);
     }
   }else{
     return Rcpp::wrap(true);
   }
   
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLAST2DBs() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLAST2DBs(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLAST2DBs() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLAST2DBs() - Unknown Exception" << std::endl <<std::flush;
 } 
 return Rcpp::wrap(false); 
}


// simple trim helpers
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char ch){ return !std::isspace(ch); }));
}
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
}
static inline void trim(std::string &s) { ltrim(s); rtrim(s); }

//' @name GetFASTAHeaders
//'
//' @title Get FASTA Headers
//'
//' @description Get the header strings of a FASTA file as a String Vector
//'
//' @seealso [QuickBLAST::MakeBLASTDB()], [QuickBLAST::isBLASTDB()]
//' @param path (std::string) Path to FASTA file
//' @param keep_gt (bool) Keep the '>' symbol? (Default: FALSE)
//' @return (SEXP) Rcpp::StringVector - on success, FALSE - Otherwise.
//' @examples
//' \dontrun{
//' QuickBLAST::GetFASTAHeaders(
//'   system.file(
//'     "extdata",
//'     "protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   keep_gt = F
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP GetFASTAHeaders(const std::string &path, bool keep_gt = false) {
 std::ifstream in(path, std::ios::in | std::ios::binary);
 if (!in.is_open()) {
   Rcpp::Rcerr << "Cannot open file: " + path << std::endl <<std::flush;
   return Rcpp::wrap(false); 
 }
 
 std::vector<std::string> headers;
 headers.reserve(1024); // sensible initial reservation
 
 std::string line;
 while (std::getline(in, line)) {
   // find first non-space character (handles files with stray leading whitespace)
   std::size_t i = 0;
   while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
   if (i < line.size() && line[i] == '>') {
     // get rest of line (including the '>' if requested)
     std::string hdr;
     if (keep_gt) {
       hdr = line.substr(i);        // include '>' and the rest
     } else {
       hdr = line.substr(i + 1);    // exclude the leading '>'
     }
     trim(hdr);
     headers.push_back(std::move(hdr));
   }
 }
 
 // move into Rcpp::StringVector
 Rcpp::StringVector out(headers.size());
 for (std::size_t j = 0; j < headers.size(); ++j) out[j] = headers[j];
 if(out.size() > 0)
   return out;
 else
   return Rcpp::wrap(false);
}

//' @name BLASTFile2DB
//'
//' @title BLAST File to a on-disk BLAST DB
//'
//' @description Runs BLAST using the ptr
//' @note Calls makeblastdb to create a BLAST DB of subject if it is not a DB
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::MakeBLASTDB()], [QuickBLAST::isBLASTDB()]
//' @param ptr (\code{Rcpp::XPtr<QuickBLAST>}) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query DB
//' @param subject (string) Subject DB
//' @param out_file (string) Ouput file (Optional)
//' @param out_format (string) Ouput Format. 'ipc'/'csv'/'parquet' (Optional) (Default: 'parquet').
// ' @param seq_limit (int) Batch Size to BLAST at a time. { -1 = Whole File, 0 - One sequence at a time or > 0 } - Size of BLAST sequences buffer (Optional)
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size - Size of file write buffer (Optional).
//' @return (SEXP) Rcpp::List - if return_values == TRUE, out_file - Otherwise.
//' @examples
//' \dontrun{
//' blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
//'   seq_type = 1,
//'   strand = 0,
//'   program = "blastp",
//'   save_sequences = F,
//'   save_hsp_sequences = F
//' )
//' QuickBLAST::BLASTFile2DB(
//'   ptr=blastp_inst,
//'   query=system.file(
//'     "extdata","protein_query.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   subject=system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   num_threads=24,
//'   out_file="test.db.arrow",
//'   return_values = T
//' )
//' QuickBLAST::MakeBLASTDB(
//'   blastp_inst,
//'   system.file(
//'     "extdata",
//'     "protein_subject.fasta",
//'     package = "QuickBLAST",
//'     mustWork = T
//'   ),
//'   "protein_subject.db"
//' )
//' QuickBLAST::BLASTFile2DB(
//'   ptr=blastp_inst,
//'   query="protein_query.fasta",
//'   subject="protein_subject.db",
//'   num_threads=24,
//'   out_file="test.db.arrow",
//'   return_values = T
//' )
//' }
//' @export
// [[Rcpp::export]]
SEXP BLASTFile2DB(SEXP ptr, SEXP query, SEXP subject, SEXP out_file = R_NilValue, SEXP out_format = R_NilValue, unsigned int num_threads = 0, bool return_values = true, unsigned int min_batch_size = 0){
 try{
   auto start = std::chrono::high_resolution_clock::now();
   unsigned int threads = DetectThreadLimit(num_threads);
   Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
   
   if(out_file == R_NilValue && return_values == false){
     Rcpp::Rcerr << "Error: Both out_file cannot be empty and return_values == FALSE" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (query == R_NilValue || TYPEOF(query) != STRSXP || Rf_length(query) != 1) {
     Rcpp::Rcerr << "query must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject == R_NilValue || TYPEOF(subject) != STRSXP || Rf_length(subject) != 1) {
     Rcpp::Rcerr << "subject must be a single string (character vector of length 1)" << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   std::string query_ = Rcpp::as<std::string>(query);
   std::string subject_ = Rcpp::as<std::string>(subject);
   
   if(Rcpp::as<bool>(isBLASTDB(ptr, query)["is_db"])){
     Rcpp::Rcerr<< "Error: Input query must be a FASTA file. Use QuickBLAST::BLAST2DBs() when inputs are BLAST DBs" << std::endl << std::flush; 
     return Rcpp::wrap(false);
   }
   if(!Rcpp::as<bool>(isBLASTDB(ptr, subject)["is_db"])){
     //Not a BLAST DB, assume FASTA and try making BLAST DB
     subject_ = subject_ + ".db";
     MakeBLASTDB(ptr, subject, Rcpp::wrap(subject_), false);
   }
   
   std::string out_file_ = out_file == R_NilValue ? "" : Rcpp::as<std::string>(out_file);
   std::string out_format_ = out_format == R_NilValue ? "parquet" : Rcpp::as<std::string>(out_format);   
   if (query_.empty()) {
     Rcpp::Rcerr << "query is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (subject_.empty()) {
     Rcpp::Rcerr << "subject is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   if (out_format_.empty()) {
     Rcpp::Rcerr << "out_format is empty." << std::endl <<std::flush;
     return Rcpp::wrap(false); 
   }
   
   if (!std::filesystem::exists(query_)) {
     Rcpp::Rcerr << "[BLASTFile2DB()] Error: Query file does not exist or is not accessible: " << query_ << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   if (!std::filesystem::exists(subject_)) {
     Rcpp::Rcerr << "[BLASTFile2DB()] Error: Subject file does not exist or is not accessible: " << subject_ << std::endl << std::flush;
     return Rcpp::wrap(false);
   }
   
   std::shared_ptr<arrow::RecordBatchVector> ret_vals = std::make_shared<arrow::RecordBatchVector>();
   
   if (return_values)
   {
     ret_vals = ptr_->BLAST_f2db(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size); 
   }
   else
   {
     static_cast<void>(ptr_->BLAST_f2db(query_, subject_, out_file_, out_format_, threads, return_values, min_batch_size));
   }
   
   PrintClock(start);
   
   if(return_values){
     if(ret_vals->size() > 0){
       auto ret_vals_ = RecordBatchVectorToFlattenedDFList_parallel(*ret_vals);
       ret_vals->clear();
       ret_vals->shrink_to_fit();
       return ret_vals_; //rm_null(ret_vals_);
     }else{
       return Rcpp::wrap(false);
     }
   }else{
     return Rcpp::wrap(true);
   }
   
 }catch(const Rcpp::exception &e){
   Rcpp::Rcerr << std::string("BLASTFile2DB() - Rcpp Exception : ") + e.what() << std::endl <<std::flush;
 }catch(const std::runtime_error &e){
   Rcpp::Rcerr << std::string("BLASTFile2DB(): C++ Runtime Error : ") + e.what() << std::endl <<std::flush;
 }catch(const std::exception &e){
   Rcpp::Rcerr << std::string("BLASTFile2DB() - C++ Exception : ") + e.what() << std::endl <<std::flush;
 }catch(...){
   Rcpp::Rcerr << "BLASTFile2DB() - Unknown Exception" << std::endl <<std::flush;
 }
 return Rcpp::wrap(false); 
}


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

// Extracts base primitives into FlatDF
static void ExtractPrimitive_worker(const std::shared_ptr<arrow::Array>& array,
                                    const std::shared_ptr<arrow::DataType>& type,
                                    const std::string& colname,
                                    int64_t nrows, FlatDF& out) {
  switch (type->id()) {
  case arrow::Type::STRING:
  case arrow::Type::LARGE_STRING: {
    add_string_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    if (type->id() == arrow::Type::STRING) {
      auto sarr = std::static_pointer_cast<arrow::StringArray>(array);
      for (int64_t i = 0; i < nrows; ++i) {
        if (sarr->IsValid(i)) {
          auto sv = sarr->GetView(i);
          cd.str_vals[i].assign(sv.data(), sv.size());
          cd.str_valid[i] = true;
        }
      }
    } else {
      auto sarr = std::static_pointer_cast<arrow::LargeStringArray>(array);
      for (int64_t i = 0; i < nrows; ++i) {
        if (sarr->IsValid(i)) {
          auto sv = sarr->GetView(i);
          cd.str_vals[i].assign(sv.data(), sv.size());
          cd.str_valid[i] = true;
        }
      }
    }
    break;
  }
  case arrow::Type::INT8: case arrow::Type::INT16: 
  case arrow::Type::INT32: case arrow::Type::INT64: {
    add_int_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    for (int64_t i = 0; i < nrows; ++i) {
      if (array->IsValid(i)) {
        if (type->id() == arrow::Type::INT32) cd.int_vals[i] = std::static_pointer_cast<arrow::Int32Array>(array)->Value(i);
        else if (type->id() == arrow::Type::INT64) cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int64Array>(array)->Value(i));
        else if (type->id() == arrow::Type::INT16) cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int16Array>(array)->Value(i));
        else cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int8Array>(array)->Value(i));
      }
    }
    break;
  }
  case arrow::Type::FLOAT: case arrow::Type::DOUBLE: {
    add_double_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    for (int64_t i = 0; i < nrows; ++i) {
      if (array->IsValid(i)) {
        if (type->id() == arrow::Type::DOUBLE) cd.dbl_vals[i] = std::static_pointer_cast<arrow::DoubleArray>(array)->Value(i);
        else cd.dbl_vals[i] = static_cast<double>(std::static_pointer_cast<arrow::FloatArray>(array)->Value(i));
      }
    }
    break;
  }
  case arrow::Type::BOOL: {
    add_logical_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    auto barr = std::static_pointer_cast<arrow::BooleanArray>(array);
    for (int64_t i = 0; i < nrows; ++i) {
      if (barr->IsValid(i)) cd.logical_vals[i] = barr->Value(i) ? 1 : 0;
    }
    break;
  }
  default: {
    add_string_column(out, colname, nrows); // Fallback to NA
    break;
  }
  }
}

// Extracts primitives from inside Lists into FlatDF
static void ExtractListPrimitive_worker(const std::shared_ptr<arrow::ListArray>& larr,
                                        const std::shared_ptr<arrow::Array>& values_arr,
                                        const std::shared_ptr<arrow::DataType>& type,
                                        const std::string& colname, int64_t p, int64_t nrows, FlatDF& out) {
  switch (type->id()) {
  case arrow::Type::STRING:
  case arrow::Type::LARGE_STRING: {
    add_string_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    auto sval = std::static_pointer_cast<arrow::StringArray>(values_arr); // Cast generalized for brevity
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (p < len && sval->IsValid(larr->value_offset(i) + p)) {
        auto sv = sval->GetView(larr->value_offset(i) + p);
        cd.str_vals[i].assign(sv.data(), sv.size());
        cd.str_valid[i] = true;
      }
    }
    break;
  }
  case arrow::Type::INT8: case arrow::Type::INT16: 
  case arrow::Type::INT32: case arrow::Type::INT64: {
    add_int_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (p < len) {
        int64_t idx = larr->value_offset(i) + p;
        if (values_arr->IsValid(idx)) {
          if (type->id() == arrow::Type::INT32) cd.int_vals[i] = std::static_pointer_cast<arrow::Int32Array>(values_arr)->Value(idx);
          else if (type->id() == arrow::Type::INT64) cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int64Array>(values_arr)->Value(idx));
          else if (type->id() == arrow::Type::INT16) cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int16Array>(values_arr)->Value(idx));
          else cd.int_vals[i] = static_cast<int>(std::static_pointer_cast<arrow::Int8Array>(values_arr)->Value(idx));
        }
      }
    }
    break;
  }
  case arrow::Type::FLOAT: case arrow::Type::DOUBLE: {
    add_double_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (p < len) {
        int64_t idx = larr->value_offset(i) + p;
        if (values_arr->IsValid(idx)) {
          if (type->id() == arrow::Type::DOUBLE) cd.dbl_vals[i] = std::static_pointer_cast<arrow::DoubleArray>(values_arr)->Value(idx);
          else cd.dbl_vals[i] = static_cast<double>(std::static_pointer_cast<arrow::FloatArray>(values_arr)->Value(idx));
        }
      }
    }
    break;
  }
  case arrow::Type::BOOL: {
    add_logical_column(out, colname, nrows);
    ColumnData& cd = out.cols.back();
    auto barr = std::static_pointer_cast<arrow::BooleanArray>(values_arr);
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (p < len) {
        int64_t idx = larr->value_offset(i) + p;
        if (barr->IsValid(idx)) cd.logical_vals[i] = barr->Value(idx) ? 1 : 0;
      }
    }
    break;
  }
  default: {
    add_string_column(out, colname, nrows); // Fallback to NA
    break;
  }
  }
}


static void CollectFlattenedColumns_worker(const std::shared_ptr<arrow::Array>& array,
                                           const std::shared_ptr<arrow::DataType>& type,
                                           const std::string& prefix,
                                           int64_t nrows,
                                           FlatDF& out) {
  if (!array) return;
  
  // 1. STRUCT LOGIC
  if (type->id() == arrow::Type::STRUCT) {
    auto sarr = std::static_pointer_cast<arrow::StructArray>(array);
    int nf = sarr->num_fields();
    for (int f = 0; f < nf; ++f) {
      std::string child_name = type->field(f)->name();
      // Using your original naming logic here
      // std::string colname = child_name; 
      std::string colname = prefix + "_" + child_name;
      CollectFlattenedColumns_worker(sarr->field(f), type->field(f)->type(), colname, nrows, out);
    }
    return;
  }
  
  // 2. LIST LOGIC
  if (type->id() == arrow::Type::LIST) {
    auto larr = std::static_pointer_cast<arrow::ListArray>(array);
    auto val_type = type->field(0)->type();
    auto values = larr->values();
    
    int64_t maxlen = 0;
    for (int64_t i = 0; i < nrows; ++i) {
      int64_t len = larr->value_length(i);
      if (len > maxlen) maxlen = len;
    }
    
    for (int64_t p = 0; p < maxlen; ++p) {
      if (val_type->id() == arrow::Type::STRUCT) {
        auto vals_struct = std::static_pointer_cast<arrow::StructArray>(values);
        int sub_nf = vals_struct->num_fields();
        for (int sf = 0; sf < sub_nf; ++sf) {
          std::string sf_name = val_type->field(sf)->name();
          std::string colname = prefix + "_" + std::to_string(p) + "_" + sf_name;
          auto subfield_array = vals_struct->field(sf);
          auto sub_type = val_type->field(sf)->type();
          ExtractListPrimitive_worker(larr, subfield_array, sub_type, colname, p, nrows, out);
        }
      } else {
        std::string colname = prefix + "_" + std::to_string(p);
        ExtractListPrimitive_worker(larr, values, val_type, colname, p, nrows, out);
      }
    }
    return;
  }
  
  // 3. BASE PRIMITIVE LOGIC
  ExtractPrimitive_worker(array, type, prefix, nrows, out);
}

// Worker: convert a single RecordBatch into FlatDF (C++ only).
FlatDF process_rb_to_flatdf_worker(const std::shared_ptr<arrow::RecordBatch>& rb) {
  FlatDF out;
  if (!rb) return out;
  
  int64_t nrows = rb->num_rows();
  if (nrows == 0) return out;
  
  out.nrows = nrows;
  
  int nfields = rb->schema()->num_fields();
  for (int i = 0; i < nfields; ++i) {
    auto field = rb->schema()->field(i);
    auto arr = rb->column(i);
    
    // Hand it off to the recursive flattener
    CollectFlattenedColumns_worker(arr, field->type(), field->name(), nrows, out);
  }
  
  return out;
}

// Convert one FlatDF into an Rcpp data.frame (must run on main thread)
static Rcpp::List convert_flatdf_to_R(const FlatDF &flat) {
  int ncols = static_cast<int>(flat.cols.size());
  int nrows = static_cast<int>(flat.nrows);
  Rcpp::List out(ncols);
  Rcpp::CharacterVector names(ncols);
  
  for (int j = 0; j < ncols; ++j) {
    const ColumnData &cd = flat.cols[j];
    names[j] = cd.name;
    
    switch(cd.kind) {
    case ColKind::INT: {
      Rcpp::IntegerVector v(nrows);
      // memcpy is safe here because NA_INTEGER matches bit-for-bit
      if(nrows > 0){ 
        std::copy(cd.int_vals.begin(), cd.int_vals.end(), v.begin());
      }
      out[j] = v;
      break;
    }
    case ColKind::DOUBLE: {
      Rcpp::NumericVector v(nrows);
      if(nrows > 0){ 
        std::copy(cd.dbl_vals.begin(), cd.dbl_vals.end(), v.begin());
      }
      out[j] = v;
      break;
    }
    case ColKind::STRING: {
      Rcpp::StringVector v(nrows);
      if(nrows > 0){ 
        for (int i = 0; i < nrows; ++i) {
          // v[i] = cd.str_vals[i].empty() ? NA_STRING : Rcpp::String(cd.str_vals[i]);
          v[i] = cd.str_valid[i] ? Rcpp::String(cd.str_vals[i]) : NA_STRING;
        }
      }
      out[j] = v;
      break;
    }
    case ColKind::LOGICAL: {
      Rcpp::LogicalVector v(nrows);
      if(nrows > 0){ 
        for (int i = 0; i < nrows; ++i) {
          int marker = cd.logical_vals[i];
          v[i] = (marker == INT_MIN) ? NA_LOGICAL : (marker ? TRUE : FALSE);
        }
      }
      out[j] = v;
      break;
    }
    }
  }
  
  out.attr("names") = names;
  out.attr("class") = Rcpp::CharacterVector::create("data.frame");
  // Fast, compact row names creation
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -nrows);
  return out;
}

// Public: parallel RecordBatchVector -> list of data.frames (R)
Rcpp::List RecordBatchVectorToFlattenedDFList_parallel(const arrow::RecordBatchVector &rbv) {
  size_t m = rbv.size();
  if (m == 0) {
    // return Rcpp::List::create(Rcpp::Named("rbv_hits") = 0);
    return R_NilValue; 
  }
  
  Rcpp::Rcout << "RecordBatchVector size: " << m << std::endl;
  
  int64_t total_rows = 0;
  for (size_t i = 0; i < m; ++i) {
    if (rbv[i]) {
      total_rows += rbv[i]->num_rows();
    }
  }
  
  if (total_rows == 0) {
    return R_NilValue;
  }
  
  Rcpp::Rcout << "Total rows across all batches: " << total_rows << std::endl;
  
  // Create a standard C++ vector to hold the background work
  std::vector<FlatDF> intermediate_results(m);
  
  // 1. OPENMP PARALLEL REGION (No R API calls allowed inside)
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < m; ++i) {
    try {
      intermediate_results[i] = process_rb_to_flatdf_worker(rbv[i]);
    } catch (...) {
      intermediate_results[i] = FlatDF(); // Empty fallback on error
    }
  }
  
  // 2. MAIN THREAD REGION (R API is safe here)
  Rcpp::List out(m);
  for (size_t i = 0; i < m; ++i) {
    if (intermediate_results[i].cols.empty()) {
      out[i] = Rcpp::List::create();
    } else {
      out[i] = convert_flatdf_to_R(intermediate_results[i]);
    }
  }
  
  arrow::default_memory_pool()->ReleaseUnused();
  return out;
}
