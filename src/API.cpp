#include <algo/blast/QuickBLAST/commons.hpp>
// #include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/API.hpp>

// namespace Log
// {
//     // #ifdef RCPP_USE_GLOBAL_ROSTREAM
//     Rcpp::Rostream<true> &Rcppcout = Rcpp::Rcpp_cout_get();
//     Rcpp::Rostream<false> &Rcppcerr = Rcpp::Rcpp_cerr_get();
//     // #endif
// } // namespace Log

using namespace Rcpp;

// // #if defined(MINGW32) || defined(WIN32)
// extern "C"
// {
//     //     QBLIBRARY_API BOOL QuickBLASTcpp_main(std::string dllPath){
//     //         //std::string dllPath = "inst/x64/QuickBLASTcpp.dll";
// 
//     //         // Load the DLL using LoadLibraryEx
//     //         HMODULE hDLL = LoadLibraryExA(dllPath.data(), NULL, 0);
// 
//     //         if (hDLL == NULL) {
//     //             std::cout << "Failed to load the DLL: " << dllPath<< std::endl;  // Handle DLL loading failure
//     //             return FALSE;
//     //         }
// 
//     //         //FreeLibrary(hDLL);
//     //         return TRUE;
//     //     }
// 
//     // QuickBLASTHandle GetQuickBLASTInstance(unsigned int id);
//     QBLIBRARY_API unsigned int libQB_GetInstanceCount();
//     unsigned int GetInstanceID(QuickBLASTHandle ptr);
//     // QBLIBRARY_API unsigned int libQB_CreateQuickBLASTInstance(const int seq_type, const int strand, const char *program, const char *options, int save_sequences);
//     QBLIBRARY_API SEXP libQB_CreateQuickBLASTInstance(const int seq_type, const int strand, const char *program, const char *options, int save_sequences, const unsigned int num_threads);
// 
//     QBLIBRARY_API SEXP libQB_DeleteQuickBLASTInstance(SEXP ptr_id);
//     //     //QBLIBRARY_API ArrowWrapper* CreateArrowWrapperInstance();
//     //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
//     //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
//     QBLIBRARY_API bool libQB_SetQuickBLASTOptions(SEXP ptr_id, SEXP program_name, SEXP options);
//     QBLIBRARY_API SEXP BLAST2Seqs(SEXP ptr_id, SEXP query, SEXP subject);
//     QBLIBRARY_API SEXP BLAST2Files(SEXP ptr_id, SEXP query, SEXP subject, SEXP out_file, SEXP seq_limit, SEXP num_threads, SEXP show_progress, SEXP return_values, SEXP min_batch_size);
//     // QBLIBRARY_API SEXP BLAST2Folders(int ptr_id, std::string query, std::string subject, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
//     // QBLIBRARY_API SEXP BLAST1Folder(int ptr_id, std::string input_folder, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
//     // // QBLIBRARY_API std::string getFilenameWithoutExtension(const std::string &filename);
//     // // QBLIBRARY_API Rcpp::List rm_null(Rcpp::List x);
//     // // QBLIBRARY_API std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension);
//     // QBLIBRARY_API SEXP libQB_isQuickBLASTLoaded();
//     QBLIBRARY_API std::string libQB_isQuickBLASTLoaded();
//     //     //QBLIBRARY_API  bool QueryOOBESupport() { return false; }
//     //     /*QBLIBRARY_API int arrow_struct_num_fields(std::shared_ptr<arrow::StructArray> arr);
//     //     QBLIBRARY_API std::shared_ptr<arrow::Array> arrow_struct_field(std::shared_ptr<arrow::StructArray> arr, const int i);
//     //     QBLIBRARY_API std::shared_ptr<arrow::Array> arrow_array_slice(std::shared_ptr<arrow::Array> arr, const int offset, const int length);
//     //     QBLIBRARY_API int arrow_schema_num_fields(std::shared_ptr<arrow::Schema> sch);
//     //     QBLIBRARY_API std::shared_ptr<arrow::Field> arrow_schema_field(std::shared_ptr<arrow::Schema> sch, const int i);
//     //     QBLIBRARY_API std::shared_ptr<arrow::DataType> arrow_schema_field_type(std::shared_ptr<arrow::Schema> sch, const int i);
//     //     QBLIBRARY_API std::string arrow_schema_field_name(std::shared_ptr<arrow::Schema> sch, const int i);
//     //     QBLIBRARY_API bool arrow_int8array_isvalid(std::shared_ptr<arrow::Int8Array> arr, int i);
//     //     QBLIBRARY_API bool arrow_strarray_isvalid(std::shared_ptr<arrow::StringArray> arr, int i);
//     //     QBLIBRARY_API bool arrow_dblarray_isvalid(std::shared_ptr<arrow::DoubleArray> arr, int i);*/
// } // extern "C"
// 
// // #endif

//' @name isQuickBLASTLoaded
//' @title Check R <-> C++ (FFI) connection
//'
//' @description This function does nothing than check the connection between the R package and C++ libraries
//'
//' @return String that successfully confirms when the package is loaded properly
//' @export
// [[Rcpp::export]]
RcppExport bool isQuickBLASTLoaded() //SEXP libQB_isQuickBLASTLoaded()
{
    // std::string ret_str = "C++ - QuickBLAST dependencies Loaded!";
    // // Rprintf("%s - R print\n", ret_str.c_str());
    //  Rcpp::Rcout << ret_str.c_str() << " - R print"<< std::endl;
    //  Rcpp::Rcout << std::flush;
    ArrowWrapper *testwrap = new ArrowWrapper();
    std::shared_ptr<ArrowWrapper> testwrap_ = std::make_shared<ArrowWrapper>();
    // // return Rcpp::wrap(ret_str);
    // return ret_str;
    return true;
}

void PrintClock(std::chrono::time_point<std::chrono::high_resolution_clock> start)
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    // Print the time in seconds
    // Rprintf("Clock : %f seconds \n", elapsed_seconds.count());
     Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;
     Rcpp::Rcout << std::flush;
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

            auto sublist_array = list_array->values()->Slice(list_array->value_offset(i), list_array->value_length(i)); // arrow_array_slice(list_array->values(), list_array->value_offset(i), list_array->value_length(i));

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
            // if (arrow_strarray_isvalid(string_array, i))
            if (string_array->IsValid(i))
            {
                strings[i] = Rcpp::String(string_array->GetString(i));
            }
            // else
            // {
            //   strings[i] = NA_STRING;
            // }
        }

        return strings;
    }
    else if (type->id() == arrow::Type::INT8)
    {
        auto int_array = std::static_pointer_cast<arrow::Int8Array>(array);

        Rcpp::IntegerVector ints(int_array->length());

        for (int i = 0; i < int_array->length(); ++i)
        {
            // if (arrow_int8array_isvalid(int_array, i))
            if (int_array->IsValid(i))
            {
                ints[i] = int_array->Value(i);
            }
            // else
            // {
            //   ints[i] = NA_INTEGER;
            // }
        }

        return ints;
    }
    else if (type->id() == arrow::Type::DOUBLE)
    { // Use arrow::Type::DOUBLE instead of FLOAT64

        auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
        Rcpp::NumericVector doubles(double_array->length());

        for (int i = 0; i < double_array->length(); ++i)
        {
            // if (arrow_dblarray_isvalid(double_array, i))
            if (double_array->IsValid(i))
            {
                doubles[field_name] = double_array->Value(i);
            }
            // else
            // {
            //   doubles[i] = NA_REAL;
            // }
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
        Rcpp::checkUserInterrupt();
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
        Rcpp::checkUserInterrupt();
        std::shared_ptr<arrow::RecordBatch> rb = rb_vector[i];
        result_list[i] = Hits2RList(rb);
    }

    return result_list;
}

std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension)
{
    std::vector<std::string> outFiles;

    if (!extension.empty())
    {
        for (const auto &entry : std::filesystem::directory_iterator(folderPath))
        {
            if (entry.is_regular_file() && entry.path().extension() == extension)
            {
                outFiles.emplace_back(entry.path().filename().string());
            }
        }
    }
    else
    {
        for (const auto &entry : std::filesystem::directory_iterator(folderPath))
        {
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

// QuickBLASTHandle GetQuickBLASTInstance(unsigned int id)
// {
//     assert(id >= 0);
//     QuickBLASTHandle handle;
//     // handle.ptr = cppObj_list[id].get();
//     handle.ptr = cppObj_list[id];
//     assert(handle.ptr != nullptr);
//     handle.id = id;
//     return handle;
// }

//' @name GetInstanceCount
//' @title Get count of QuickBLAST instances stored in C++ side
//'
//' @description This function gives the size of the list of QuickBLAST C++ object list
//'
//' @return Count of QuickBLAST instances
//' @export
// [[Rcpp::export]]
RcppExport unsigned int GetInstanceCount()
{
   //  // Rprintf("testing c++ side");
   //  // Rprintf("\n%d\n", (int)cppObj_list.size());
   //   Rcpp::Rcout << "testing c++ side" << std::endl;
   // Rcpp::Rcout << (int)cppObj_list.size() << std::endl;
   // Rcpp::Rcout << std::flush;
    // Rprintf("\n%s\n", (int)cppObj_list.size());
    // if (!cppObj_list.empty())
    // {
    return (unsigned int)cppObj_list.size();
    // }
    // else
    // {
    //     return Rcpp::wrap(0);
    // }
    // return Rcpp::wrap((int)cppObj_list.size());
}

//' @name GetInstanceID
//' @title Get ID/Index of a QuickBLAST instance stored in C++ side
//'
//' @description This function fetches the ID/Index of a QuickBLAST instance of a Rcpp::XPtr<QuickBLAST> stored in C++ side.
//'
//' @param ptr (Rcpp::XPtr<QuickBLAST>) Pointer of QuickBLAST instance from the R side.
//'
//' @return (unsigned int) ID/Index of the QuickBLAST instance pointer, FALSE otherwise
//' @export
// [[Rcpp::export]]
RcppExport SEXP GetInstanceID(SEXP ptr)
{
  try
  {
    // return ptr.ptr->obj_id;
    Rcpp::XPtr<QuickBLAST> ptr_ = Rcpp::as<Rcpp::XPtr<QuickBLAST>>(ptr);
    // Find the first key associated with the target value
    auto it = std::find_if(cppObj_list.begin(), cppObj_list.end(),
                           [&](const std::pair<unsigned int, Rcpp::XPtr<QuickBLAST>>& pair) {
                             return pair.second.get() == ptr_.get();
                           });
    
    if (it != cppObj_list.end()) {
      return Rcpp::wrap(it->first);
    }
  }catch (const std::exception &e)
  {
    Rprintf("Error fetching ID/Index of QuickBLAST instance: %s\n", e.what());
  }
  return Rcpp::wrap(false);
  // return ptr.ptr->GetObjectID();
}

//' @name GetQuickBLASTInstance
//' @title Get QuickBLAST instance stored in C++ side at ID/Index
//'
//' @description This function fetches the QuickBLAST instance of a Rcpp::XPtr<QuickBLAST> at ID/Index stored in C++ side.
//'
//' @param ptr_id (unsigned int) ID/Index of Pointer to a QuickBLAST instance (in C++ side).
//'
//' @return (Rcpp::XPtr<QuickBLAST>) Pointer to a QuickBLAST instance, FALSE otherwise
//' @export
// [[Rcpp::export]]
RcppExport SEXP GetQuickBLASTInstance(unsigned int ptr_id)
{
  try {
    return cppObj_list.at(ptr_id);
  } catch (const std::out_of_range& e) {
    std::cerr << "Error fetching QuickBLAST instance (" << ptr_id << "): " << e.what() << std::endl;
  }
  return Rcpp::wrap(false);
}

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
        REprintf("Only named list or string allowed for BLAST options : Check QuickBLAST::GetAvailableBLASTOptions()\n");
        // return Rcpp::wrap(false);
        break;
    }
    return options_;
}



unsigned int DetectThreadLimit(unsigned int num_threads){
  unsigned int hw = std::thread::hardware_concurrency();
  if (hw == 0) hw = 1;                      // fallback if implementation returns 0
  unsigned int threads = num_threads ? num_threads > hw ? hw : num_threads : hw;
  
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
//' @param save_sequences (bool) Save Sequences in the arrow::RecordBatch?.
//' @param num_threads (int) Number of threads.
//' @return (Rcpp::XPtr<QuickBLAST>) Pointer to a QuickBLAST Instance (Cannot be used in R)
//' @export
// [[Rcpp::export]]
RcppExport SEXP CreateQuickBLASTInstance(const int seq_type, const int strand, SEXP program, SEXP options = R_NilValue, const bool save_sequences = false, const unsigned int num_threads = 0)
{
    try
    {
      // Auto-detect threads when 0 is passed
        unsigned int threads = DetectThreadLimit(num_threads);
      
        // return std::make_shared<QuickBLAST>(seq_type, strand, program, options, save_sequences);
        //  new QuickBLAST(seq_type, strand, program, options, save_sequences);

        QuickBLAST::ESeqType seq_type_ = static_cast<QuickBLAST::ESeqType>(seq_type);
        QuickBLAST::EStrand strand_ = static_cast<QuickBLAST::EStrand>(strand);
        std::string program_ = Rcpp::as<std::string>(program);

        std::string options_ = options == R_NilValue ? "" : ConvertBLASTOptions2String(options);
        // std::string options_ = options;
        bool save_sequences_ = save_sequences; //static_cast<bool>(save_sequences);
        // Rprintf("DBG: %d %d %s %s %d \n", seq_type_, strand_, program_.c_str(), options_.c_str(), save_sequences_);
        // QuickBLASTHandle handle;
        // std::shared_ptr<QuickBLAST> objPtr = std::make_shared<QuickBLAST>(seq_type_, strand_, program_, options_, save_sequences_);
        QuickBLAST *objPtr = new QuickBLAST(seq_type_, strand_, program_, options_, save_sequences_);
        objPtr->SetThreadCount(threads);
        Rcpp::XPtr<QuickBLAST> ptr(objPtr, true);
        // ptr.attr("class") = "QuickBLAST_XPtr";
        Rf_classgets(ptr, Rf_mkString("QuickBLAST_XPtr"));
        // Rprintf("DBG1: %d \n", cppObj_list.size());
        unsigned int list_size = cppObj_list.size();
        // handle.id = list_size;
        // // handle.ptr = objPtr.get();
        // handle.ptr = objPtr;
        // Log::Rcppcout << "dbg2" << std::endl;
        // cppObj_list.insert(std::make_pair(list_size, objPtr));
        cppObj_list.insert(std::make_pair(list_size, ptr));
        // Rprintf("DBG2: %d \n", cppObj_list.size());
        // // Log::Rcppcout << "dbg3" << std::endl;
        // return list_size;
        return ptr;
    }
    catch (const std::exception &e)
    {
        Rprintf("C++ side Exception: %s\n", e.what());
    }
    return Rcpp::wrap(false);
    // return 0;
}

//' @name DeleteQuickBLASTInstance
//' @title Delete a QuickBLAST instance stored in C++ side
//'
//' @description This function deletes a QuickBLAST instance based on the instance ID
//'
//' @param ptr_id (unsigned int) ID/Index of stored QuickBLAST instance.
//'
//' @return TRUE - if the instance is deleted successfully, FALSE otherwise
//' @export
// [[Rcpp::export]]
RcppExport bool DeleteQuickBLASTInstance(const unsigned int ptr_id)
{
  try
    {
    // unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    //// ptr.reset();
    // fetch the object from the map
    Rprintf("Deleting QuickBLAST Instance : %d\n", ptr_id);
    Rcpp::as<Rcpp::XPtr<QuickBLAST>>(GetQuickBLASTInstance(ptr_id))->~QuickBLAST();
    cppObj_list.erase(ptr_id);
    // // ptr.ptr->~QuickBLAST();
    // //// delete ptr.ptr.get();
    // Rcpp::XPtr<QuickBLAST> ptr = Rcpp::as<Rcpp::XPtr<QuickBLAST>>(ptr_id);
    // ptr->~QuickBLAST();
    // // return true;
    // return Rcpp::wrap(true);
    return true;
    }
  catch (const std::exception &e)
    {
      Rprintf("Error deleting QuickBLAST instance (%d): %s\n", ptr_id, e.what());
    }
  return false;
}

//' @name SetQuickBLASTOptions
//'
//' @title Set BLAST options for a QuickBLAST instance.
//'
//' @description Set/Modify the BLAST options for a QuickBLAST instance.
//'
//' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
//' @param ptr (Rcpp::XPtr<QuickBLAST>) Pointer to a QuickBLAST Instance
//' @param program_name (string) Name of the BLAST program
//' @param options (string (or) Named List) List of BLAST options - check QuickBLAST::GetAvailableBLASTOptions(). String should be of the format "-option1 value1 -option2 value2"
//' @return (bool) TRUE - if options set for the QuickBLAST instance, FALSE otherwise
//' @export
// [[Rcpp::export]]
RcppExport bool SetQuickBLASTOptions(SEXP ptr, SEXP program_name, SEXP options)
{
  try
  {
    // unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    std::string program_name_ = Rcpp::as<std::string>(program_name);
    std::string options_ = ConvertBLASTOptions2String(options);

    // auto ptr = GetQuickBLASTInstance(ptr_id_);
    Rcpp::XPtr<QuickBLAST> ptr_ = Rcpp::as<Rcpp::XPtr<QuickBLAST>>(ptr);
    // ptr.ptr->GetQuickBLASTOptions() = *ptr.ptr->SetQuickBLASTOptions(program_name_, options_);
    ptr_->GetQuickBLASTOptions() = *ptr_->SetQuickBLASTOptions(program_name_, options_);
    return true;
  }catch (const std::exception &e)
  {
    Rprintf("Error setting options for QuickBLAST instance: %s\n", e.what());
  }
  return false;
}

Rcpp::XPtr<QuickBLAST> ResolveQuickBLASTInstance(SEXP inst)
{
  switch (TYPEOF(inst)) {
  case INTSXP: 
  case REALSXP: {
    // numeric 
    // return GetQuickBLASTInstance(Rcpp::as<unsigned int>(inst));
    unsigned int ptr_id = Rcpp::as<unsigned int>(inst);
    try {
      return cppObj_list.at(ptr_id);
    } catch (const std::out_of_range& e) {
      Rcpp::stop("ResolveQuickBLASTInstance() - Error fetching QuickBLAST instance (%d): %s", ptr_id, e.what());
    }
    break;
  }
  case EXTPTRSXP: {
    // external pointer (XPtr)
    if (!Rf_inherits(inst, "QuickBLAST_XPtr")) {
      // helpful error message
      Rcpp::stop("Expected external pointer of class 'Rcpp::XPtr<QuickBLAST>' :: %s :: %s", Rf_inherits(inst, "QuickBLAST_XPtr"), Rf_type2char(TYPEOF(inst)));
    }
    // construct XPtr<QuickBLAST> from SEXP
    Rcpp::XPtr<QuickBLAST> xp(inst);
    if (!xp) Rcpp::stop("Received NULL QuickBLAST pointer");
    return xp;
    break;
  }
  default:
    Rcpp::stop("Unsupported SEXP type: %s", Rf_type2char(TYPEOF(inst)));
  }
}

//' @name BLAST2Seqs
//'
//' @title BLAST 2 Sequence strings
//'
//' @description BLAST 2 nucleotide or protein strings with a QuickBLAST instance.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (Rcpp::XPtr<QuickBLAST>) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query sequence
//' @param subject (string) Subject sequence.
//' @return (Rcpp::XPtr<QuickBLAST>) Pointer to a QuickBLAST Instance (Cannot be used in R)
//' @export
// [[Rcpp::export]]
RcppExport SEXP BLAST2Seqs(SEXP ptr, SEXP query, SEXP subject)
{
  
    assert(TYPEOF(query) == CHARSXP);
    assert(TYPEOF(subject) == CHARSXP);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
    // // unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    // // QuickBLASTHandle ptr = GetQuickBLASTInstance(ptr_id_);
    // Rcpp::XPtr<QuickBLAST> ptr = Rcpp::as<Rcpp::XPtr<QuickBLAST>>(ptr_id);
    std::string query_ = Rcpp::as<std::string>(query);
    std::string subject_ = Rcpp::as<std::string>(subject);

    assert(!query_.empty());
    assert(!subject_.empty());
    // ArrowRBVHandle ret_val;
    std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
    // std::shared_ptr<arrow::RecordBatch> ret_rb = ptr.ptr->BLAST_seqs(query_, subject_); //
    std::shared_ptr<arrow::RecordBatch> ret_rb = ptr_->BLAST_seqs(query_, subject_); //
    // ret_val.ptr.emplace_back(ptr.ptr->BLAST_seqs(query_, subject_)); //
    if (ret_rb)
    {
      try{
        arrow::Status rb_sts = ret_rb->ValidateFull();
        // Rcpp::Rcout << "here18.1.1:" << rb_sts.message()  << std::endl << rb_sts.ToString() << std::endl << "rows:" << ret_rb->num_rows() << "\ncols:" << ret_rb->num_columns()  << std::endl << std::flush; //DEBUG
        if (!rb_sts.ok())
        {
            REprintf("ERR : Invalid RB : %s \n %s \n", rb_sts.message().c_str(), rb_sts.detail()->ToString().c_str());
            // PrintClock(start);
            // ret_val.ptr.clear();
            // ret_val.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
            ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
            // return arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie();
        }
        else
        {
          // Rcpp::Rcout << "here18.1.2:" << rb_sts.ok()  << std::endl << std::flush; //DEBUG
          if(ret_rb->num_rows() <= 0){
            Rcpp::stop("No alignments could be computed");
          }
            ret_val->emplace_back(ret_rb); //
            // Rcpp::Rcout << "here18.1.3"  << std::endl << std::flush; //DEBUG
        }
      }catch(...){
        ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
        throw;
      }
    }
    else
    {
        // PrintClock(start);
        // ret_val.ptr.clear();
        // ret_val.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
        ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
        // return arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie();
    }

    Rcpp::List ret_vals_ = as<List>(Hits2RList(*ret_val));

    PrintClock(start);
    return rm_null(ret_vals_);
}



//' @name BLAST2Folders
//'
//' @title BLAST 2 Folders
//'
//' @description BLAST 2 Folders with FASTA files containing nucleotide or protein sequences with a QuickBLAST instance. The files from query and subject folders are selected with the extension parameter
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (Rcpp::XPtr<QuickBLAST>) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query folder
//' @param subject (string) Subject folder.
//' @param extension (string) Extension of files in folder.
//' @param out_folder (string) Ouput Folder (Must be empty).
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param reciprocal_hits (bool) Perform Bi-directional (Reciprocal => query <-> subject) BLAST? (Default: FALSE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size (Optional).
//' @return (bool) TRUE - on success, FALSE - Otherwise. (Results are not returned as R Lists to reduce overhead)
//' @export
// [[Rcpp::export]]
RcppExport bool BLAST2Folders(SEXP ptr, SEXP query, SEXP subject, SEXP extension, SEXP out_folder, unsigned int num_threads = 0, bool reciprocal_hits = false, unsigned int min_batch_size = 0)
{
    // int typ1 = TYPEOF(query);
    // int typ2 = TYPEOF(subject);
    auto start = std::chrono::high_resolution_clock::now();
    unsigned int threads = DetectThreadLimit(num_threads);
    Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
    assert(TYPEOF(out_folder) == CHARSXP);
    assert(TYPEOF(query) == CHARSXP);
    assert(TYPEOF(subject) == CHARSXP);
    assert(TYPEOF(extension) == CHARSXP);
    // assert(typ1 == LISTSXP || typ1 == VECSXP || typ1 == CHARSXP);
    // assert(typ2 == LISTSXP || typ2 == VECSXP || typ2 == CHARSXP);
    // assert(ptr.ptr != nullptr);
    
    std::string query_ = Rcpp::as<std::string>(query);
    std::string subject_ = Rcpp::as<std::string>(subject);
    std::string out_folder_ = Rcpp::as<std::string>(out_folder);
    std::string extension_ = Rcpp::as<std::string>(extension);
    
    assert(!query_.empty());
    assert(!subject_.empty());
    assert(!out_folder_.empty());
    
    int seq_limit = -1;

    std::filesystem::path outPath(out_folder_);
    std::filesystem::create_directory(outPath);
    if (!std::filesystem::is_empty(outPath))
    {
      Rcpp::stop("out_folder : %s - Folder must be empty.", out_folder_);
      // std::cerr << "out_folder : Folder must be empty.";
      return false;
    }
    // assert(min_batch_size > 0);

    if(min_batch_size == 0){
      min_batch_size = threads;
    }
    
    std::filesystem::path query_path(query_);
    std::filesystem::path subject_path(subject_);
    std::vector<std::string> queryFiles = getFilesInDir(query_, extension_);
    std::vector<std::string> subjectFiles = getFilesInDir(subject_, extension_);

    assert(!queryFiles.empty());
    assert(!subjectFiles.empty());
    assert(queryFiles.size() > 0);
    assert(subjectFiles.size() > 0);
    assert(queryFiles.size() * subjectFiles.size() > 0);
    
    int iterations = queryFiles.size() * subjectFiles.size();

    // assert(iterations > 0);

    // Rcpp::List ret_lst(iterations);
    // Rcpp::CharacterVector names(iterations);
    Progress progress_bar(iterations, true);

    for (int i = 0; i < (int)queryFiles.size(); i++)
    {
        for (int j = 0; j < (int)subjectFiles.size(); j++)
        {
            Rcpp::checkUserInterrupt();
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
                std::cerr << "Warn : File not found : " << query_input.string() << " or " << subject_input.string() << std::endl;
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
            std::cout << "BLASTing : " << base_name << std::endl;
            static_cast<void>(ptr_->BLAST_files(query_input.string(), subject_input.string(), outFileStr, seq_limit, threads, true, false, min_batch_size));

            // static_cast<void>(cpp_BLAST2Files(sh_ptr, query_input.string(), subject_input.string(), outFile_.string(), seq_limit, num_threads_, true, false, min_batch_size_));
            // ret_lst[one_dim_index] = outFile_.string();
            // names[one_dim_index] = base_name;
        }
    }
    // ret_lst.names() = names;
    PrintClock(start);
    // return ret_lst;
    return true;
}

//' @name BLAST1Folder
//'
//' @title BLAST Files within a Folder
//'
//' @description BLAST FASTA files containing nucleotide or protein sequences, within a folder with a QuickBLAST instance. The files from the folder are selected with the extension parameter and BLAST'd against each other.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (Rcpp::XPtr<QuickBLAST>) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param input_folder (string) Input folder
//' @param extension (string) Extension of files in folder.
//' @param out_folder (string) Ouput Folder (Must be empty).
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param reciprocal_hits (bool) Perform Bi-directional (Reciprocal => query <-> subject) BLAST? (Default: FALSE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size (Optional).
//' @return (bool) TRUE - on success, FALSE - Otherwise. (Results are not returned as R Lists to reduce overhead)
//' @export
// [[Rcpp::export]]
RcppExport bool BLAST1Folder(SEXP ptr, SEXP input_folder, SEXP extension, SEXP out_folder, unsigned int num_threads = 0, bool reciprocal_hits = false, unsigned int min_batch_size = 0)
{
    // int typ1 = TYPEOF(input_folder);
    // assert(typ1 == LISTSXP || typ1 == VECSXP);

    auto start = std::chrono::high_resolution_clock::now();
    unsigned int threads = DetectThreadLimit(num_threads);
    Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
    assert(TYPEOF(out_folder) == CHARSXP);
    assert(TYPEOF(input_folder) == CHARSXP);
    assert(TYPEOF(extension) == CHARSXP);
    // assert(typ1 == LISTSXP || typ1 == VECSXP || typ1 == CHARSXP);
    // assert(typ2 == LISTSXP || typ2 == VECSXP || typ2 == CHARSXP);
    // assert(ptr.ptr != nullptr);
    
    std::string input_folder_ = Rcpp::as<std::string>(input_folder);
    std::string out_folder_ = Rcpp::as<std::string>(out_folder);
    std::string extension_ = Rcpp::as<std::string>(extension);
    
    assert(!input_folder_.empty());
    assert(!out_folder_.empty());
    
    // assert(ptr.ptr != nullptr);
    int seq_limit = -1;
    // auto start = std::chrono::high_resolution_clock::now();

    std::filesystem::path outPath(out_folder_);
    std::filesystem::create_directory(outPath);
    if (!std::filesystem::is_empty(outPath))
    {
        // std::cerr << "out_folder : Folder must be empty.";
        // return false;
        Rcpp::stop("out_folder : %s - Folder must be empty.", out_folder_);
        return false;
    }
    // assert(min_batch_size > 0);
    if(min_batch_size == 0){
      min_batch_size = threads;
    }

    std::vector<std::string> inputFiles = getFilesInDir(input_folder_, extension_);

    assert(!inputFiles.empty());

    std::filesystem::path folder = input_folder_;

    assert(inputFiles.size() * inputFiles.size() > 0);
    
    int iterations = inputFiles.size() * inputFiles.size();
    Progress progress_bar(iterations, true);
    // assert(iterations > 0);

    // std::list<std::string> ret_lst(iterations);
    // std::list<std::string> names(iterations);
    for (int i = 0; i < (int)inputFiles.size(); i++)
    {
        for (int j = 0; j < (int)inputFiles.size(); j++)
        {
            Rcpp::checkUserInterrupt();
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
            std::cout << "BLASTing : " << base_name << std::endl;
            static_cast<void>(ptr_->BLAST_files(qry_filePath.string(), subj_filePath.string(), outFileStr, seq_limit, threads, true, false, min_batch_size));
            // ret_lst[one_dim_index] = outFile_.string();
            // names[one_dim_index] = base_name;
        }
    }
    // ret_lst.names() = names;
    PrintClock(start);

    // return rm_null(ret_lst);
    return true;
}



//' @name BLAST2Files
//'
//' @title BLAST 2 Files
//'
//' @description BLAST 2 FASTA files containing nucleotide or protein sequences with a QuickBLAST instance.
//'
//' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
//' @param ptr (Rcpp::XPtr<QuickBLAST>) or (unsigned int) Pointer/ID of QuickBLAST instance
//' @param query (string) Query file
//' @param subject (string) Subject file
//' @param out_file (string) Ouput file (Optional - Intermediate temporary file is created if this option is not provided. Make sure tempdir has enough space)
//' @param seq_limit (int) Batch Size to BLAST at a time. { -1 = Whole File, 0 - One sequence at a time or > 0 } (Optional)
//' @param num_threads (unsigned int) Number of threads. (Optional)
//' @param show_progress (bool) Show progress (Default: TRUE) (Optional)
//' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
//' @param min_batch_size (unsigned int) Minimum batch size (Optional).
//' @return (SEXP) Rcpp::List - if return_values == TRUE, out_file - Otherwise.
//' @export
// [[Rcpp::export]]
RcppExport SEXP BLAST2Files(SEXP ptr, SEXP query, SEXP subject, SEXP out_file = R_NilValue, int seq_limit = -1, unsigned int num_threads = 0, bool show_progress = true, bool return_values = true, unsigned int min_batch_size = 0)
{
    auto start = std::chrono::high_resolution_clock::now();
    unsigned int threads = DetectThreadLimit(num_threads);
    Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
    assert(out_file == R_NilValue || return_values == true);
    assert(out_file != R_NilValue || return_values == false);
    assert(TYPEOF(out_file) == CHARSXP || out_file == R_NilValue);
    assert(TYPEOF(query) == CHARSXP);
    assert(TYPEOF(subject) == CHARSXP);
    // unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    // QuickBLASTHandle ptr = GetQuickBLASTInstance(ptr_id_);

    std::string query_ = Rcpp::as<std::string>(query);
    std::string subject_ = Rcpp::as<std::string>(subject);
    std::string out_file_ = out_file == R_NilValue ? std::tmpnam(nullptr) : Rcpp::as<std::string>(out_file);
  
    if(min_batch_size == 0){
      min_batch_size = threads;
    }    

    assert(!query_.empty());
    assert(!subject_.empty());
    assert(!out_file_.empty()); //|| return_values == true);

    assert(std::filesystem::exists(query_));
    assert(std::filesystem::exists(subject_));
    // ArrowRBVHandle ret_vals;
    std::shared_ptr<arrow::RecordBatchVector> ret_vals;
    if (return_values)
    {
        // std::shared_ptr<arrow::RecordBatchVector> ret_vals = ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_);
        // ret_vals.ptr = ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_); //.get();
        ret_vals = ptr_->BLAST_files(query_, subject_, out_file_, seq_limit, threads, show_progress, return_values, min_batch_size); //.get();
        // PrintClock(start);
        // return ret_vals;
    }
    else
    {
        // static_cast<void>(ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_));
        static_cast<void>(ptr_->BLAST_files(query_, subject_, out_file_, seq_limit, threads, show_progress, return_values, min_batch_size));

        // PrintClock(start);
        // return std::make_shared<arrow::RecordBatchVector>();
        // ret_vals.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
        ret_vals->emplace_back(arrow::RecordBatch::MakeEmpty(ptr_->GetSchema()).ValueOrDie());
    }

    Rcpp::List ret_vals_ = as<List>(Hits2RList(*ret_vals));

    PrintClock(start);
    if(return_values){
      return rm_null(ret_vals_);
    }
}


//' @name RemoteBLAST
 //'
 //' @title BLAST query against remote NCBI DBs
 //'
 //' @description BLAST the input query against remote NCBI DBs (one sequence at a time - to respect rate limits)
 //'
 //' @seealso  [QuickBLAST::GetInstanceID()], [QuickBLAST::GetQuickBLASTInstance()], [QuickBLAST::BLAST2Files()], [QuickBLAST::BLAST2Seqs()], [QuickBLAST::BLAST2Folders()], [QuickBLAST::BLAST1Folder()], [QuickBLAST::RemoteBLAST()]
 //' @param ptr (Rcpp::XPtr<QuickBLAST>) or (unsigned int) Pointer/ID of QuickBLAST instance
 //' @param database (string) Name of the remote NCBI DB
 //' @param query_input (Rcpp::List) (Named) List of input queries (Sequences, Files, Folders - type is determined by input_type parameter)
 //' @param input_type (QuickBLAST::EInputType) Input type (Check [QuickBLAST::GetQuickBLASTEnums()])
 //' @param outFile (string) Output file name (Optional)
 //' @param num_threads (unsigned int) Number of threads. (Optional)
 //' @param return_values (bool) Return BLAST Hits as Rcpp::List (Default: TRUE) (Optional)
 //' @return (SEXP) Rcpp::List - if return_values == TRUE, outFile - Otherwise.
 //' @export
 // [[Rcpp::export]]
 RcppExport SEXP RemoteBLAST(SEXP ptr, SEXP database, SEXP query_input, int input_type, SEXP outFile = R_NilValue, bool return_values = true)
 {
      auto start = std::chrono::high_resolution_clock::now();
      
      assert(outFile == R_NilValue || return_values == true);
      assert(outFile != R_NilValue || return_values == false);
      assert(TYPEOF(outFile) == CHARSXP || outFile == R_NilValue);
      assert(TYPEOF(query_input) == CHARSXP);
      // assert(TYPEOF(program) == CHARSXP);
      
      Rcpp::XPtr<QuickBLAST> ptr_ = ResolveQuickBLASTInstance(ptr);
      std::string program_ = ptr_->GetProgram();//Rcpp::as<std::string>(program);
      std::string database_ = Rcpp::as<std::string>(database);
      Rcpp::List query_input_ = Rcpp::as<Rcpp::List>(query_input);
      QuickBLAST::EInputType input_type_ = static_cast<QuickBLAST::EInputType>(input_type); //Rcpp::as<int>(input_type);
      std::string outFile_ = outFile == R_NilValue ? std::tmpnam(nullptr) : Rcpp::as<std::string>(outFile);
      
      std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
      
      switch(input_type_){
      case QuickBLAST::EInputType::eFile: {
        break;
      }
      case QuickBLAST::EInputType::eSequenceString:{
        std::shared_ptr<arrow::RecordBatchVector> ret_rb = ptr_->BLAST_remote(program_, database_, query_input_, input_type_, outFile_, return_values, 120, 4000);
        // ret_val->emplace_back(ret_rb);
        break;
      }
      case QuickBLAST::EInputType::eFolder:{
        break;
      }
      }
      
      Rcpp::List ret_vals_ = Rcpp::as<List>(Hits2RList(*ret_val));
      
      PrintClock(start);
      if(return_values){
        return rm_null(ret_vals_);
      }
 }