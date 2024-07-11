#include <RcppCommon.h>
#include <Rcpp.h>
#include <algo/blast/QuickBLAST/commons.hpp>
// #include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
// #include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/API.hpp>

using namespace Rcpp;

// #if defined(MINGW32) || defined(WIN32)
extern "C"
{
    //     QBLIBRARY_API BOOL QuickBLASTcpp_main(std::string dllPath){
    //         //std::string dllPath = "inst/x64/QuickBLASTcpp.dll";

    //         // Load the DLL using LoadLibraryEx
    //         HMODULE hDLL = LoadLibraryExA(dllPath.data(), NULL, 0);

    //         if (hDLL == NULL) {
    //             std::cout << "Failed to load the DLL: " << dllPath<< std::endl;  // Handle DLL loading failure
    //             return FALSE;
    //         }

    //         //FreeLibrary(hDLL);
    //         return TRUE;
    //     }

    QuickBLASTHandle GetQuickBLASTInstance(int id);
    QBLIBRARY_API int GetInstanceCount();
    int GetInstanceID(QuickBLASTHandle ptr);
    QBLIBRARY_API SEXP CreateQuickBLASTInstance(SEXP seq_type, SEXP strand, SEXP program, SEXP options, SEXP save_sequences);
    QBLIBRARY_API SEXP DeleteQuickBLASTInstance(SEXP ptr_id);
    //     //QBLIBRARY_API ArrowWrapper* CreateArrowWrapperInstance();
    //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
    //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
    QBLIBRARY_API SEXP SetQuickBLASTOptions(SEXP ptr_id, SEXP program_name, SEXP options);
    QBLIBRARY_API SEXP BLAST2Seqs(SEXP ptr_id, SEXP query, SEXP subject);
    QBLIBRARY_API SEXP BLAST2Files(SEXP ptr_id, SEXP query, SEXP subject, SEXP out_file, SEXP seq_limit, SEXP num_threads, SEXP show_progress, SEXP return_values, SEXP min_batch_size);
    // QBLIBRARY_API SEXP BLAST2Folders(int ptr_id, std::string query, std::string subject, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
    // QBLIBRARY_API SEXP BLAST1Folder(int ptr_id, std::string input_folder, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
    // // QBLIBRARY_API std::string getFilenameWithoutExtension(const std::string &filename);
    // // QBLIBRARY_API Rcpp::List rm_null(Rcpp::List x);
    // // QBLIBRARY_API std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension);
    QBLIBRARY_API void test_QBcpp();
    //     //QBLIBRARY_API  bool QueryOOBESupport() { return false; }
    //     /*QBLIBRARY_API int arrow_struct_num_fields(std::shared_ptr<arrow::StructArray> arr);
    //     QBLIBRARY_API std::shared_ptr<arrow::Array> arrow_struct_field(std::shared_ptr<arrow::StructArray> arr, const int i);
    //     QBLIBRARY_API std::shared_ptr<arrow::Array> arrow_array_slice(std::shared_ptr<arrow::Array> arr, const int offset, const int length);
    //     QBLIBRARY_API int arrow_schema_num_fields(std::shared_ptr<arrow::Schema> sch);
    //     QBLIBRARY_API std::shared_ptr<arrow::Field> arrow_schema_field(std::shared_ptr<arrow::Schema> sch, const int i);
    //     QBLIBRARY_API std::shared_ptr<arrow::DataType> arrow_schema_field_type(std::shared_ptr<arrow::Schema> sch, const int i);
    //     QBLIBRARY_API std::string arrow_schema_field_name(std::shared_ptr<arrow::Schema> sch, const int i);
    //     QBLIBRARY_API bool arrow_int8array_isvalid(std::shared_ptr<arrow::Int8Array> arr, int i);
    //     QBLIBRARY_API bool arrow_strarray_isvalid(std::shared_ptr<arrow::StringArray> arr, int i);
    //     QBLIBRARY_API bool arrow_dblarray_isvalid(std::shared_ptr<arrow::DoubleArray> arr, int i);*/
} // extern "C"

// #endif

void test_QBcpp()
{
    std::cout << "Hello from QuickBLASTcpp" << std::endl;
}

void PrintClock(std::chrono::time_point<std::chrono::high_resolution_clock> start)
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    // Print the time in seconds
    std::cout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;
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

QuickBLASTHandle GetQuickBLASTInstance(int id)
{
    QuickBLASTHandle handle;
    // handle.ptr = obj_list[id].get();
    handle.ptr = obj_list[id];
    assert(handle.ptr != nullptr);
    handle.id = id;
    return handle;
}

int GetInstanceCount()
{
    return (int)obj_list.size();
}

int GetInstanceID(QuickBLASTHandle ptr)
{
    return ptr.ptr->obj_id;
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
        Rcpp::Rcerr << "Only named list or string allowed for BLAST options : Check QuickBLAST::GetAvailableBLASTOptions()";
        // return Rcpp::wrap(false);
        break;
    }
    return options_;
}

SEXP CreateQuickBLASTInstance(SEXP seq_type, SEXP strand, SEXP program, SEXP options, SEXP save_sequences)
{
    // return std::make_shared<QuickBLAST>(seq_type, strand, program, options, save_sequences);
    //  new QuickBLAST(seq_type, strand, program, options, save_sequences);

    QuickBLAST::ESeqType seq_type_ = static_cast<QuickBLAST::ESeqType>(as<int>(seq_type));
    QuickBLAST::EStrand strand_ = static_cast<QuickBLAST::EStrand>(as<int>(strand));
    std::string program_ = as<std::string>(program);

    std::string options_ = ConvertBLASTOptions2String(options);

    // std::string options_ = as<std::string>(options);
    bool save_sequences_ = as<bool>(save_sequences);

    QuickBLASTHandle handle;
    std::shared_ptr<QuickBLAST> objPtr = std::make_shared<QuickBLAST>(seq_type_, strand_, program_, options_, save_sequences_);
    unsigned int list_size = obj_list.size();
    // handle.id = list_size;
    // // handle.ptr = objPtr.get();
    // handle.ptr = objPtr;
    obj_list.insert(std::make_pair(list_size, objPtr));
    return Rcpp::wrap(list_size);
}

SEXP DeleteQuickBLASTInstance(SEXP ptr_id)
{
    unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    //// ptr.reset();
    // fetch the object from the map
    GetQuickBLASTInstance(ptr_id_).ptr->~QuickBLAST();
    obj_list.erase(ptr_id_);
    // ptr.ptr->~QuickBLAST();
    //// delete ptr.ptr.get();
    return Rcpp::wrap(true);
}

SEXP SetQuickBLASTOptions(SEXP ptr_id, SEXP program_name, SEXP options)
{
    unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    std::string program_name_ = as<std::string>(program_name);
    std::string options_ = ConvertBLASTOptions2String(options);

    auto ptr = GetQuickBLASTInstance(ptr_id_);
    ptr.ptr->GetQuickBLASTOptions() = *ptr.ptr->SetQuickBLASTOptions(program_name_, options_);
    return Rcpp::wrap(true);
}

// #ifndef RCPP
SEXP BLAST2Seqs(SEXP ptr_id, SEXP query, SEXP subject)
{
    auto start = std::chrono::high_resolution_clock::now();

    unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    QuickBLASTHandle ptr = GetQuickBLASTInstance(ptr_id_);
    std::string query_ = as<std::string>(query);
    std::string subject_ = as<std::string>(subject);

    assert(!query_.empty());
    assert(!subject_.empty());
    ArrowRBVHandle ret_val;
    std::shared_ptr<arrow::RecordBatch> ret_rb = ptr.ptr->BLAST_seqs(query_, subject_); //
    // ret_val.ptr.emplace_back(ptr.ptr->BLAST_seqs(query_, subject_)); //

    if (ret_rb)
    {
        arrow::Status rb_sts = ret_rb->ValidateFull();
        if (!rb_sts.ok())
        {
            std::cerr << "ERR : Invalid RB : " << rb_sts.message() << rb_sts.detail() << std::endl;
            // PrintClock(start);
            // ret_val.ptr.clear();
            ret_val.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
            // return arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie();
        }
        else
        {
            ret_val.ptr->emplace_back(ret_rb); //
        }
    }
    else
    {
        // PrintClock(start);
        // ret_val.ptr.clear();
        ret_val.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
        // return arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie();
    }

    Rcpp::List ret_vals_ = as<List>(Hits2RList(*ret_val.ptr));

    PrintClock(start);

    return rm_null(ret_vals_);
}

/*

bool BLAST2Folders(QuickBLASTHandle ptr, std::string query, std::string subject, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size)
{
    // int typ1 = TYPEOF(query);
    // int typ2 = TYPEOF(subject);

    // assert(typ1 == LISTSXP || typ1 == VECSXP);
    // assert(typ2 == LISTSXP || typ2 == VECSXP);
    assert(ptr.ptr != nullptr);
    int seq_limit = -1;
    auto start = std::chrono::high_resolution_clock::now();

    std::filesystem::path outPath(out_folder);
    std::filesystem::create_directory(outPath);
    if (!std::filesystem::is_empty(outPath))
    {
        std::cerr << "out_folder : Folder must be empty.";
        return false;
    }
    assert(min_batch_size > 0);
    assert(!query.empty());
    assert(!subject.empty());

    std::filesystem::path query_path(query);
    std::filesystem::path subject_path(subject);
    std::vector<std::string> queryFiles = getFilesInDir(query, extension);
    std::vector<std::string> subjectFiles = getFilesInDir(subject, extension);

    assert(!queryFiles.empty());
    assert(!subjectFiles.empty());

    int iterations = queryFiles.size() * subjectFiles.size();

    assert(iterations > 0);

    // Rcpp::List ret_lst(iterations);
    // Rcpp::CharacterVector names(iterations);
    //  Progress progress_bar(iterations, true);

    for (int i = 0; i < (int)queryFiles.size(); i++)
    {
        for (int j = 0; j < (int)subjectFiles.size(); j++)
        {
            // progress_bar.increment();
            if (queryFiles[i] == subjectFiles[j])
            {
                continue;
            }

            std::string base_name = getFilenameWithoutExtension(queryFiles[i]) + "-" + getFilenameWithoutExtension(subjectFiles[j]);
            std::filesystem::path outFile_ = outPath / (base_name + ".hits");
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
            static_cast<void>(ptr.ptr->BLAST_files(query_input.string(), subject_input.string(), outFile_.string(), seq_limit, num_threads, true, false, min_batch_size));

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

bool BLAST1Folder(QuickBLASTHandle ptr, std::string input_folder, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size)
{

    // int typ1 = TYPEOF(input_folder);

    // assert(typ1 == LISTSXP || typ1 == VECSXP);

    assert(ptr.ptr != nullptr);
    int seq_limit = -1;
    auto start = std::chrono::high_resolution_clock::now();

    std::filesystem::path outPath(out_folder);
    std::filesystem::create_directory(outPath);
    if (!std::filesystem::is_empty(outPath))
    {
        std::cerr << "out_folder : Folder must be empty.";
        return false;
    }
    assert(min_batch_size > 0);
    assert(!input_folder.empty());

    std::vector<std::string> inputFiles = getFilesInDir(input_folder, extension);

    assert(!inputFiles.empty());

    std::filesystem::path folder = input_folder;

    int iterations = inputFiles.size() * inputFiles.size();

    assert(iterations > 0);

    // std::list<std::string> ret_lst(iterations);
    // std::list<std::string> names(iterations);
    for (int i = 0; i < (int)inputFiles.size(); i++)
    {
        for (int j = 0; j < (int)inputFiles.size(); j++)
        {

            if (inputFiles[i] == inputFiles[j])
            {
                continue;
            }

            std::string base_name = getFilenameWithoutExtension(inputFiles[i]) + "-" + getFilenameWithoutExtension(inputFiles[j]);
            std::filesystem::path outFile_ = outPath / (base_name + ".hits");

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
            static_cast<void>(ptr.ptr->BLAST_files(qry_filePath.string(), subj_filePath.string(), outFile_.string(), seq_limit, num_threads, true, false, min_batch_size));
            // ret_lst[one_dim_index] = outFile_.string();
            // names[one_dim_index] = base_name;
        }
    }
    // ret_lst.names() = names;
    PrintClock(start);

    // return rm_null(ret_lst);
    return true;
}

*/

SEXP BLAST2Files(SEXP ptr_id, SEXP query, SEXP subject, SEXP out_file, SEXP seq_limit, SEXP num_threads, SEXP show_progress, SEXP return_values, SEXP min_batch_size)
{
    auto start = std::chrono::high_resolution_clock::now();

    unsigned int ptr_id_ = as<unsigned int>(ptr_id);
    QuickBLASTHandle ptr = GetQuickBLASTInstance(ptr_id_);
    std::string query_ = as<std::string>(query);
    std::string subject_ = as<std::string>(subject);
    std::string out_file_ = as<std::string>(out_file);
    unsigned int seq_limit_ = as<unsigned int>(seq_limit);
    unsigned int num_threads_ = as<unsigned int>(num_threads);
    unsigned int min_batch_size_ = as<unsigned int>(min_batch_size);
    bool show_progress_ = as<bool>(show_progress);
    bool return_values_ = as<bool>(return_values);

    assert(min_batch_size_ > 0);
    assert(!query_.empty());
    assert(!subject_.empty());
    assert(!out_file_.empty());

    assert(std::filesystem::exists(query_));
    assert(std::filesystem::exists(subject_));
    ArrowRBVHandle ret_vals;
    if (return_values_)
    {
        // std::shared_ptr<arrow::RecordBatchVector> ret_vals = ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_);
        ret_vals.ptr = ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_); //.get();
                                                                                                                                                     // PrintClock(start);
                                                                                                                                                     // return ret_vals;
    }
    else
    {
        static_cast<void>(ptr.ptr->BLAST_files(query_, subject_, out_file_, seq_limit_, num_threads_, show_progress_, return_values_, min_batch_size_));
        // PrintClock(start);
        // return std::make_shared<arrow::RecordBatchVector>();
        ret_vals.ptr->emplace_back(arrow::RecordBatch::MakeEmpty(ptr.ptr->GetSchema()).ValueOrDie());
    }

    Rcpp::List ret_vals_ = as<List>(Hits2RList(*ret_vals.ptr));

    PrintClock(start);

    return rm_null(ret_vals_);
}
