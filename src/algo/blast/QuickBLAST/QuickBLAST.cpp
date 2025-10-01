#include <chrono>
#include <iostream>
#include <string_view>
#include <map>
#include <tuple>
#include <future>
#include <iomanip>
#include <sys/mman.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unistd.h>
#include <cassert>
#include <thread>
#include <filesystem>

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

// #include <ArrowWrapper.cpp>
// #include <algo/blast/QuickBLAST/quick_blast-functions.hpp>
// #include <algo/blast/QuickBLAST/arrow_wrapper-functions.hpp>

/* void QuickBLAST::PrintProgressBar(int current, int total, int barWidth)
{
  Rcpp::Rcout << "\033[2J\033[1;1H";
  float progress = static_cast<float>(current) / total;
  int filledWidth = static_cast<int>(progress * barWidth);

  Rcpp::Rcout << "[" << std::setw(barWidth) << std::left << std::string(filledWidth, '=') << std::right << "] ";
  Rcpp::Rcout << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
  Rcpp::Rcout << "  (" << current << " / " << total << ")";
  Rcpp::Rcout << std::flush;
} */

// struct QuickBLAST::Impl {

//     int num_threads = 4;
//     std::string_view run_name;
//     unsigned int obj_id;

//     std::string program;
//     ncbi::CRef<ncbi::blast::CBlastOptionsHandle> opts;
//     // Rcpp::List blast_options_list;
//     // std::string blast_options_str;
//     // SEXP blast_options;
//     std::string blast_options;
//     ESeqType seq_type;
//     EStrand strand;
//     std::shared_ptr<ArrowWrapper> arrow_wrapper;
//     // ArrowWrapper arrow_wrapper;
//     // Rcpp::XPtr<ArrowWrapper> arrow_wrapper;
//     int hit_count = 0;
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//     omp_lock_t hit_countLock;
// #endif
//     bool save_sequences = false;
//     int blast_sequence_limit = 1000;
//     // bool db_scan_mode = false;
//     // std::promise<arrow::Status> ok_promise;

//     std::shared_ptr<ArrowWrapper> GetArrowWrapper(){ return arrow_wrapper; }

//     std::shared_ptr<arrow::Schema> GetSchema() { return arrow_wrapper->GetSchema(); };
//     void SetThreadCount(int num_threads)
//     {
//         this->num_threads = num_threads;
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//         omp_set_num_threads(num_threads);
// #endif
//         arrow_wrapper->SetThreadCount(num_threads);
//     }
//     int GetHitCount()
//     {
//         return hit_count;
//     }
//     void AddHitCount(int val = 1)
//     {
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//         omp_set_lock(&hit_countLock);
// #endif
//         hit_count += val;
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//         omp_unset_lock(&hit_countLock);
// #endif
//     }
//     ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions()
//     {
//         return *opts;
//     }
//     void ResetHitCount() { hit_count = 0; }
//     unsigned int GetObjectID(){
//       return obj_id;
//     }
//     void SetObjectID(unsigned int id){
//         obj_id = id;
//     }
// // QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, Rcpp::List options, bool save_sequences)
// // {
// // #ifdef _OPENMP
// //   this->num_threads = omp_get_num_threads();
// // #else
// //   this->num_threads = 1;
// // #endif
// //   arrow_wrapper = std::make_shared<ArrowWrapper>();
// //   this->save_sequences = save_sequences;
// //   this->program = program;
// //   this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<Rcpp::List>(program, options));
// //   // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
// //   this->strand = strand;
// //   this->seq_type = seq_type;
// //   ok_promise.set_value(arrow::Status::OK());
// // #ifdef _OPENMP
// //   omp_init_lock(&hit_countLock);
// // #endif
// // }

// // QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences)
// Impl(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences)
//         : seq_type(seq_type), strand(strand), program(program), blast_options(options), save_sequences(save_sequences)
// {
//   try
//   {
//     Rprintf("DBG0.1 QB \n");
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//     this->num_threads = omp_get_num_threads();
// #else
//     this->num_threads = 1;
// #endif
//     // SetThreadCount(omp_get_max_threads());
//     // SetThreadCount(std::thread::hardware_concurrency());
//     Rprintf("DBG1 QB \n");
//     arrow_wrapper = std::make_shared<ArrowWrapper>();
//     // arrow_wrapper = ArrowWrapper();
//     // ArrowWrapper *arrow_ptr = new ArrowWrapper();
//     // Rcpp::XPtr<ArrowWrapper> arrow_ptr_(arrow_ptr, true);
//     // arrow_wrapper.reset(new ArrowWrapper(), true);
//     this->save_sequences = save_sequences;
//     this->program = program;
//     Rprintf("DBG2 QB \n");
//     this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<std::string>(program, options));
//     // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
//     this->strand = strand;
//     this->seq_type = seq_type;
//     // ok_promise.set_value(arrow::Status::OK());
//     Rprintf("DBG3 QB \n");
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//     omp_init_lock(&hit_countLock);
// #endif
//   }
//   catch (const std::exception &e)
//   {
//     Rprintf("C++ side Exception: QB : %s\n", e.what());
//   }
// }

// // QuickBLAST::~QuickBLAST()
// ~Impl()
// {
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//   omp_destroy_lock(&hit_countLock);
// #endif

//   // DO NOT DELETE NCBI C++ OBJECTs or PTRs or face Corruption
//   //  delete self;
//   //  opts->ReleaseReference();
//   // delete opts;

//   // delete arrow_wrapper;

//   cout << "~QuickBLAST " << std::endl;
// }

// // // Function to process a single FASTA block
// // void QuickBLAST::PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream)
// // {
// //   if (outputStream != nullptr)
// //   {

// //     // Print FastaSequenceData object
// //     (*outputStream) << "No: " << data->rec_no << std::endl;
// //     (*outputStream) << "Header: " << data->header << std::endl;
// //     (*outputStream) << "Sequence: " << data->seq << std::endl;
// //     (*outputStream) << std::endl;
// //     outputStream->flush();
// //   }
// // }

// // int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
// int GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
// {
//   int frame = 0;
//   if (strand == eNa_strand_plus)
//   {
//     frame = (start % 3) + 1;
//   }
//   else if (strand == eNa_strand_minus)
//   {
//     frame = -(((int)length - start - 1) % 3 + 1);
//   }
//   return frame;
// }

// // std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data)
// std::shared_ptr<arrow::RecordBatch> ExtractFASTA(const FastaSequenceData &fasta_data)
// {
//   // Rcpp::checkUserInterrupt();
//   std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr;
//   std::shared_ptr<arrow::Int32Builder> rec_no_builder;
//   std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
//   rec_no_builder = std::make_shared<arrow::Int32Builder>();
//   fasta_seq_builder = std::make_shared<arrow::StringBuilder>();
//   fasta_h_builder = std::make_shared<arrow::StringBuilder>();
//   static_cast<void>(rec_no_builder->Append(fasta_data.rec_no));
//   static_cast<void>(fasta_h_builder->Append(fasta_data.header));
//   static_cast<void>(fasta_seq_builder->Append(fasta_data.seq));
//   static_cast<void>(fasta_seq_builder->Finish(&seqArr));
//   static_cast<void>(fasta_h_builder->Finish(&hArr));
//   static_cast<void>(rec_no_builder->Finish(&recnoArr));
//   return arrow::RecordBatch::Make(arrow_wrapper->GetFASTASchema(), 1, {recnoArr, hArr, seqArr});
// }

// // std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc)
// std::string GetSSeqLocSequence(const SSeqLoc &seq_loc)
// {
//   const CSeq_id &id = *(seq_loc.seqloc->GetId());

//   // Get the Bioseq using the Seq-id.
//   CBioseq_Handle bioseq_handle = seq_loc.scope->GetBioseqHandle(id);

//   // Terminate the program if the GI cannot be resolved to a Bioseq.
//   if (!bioseq_handle)
//   {
//     ERR_POST(Fatal << "Bioseq not found");
//   }

//   // Get the sequence using CSeqVector.
//   // Use Iupac encoding: CSeq_data::e_Iupacna or CSeq_data::e_Iupacaa.
//   // const auto &length = bioseq_handle.GetBioseqLength();
//   const auto &seq_vect_begin = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).begin();
//   const auto &seq_vect_end = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).end();

//   std::string str(seq_vect_begin, seq_vect_end);

//   return NStr::PrintableString(str);
// }

// // SEXP QuickBLAST::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
// SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
// {

//   // Dispatch based on the data type of the array
//   if (type->id() == arrow::Type::STRUCT)
//   {

//     auto struct_array = std::static_pointer_cast<arrow::StructArray>(array);
//     int num_fields = struct_array->num_fields();

//     // Create an Rcpp list to hold the data frames representing each field of the struct
//     Rcpp::List struct_list(num_fields);
//     Rcpp::CharacterVector names(num_fields);

//     for (int i = 0; i < num_fields; i++)
//     {

//       auto field_array = struct_array->field(i);
//       auto field_type = type->field(i)->type();
//       auto field_name = type->field(i)->name();
//       names[i] = field_name;
//       struct_list[i] = Hits2RList_internal(field_array, field_type, field_name);
//     }

//     struct_list.names() = names;

//     return struct_list;
//   }
//   else if (type->id() == arrow::Type::LIST)
//   {

//     auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
//     auto value_type = type->field(0)->type();

//     // Convert the list array to an Rcpp list
//     Rcpp::List list_values(list_array->length());
//     Rcpp::CharacterVector names(list_array->length());

//     for (int i = 0; i < list_array->length(); i++)
//     {

//       auto sublist_array = list_array->values()->Slice(list_array->value_offset(i), list_array->value_length(i));

//       names[i] = field_name + "[" + std::to_string(i) + "]";
//       auto sublist_name = field_name + "[" + std::to_string(i) + "]";
//       list_values[i] = Hits2RList_internal(sublist_array, value_type, sublist_name);
//     }

//     list_values.names() = names;

//     return list_values;
//   }
//   else if (type->id() == arrow::Type::STRING || type->id() == arrow::Type::LARGE_STRING)
//   {

//     auto string_array = std::static_pointer_cast<arrow::StringArray>(array);

//     Rcpp::StringVector strings(string_array->length());

//     for (int i = 0; i < string_array->length(); ++i)
//     {

//       if (string_array->IsValid(i))
//       {
//         strings[i] = Rcpp::String(string_array->GetString(i));
//       }
//       // else
//       // {
//       //   strings[i] = NA_STRING;
//       // }
//     }

//     return strings;
//   }
//   else if (type->id() == arrow::Type::INT8)
//   {

//     auto int_array = std::static_pointer_cast<arrow::Int8Array>(array);

//     Rcpp::IntegerVector ints(int_array->length());

//     for (int i = 0; i < int_array->length(); ++i)
//     {

//       if (int_array->IsValid(i))
//       {
//         ints[i] = int_array->Value(i);
//       }
//       // else
//       // {
//       //   ints[i] = NA_INTEGER;
//       // }
//     }

//     return ints;
//   }
//   else if (type->id() == arrow::Type::DOUBLE)
//   { // Use arrow::Type::DOUBLE instead of FLOAT64

//     auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
//     Rcpp::NumericVector doubles(double_array->length());

//     for (int i = 0; i < double_array->length(); ++i)
//     {

//       if (double_array->IsValid(i))
//       {
//         doubles[field_name] = double_array->Value(i);
//       }
//       // else
//       // {
//       //   doubles[i] = NA_REAL;
//       // }
//     }

//     return doubles;
//   }
//   else
//   {
//     // For other data types that don't have a direct conversion, return R_NilValue (NA)
//     return R_NilValue;
//   }
// }

// // SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
// SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
// {
//   // Assuming the schema of the RecordBatch is accessible here

//   auto rb_schema = rb->schema();

//   // Convert each column of the RecordBatch to R objects and store in a list
//   Rcpp::List result_list(rb_schema->num_fields());

//   for (int i = 0; i < rb_schema->num_fields(); ++i)
//   {

//     auto array = rb->column(i);
//     auto field_type = rb_schema->field(i)->type();
//     auto field_name = rb_schema->field(i)->name();
//     result_list[i] = Hits2RList_internal(array, field_type, field_name);
//   }

//   return result_list;
// }

// // SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
// SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector)
// {
//   Rcpp::List result_list(rb_vector.size());

//   // Traverse the vector of RecordBatches and convert each RecordBatch
//   for (size_t i = 0; i < rb_vector.size(); ++i)
//   {
//     std::shared_ptr<arrow::RecordBatch> rb = rb_vector[i];
//     result_list[i] = Hits2RList(rb);
//   }

//   return result_list;
// }

// // std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input)
// std::vector<std::pair<std::string, std::string>> BLASTOptionsFromString(const std::string &input)
// {
//   std::vector<std::pair<std::string, std::string>> keyValuePairs;
//   std::istringstream iss(input);
//   std::string token;

//   while (iss >> token)
//   {
//     if (token[0] == '-')
//     {
//       // Extract key-value pair
//       std::string key = token.substr(1);
//       std::string value;

//       if (iss >> value)
//       {
//         keyValuePairs.emplace_back(key, value);
//       }
//       else
//       {
//         // Handle error: Missing value for key
//         cerr << "Error: Missing value for key '" << key << "'." << std::endl;
//         break;
//       }
//     }
//     else
//     {
//       // Handle error: Invalid token (not starting with '-')
//       cerr << "Error: Invalid token '" << token << "'." << std::endl;
//       break;
//     }
//   }

//   return keyValuePairs;
// }

// // template <typename T1>
// // std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
// template <typename T1>
// std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
// {

//   if constexpr (!std::is_same_v<T1, std::string> || !std::is_same_v<T1, FastaSequenceData>)
//   {
//     static_assert(std::is_same_v<T1, T1>, "Unsupported type, only std::string & FastaSequenceData are supported");
//   }

//   return arrow_wrapper->SplitFilesIntoEntries<T1>(filename, delim, num_threads, Entry_callback, return_values);
// }

// // template <typename OptionsType>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
// template <typename OptionsType>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
// {
//   this->blast_options = options;
//   this->program = program_name;
//   if constexpr (std::is_same_v<OptionsType, std::string> || std::is_same_v<OptionsType, Rcpp::String>)
//   {
//     // this->blast_options_str = options;
//     return SetQuickBLASTOptions<std::string>(program_name, options);
//   }
//   else if constexpr (std::is_same_v<OptionsType, Rcpp::List>)
//   {
//     // this->blast_options_list = options;
//     return SetQuickBLASTOptions<Rcpp::List>(program_name, options);
//   }
//   else
//   {
//     static_assert(std::is_same_v<OptionsType, OptionsType>, "Unsupported type");
//   }
// }

// // template <>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
// template <>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
// {
//   assert(!program_name.empty());

//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options.empty())
//   {
//     cout << "Using " << program_name << " Defaults..." << std::endl;
//     return opts;
//   }

//   std::vector<std::pair<std::string, std::string>> keyValuePairs = BLASTOptionsFromString(options);

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (const auto &pair : keyValuePairs)
//   {
//     std::string key_str = pair.second;
//     std::size_t key = std::hash<std::string>{}(pair.first);

//     if (key == hashMap["evalue"])
//     {
//       double val = std::stod(key_str);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = std::stod(key_str);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = (key_str == "TRUE" || key_str == "True" || key_str == "true" || key_str == "1");
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = key_str;
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = std::stod(key_str);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = std::stod(key_str);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// }

// // template <>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
// template <>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
// {
//   assert(!program_name.empty());
//   Rcpp::List options_(options);
//   for (int i = 0; i < options_.size(); ++i)
//   {
//     // Check if the element is empty
//     if (Rf_isNull(options_[i]))
//     {
//       // Remove the empty element
//       options_.erase(i);
//       --i; // Decrement the index since the list size has changed
//     }
//   }
//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options_.size() == 0 || options_.isNULL())
//   {
//     // Rcpp::Rcout << "Using " << program_name << " Defaults..." << std::endl;
//     Rprintf("Using %s Defaults...\n", program_name);
//     return opts;
//   }

//   Rcpp::CharacterVector keys = options_.names();

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (int i = 0; i < keys.size(); ++i)
//   {
//     std::string key_str = Rcpp::as<std::string>(keys[i]);
//     std::size_t key = std::hash<std::string>{}(key_str);

//     if (key == hashMap["evalue"])
//     {
//       double val = Rcpp::as<double>(options_["evalue"]);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = Rcpp::as<double>(options_["pident"]);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = Rcpp::as<bool>(options_["gapped_mode"]);
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = Rcpp::as<std::string>(options_["filter_string"]);
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = Rcpp::as<int>(options_["effective_search_space"]);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = Rcpp::as<int>(options_["cutoff_score"]);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = Rcpp::as<double>(options_["gap_trigger"]);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff"]);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff_final"]);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = Rcpp::as<int>(options_["hit_list_size"]);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = Rcpp::as<double>(options_["low_score_percentage"]);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_subject"]);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_sequence"]);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = Rcpp::as<double>(options_["qcovhsp_perc"]);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = Rcpp::as<int>(options_["window_size"]);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// }

// template <typename T>
// std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope)
// {
//     if constexpr (std::is_same_v<T, TSeqLocVector>)
//     {
//         arrow::RecordBatchVector recBth_vec;

//         for (const auto &s_it : sloc)
//         {
//             // Rcpp::checkUserInterrupt();

//             std::shared_ptr<arrow::RecordBatch> rb = ExtractHits<SSeqLoc>(alignments, qloc, s_it, scope);

//             if (rb)
//             {
//                 recBth_vec.emplace_back(std::move(rb));
//             }
//         }

//         const auto &wrt_sts = arrow_wrapper->AddRBV2Batch(recBth_vec);
//         if (wrt_sts.ok())
//         {
//             return wrt_sts.ValueOrDie();
//         }
//         else
//         {
//             cerr << "Warn : Invalid Alignment RBV (Returning Empty) : " << wrt_sts.status().detail() << std::endl
//                  << wrt_sts.status().message() << std::endl;
//         }
//         return std::make_shared<arrow::RecordBatchVector>();
//     }
//     // For SSeqLoc
//     else if constexpr (std::is_same_v<T, SSeqLoc>)
//     {
//         std::shared_ptr<arrow::RecordBatch> rb = ExtractHits<SSeqLoc>(alignments, qloc, sloc, scope);

//         const auto &wrt_sts = arrow_wrapper->AddRB2Batch(rb);
//         if (wrt_sts.ok())
//         {
//             return wrt_sts.ValueOrDie();
//         }
//         else
//         {
//             cerr << "Warn : Invalid Alignment RB (Returning Empty) : " << wrt_sts.status().detail() << std::endl
//                  << wrt_sts.status().message() << std::endl;
//         }

//         return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
//     }
//     else
//     {
//         static_assert(std::is_same_v<T, T>, "Unsupported type, only ncbi::blast::TSeqLocVector & ncbi::blast::SSeqLoc are supported");
//     }
// }

// // std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
// std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
// {
//   // assert(num_threads > 0);
//   /*   if (!arrow_wrapper || arrow_wrapper.get() == nullptr)
//     {
//       arrow_wrapper = std::make_shared<ArrowWrapper>();
//     } */
//   /* if (this->opts.Empty() || this->opts.IsNull())
//   {
//     // if (this->blast_options_list.size() > 0)
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_list);
//     // }
//     // else if (!this->blast_options_str.empty())
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_str);
//     // }
//     // else
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, "");
//     // }
//     this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, this->blast_options));
//   } */

// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//   unsigned int n_threads = num_threads > omp_get_num_threads() ? omp_get_num_threads() : num_threads;
// #else
//   unsigned int n_threads = 1;
// #endif

//   n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
//   arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile);
//   if (!outfile_sts.ok())
//   {
//     /* Rcpp::Rcerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
//                 << outfile_sts.message() << std::endl; */
//     // cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
//     //      << outfile_sts.message() << std::endl;
//     REprintf("ERROR : Could not create output file stream : %s \n %s \n", outfile_sts.detail()->ToString().c_str(), outfile_sts.message().c_str());
//     return std::make_shared<arrow::RecordBatchVector>();
//   }

//   SetThreadCount(n_threads);

//   unsigned int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);

//   unsigned int s_seq_count = arrow_wrapper->CountCharacter(subjectFile, '>', n_threads);

//   const unsigned int totalIterations = q_seq_count * s_seq_count;
//   if (blast_sequence_limit > 0)
//   {
//     blast_sequence_limit = blast_sequence_limit > totalIterations ? totalIterations : blast_sequence_limit;
//     blast_sequence_limit = blast_sequence_limit > s_seq_count ? s_seq_count - 1 : blast_sequence_limit; // - 1;
//   }
//   else
//   {
//     blast_sequence_limit = s_seq_count - 1;
//   }
//   assert(totalIterations > 0);

//   //// int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
//   ////  batch_size = 32 * n_threads; // batch_size > 0 ? batch_size : 1024;
//   arrow_wrapper->SetBatchSize(batch_size);

//   // Progress progress_bar(totalIterations, show_progress);

//   std::shared_ptr<arrow::RecordBatchVector> final_ret = StreamFile<FastaSequenceData>(
//       queryFile, ">", n_threads, [this, n_threads, subjectFile, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_q) //&progress_bar
//       {
//             CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
//             const std::shared_ptr<SSeqLoc> query_loc(std::move(CreateSSeqLocFromType<FastaSequenceData>(*data_q, scope)));

//             _ASSERT(query_loc->seqloc.NotEmpty());

//             // Rcpp::checkUserInterrupt();

//             std::shared_ptr<TSeqLocVector> subjects_loc_vec(new TSeqLocVector());

// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//             omp_lock_t query_locLock;
//             omp_lock_t subjects_loc_vecLock;
//             omp_init_lock(&query_locLock);
//             omp_init_lock(&subjects_loc_vecLock);
// #endif

//             std::shared_ptr<arrow::RecordBatchVector> ret_results = StreamFile<FastaSequenceData>(
//                 subjectFile, ">", n_threads, [this, query_loc,
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                 &query_locLock, &subjects_loc_vecLock,
// #endif
//                 & scope, &subjects_loc_vec, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_s) // &progress_bar
//                 {
//                     const std::unique_ptr<SSeqLoc> subject_loc(CreateSSeqLocFromType<FastaSequenceData>(*data_s, scope));
//                     _ASSERT(subject_loc->seqloc.NotEmpty());

//                     if (strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) != 0)
//                     {
//                         // Rcpp::checkUserInterrupt();

//                         CBl2Seq* blaster;

//                         try
//                         {

//                             switch (blast_sequence_limit)
//                             {
//                             case 0:
//                             {
//                                 blaster = new CBl2Seq(*query_loc, *subject_loc, this->GetQuickBLASTOptions());
//                                 arrow::RecordBatchVector tmp_rbv = { ExtractHits<SSeqLoc>(blaster->Run(), *query_loc, *subject_loc, *scope) };
//                                 // progress_bar.increment();
//                                 if (return_values)
//                                 {
//                                     return std::make_shared<arrow::RecordBatchVector>(tmp_rbv);
//                                 }
//                                 else
//                                 {
//                                     tmp_rbv.clear();
//                                     return std::make_shared<arrow::RecordBatchVector>();
//                                 }
//                             }
//                             break;
//                             default:
//                             {
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                                 omp_set_lock(&subjects_loc_vecLock);
// #endif
//                                 subjects_loc_vec->emplace_back(*subject_loc);
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                                 omp_unset_lock(&subjects_loc_vecLock);
// #endif

//                                 // progress_bar.increment();
//                                 if (subjects_loc_vec->size() >= blast_sequence_limit)
//                                 {
//                                     TSeqLocVector subjects_buffer_vec;
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                                     omp_set_lock(&subjects_loc_vecLock);

// #endif
//                                     subjects_buffer_vec.swap(*subjects_loc_vec);
//                                     subjects_loc_vec->clear();
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                                     omp_unset_lock(&subjects_loc_vecLock);
// #endif
//                                     blaster = new CBl2Seq(*query_loc, subjects_buffer_vec, this->GetQuickBLASTOptions(), true);

//                                     AddHitCount(subjects_buffer_vec.size());
//                                     std::shared_ptr<arrow::RecordBatchVector> tmp_rbv = ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
//                                     //return ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
//                                     if (return_values)
//                                     {
//                                         return tmp_rbv;
//                                     }
//                                     else
//                                     {
//                                         tmp_rbv->clear();
//                                         return std::make_shared<arrow::RecordBatchVector>();
//                                     }
//                                 }
//                             }
//                             break;
//                             }
//                         }
//                         catch (const CException& e)
//                         {
//                             // Handle exception ...
//                             cout << e.GetFunction() << std::endl;
//                             cout << e.GetErrCodeString() << std::endl;
//                             cout << e.GetErrCode() << std::endl;
//                             cout << e.GetModule() << std::endl;
//                             cout << e.GetPredecessor() << std::endl;
//                             cout << e.GetFile() << std::endl;
//                             cout << e.GetLine() << std::endl;
//                             cout << e.GetMsg() << std::endl;
//                             cout << e.GetStackTrace() << std::endl;
//                             cout << e.GetStackTraceLevel() << std::endl;
//                             cout << e.GetClass() << std::endl
//                                 << std::flush;
//                         }
//                         catch (const std::exception& e)
//                         {
//                             cout << e.what() << std::endl
//                                 << std::flush;
//                         }
//                     }

//                     return std::make_shared<arrow::RecordBatchVector>(); // EMPTY ERROR Return
//                 },
//                 return_values);

//             if (subjects_loc_vec->size() > 0)
//             {
//                 CBl2Seq blaster(*query_loc, *subjects_loc_vec, this->GetQuickBLASTOptions(), true);
//                 AddHitCount(subjects_loc_vec->size());
//                 std::shared_ptr<arrow::RecordBatchVector> ret_vec = ExtractHits(blaster.Run(), *query_loc, *subjects_loc_vec, *scope);

//                 if (return_values)
//                 {
//                     ret_results->insert(ret_results->end(), ret_vec->begin(), ret_vec->end());
//                 }
//                 else
//                 {
//                     ret_vec->clear();
//                 }
//             }

// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//             omp_destroy_lock(&query_locLock);
//             omp_destroy_lock(&subjects_loc_vecLock);
// #endif
//             if (return_values) {
//                 return ret_results;
//             }
//             else
//             {
//                 ret_results->clear();
//                 return std::make_shared<arrow::RecordBatchVector>();
//             } },
//       return_values);

// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
// #pragma omp barrier
// #endif

//   arrow_wrapper->FinishOutputStream();
//   if (return_values)
//   {
//     return final_ret;
//   }
//   else
//   {
//     final_ret->clear();
//     return std::make_shared<arrow::RecordBatchVector>();
//   }
// }

// // std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject)
// std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject)
// {
//   // Rcpp::checkUserInterrupt();

//   CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));

//   std::unique_ptr<SSeqLoc>
//       query_seqloc(CreateSSeqLocFromType<std::string>(query, scope));
//   std::unique_ptr<SSeqLoc> subject_seqloc(CreateSSeqLocFromType<std::string>(subject, scope));

//   CBl2Seq blaster(*query_seqloc, *subject_seqloc, GetQuickBLASTOptions());

//   return ExtractHits<SSeqLoc>(blaster.Run(), *query_seqloc, *subject_seqloc, *scope);
// }

// //' @name BLAST C++ Call
// //' @title BLAST C++ Call
// //'
// //' @description BLAST 2 Files/Seqs. This is for the QuickBLAST C++ object exposed in R
// //'
// //' @param query (string) Query FASTA File/Seq
// //' @param subject (string) Subject FASTA File/Seq
// //' @param outputFile (string) Output Filename (Arrow Feather/IPC Format)  - Not used for Sequence BLAST
// //' @param input_type - (QuickBLAST::EInputType) 0 - eFile, 1 - eSequenceString
// //' @param blast_sequence_limit (int) Batch Size to BLAST at a time  - Not used for Sequence BLAST
// //' @param show_progress (bool) TRUE - Show progress, Set FALSE for multiple instances - Not used for Sequence BLAST
// //' @return Nested List of BLAST Hits
// // auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
// auto BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
// {

//   assert(std::filesystem::exists(query));
//   assert(std::filesystem::exists(subject));

//   switch (input_type)
//   {
//   case QuickBLAST::EInputType::eFile:
//   {
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//     unsigned int n_threads = omp_get_num_threads();
// #else
//     unsigned int n_threads = 1;
// #endif
//     int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
//     batch_size = 32 * n_threads;
//     batch_size = batch_size > 0 ? batch_size : 1024;
//     return BLAST_files(query, subject, outputFile, blast_sequence_limit, n_threads, show_progress, true, batch_size);
//     // return Hits2RList(*ret_val);
//   }
//   break;
//   case QuickBLAST::EInputType::eSequenceString:
//   {
//     arrow::RecordBatchVector ret_val;
//     ret_val.emplace_back(BLAST_seqs(query, subject));
//     return std::make_shared<arrow::RecordBatchVector>(ret_val);
//     // return Hits2RList(ret_val);
//   }
//   break;
//   default:
//   {
//     // Rcpp::Rcerr << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
//     // cout << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
//     REprintf("input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !");
//     // return false; //Rcpp::wrap(false);
//     return std::make_shared<arrow::RecordBatchVector>();
//   }
//   break;
//   }
// }

// };
/////////END of Impl struct

std::shared_ptr<ArrowWrapper> QuickBLAST::Impl::GetArrowWrapper() { return arrow_wrapper; }

std::shared_ptr<arrow::Schema> QuickBLAST::Impl::GetSchema() { 
  if (!arrow_wrapper) {
    // initialize default options or throw helpful error
    Rcpp::stop("QuickBLAST::GetSchema(): arrow_wrapper is not initialised");
  }
  return arrow_wrapper->GetSchema(); 
};

void QuickBLAST::Impl::SetThreadCount(unsigned int num_threads)
{
  this->num_threads = num_threads;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_set_num_threads(num_threads);
#endif
  arrow_wrapper->SetThreadCount(num_threads);
}
unsigned int QuickBLAST::Impl::GetThreadCount()
{
  return num_threads;
}

int QuickBLAST::Impl::GetHitCount()
{
  return hit_count;
}
void QuickBLAST::Impl::AddHitCount(int val)
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_set_lock(&hit_countLock);
#endif
  hit_count += val;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_unset_lock(&hit_countLock);
#endif
}
ncbi::blast::CBlastOptionsHandle &QuickBLAST::Impl::GetQuickBLASTOptions()
{
  if (!opts) {
    // initialize default options or throw helpful error
    Rcpp::stop("GetQuickBLASTOptions: options handle 'opts' is not initialised");
  }
  return *opts;
}
void QuickBLAST::Impl::ResetHitCount() { hit_count = 0; }
unsigned int QuickBLAST::Impl::GetObjectID()
{
  return obj_id;
}
void QuickBLAST::Impl::SetObjectID(unsigned int id)
{
  obj_id = id;
}

ncbi::blast::CBlastOptionsHandle *QuickBLAST::Impl::SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality = CBlastOptions::eLocal)
{
  assert(!program_name.empty());

  ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
  // Create a CBlastOptionsHandle object
  ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program, locality);
  opts->SetDefaults();
  // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
  // Example: Extracting and setting the BLAST database

  if (options.empty())
  {
    cout << "Using " << program_name << " Defaults..." << std::endl;
    return opts;
  }

  std::vector<std::pair<std::string, std::string>> keyValuePairs = BLASTOptionsFromString(options);

  std::unordered_map<std::string, std::size_t> hashMap;

  hashMap["evalue"] = std::hash<std::string>{}("evalue");
  hashMap["pident"] = std::hash<std::string>{}("pident");
  hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
  hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
  hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
  hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
  hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
  hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
  hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
  hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
  hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
  hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
  hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
  hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
  hashMap["window_size"] = std::hash<std::string>{}("window_size");

  for (const auto &pair : keyValuePairs)
  {
    std::string key_str = pair.second;
    std::size_t key = std::hash<std::string>{}(pair.first);

    if (key == hashMap["evalue"])
    {
      double val = std::stod(key_str);
      opts->SetEvalueThreshold(val);
    }
    else if (key == hashMap["pident"])
    {
      double val = std::stod(key_str);
      opts->SetPercentIdentity(val);
    }
    else if (key == hashMap["gapped_mode"])
    {
      bool val = (key_str == "TRUE" || key_str == "True" || key_str == "true" || key_str == "1");
      opts->SetGappedMode(val);
    }
    else if (key == hashMap["filter_string"])
    {
      std::string val = key_str;
      opts->SetFilterString(val.c_str());
    }
    else if (key == hashMap["effective_search_space"])
    {
      int val = std::stoi(key_str);
      opts->SetEffectiveSearchSpace(val);
    }
    else if (key == hashMap["cutoff_score"])
    {
      int val = std::stoi(key_str);
      opts->SetCutoffScore(val);
    }
    else if (key == hashMap["gap_trigger"])
    {
      double val = std::stod(key_str);
      opts->SetGapTrigger(val);
    }
    else if (key == hashMap["gap_x_dropoff"])
    {
      double val = std::stod(key_str);
      opts->SetGapXDropoff(val);
    }
    else if (key == hashMap["gap_x_dropoff_final"])
    {
      double val = std::stod(key_str);
      opts->SetGapXDropoffFinal(val);
    }
    else if (key == hashMap["hit_list_size"])
    {
      int val = std::stoi(key_str);
      opts->SetHitlistSize(val);
    }
    else if (key == hashMap["low_score_percentage"])
    {
      double val = std::stod(key_str);
      opts->SetLowScorePerc(val);
    }
    else if (key == hashMap["max_hsp_per_subject"])
    {
      int val = std::stoi(key_str);
      opts->SetMaxHspsPerSubject(val);
    }
    else if (key == hashMap["max_hsp_per_sequence"])
    {
      int val = std::stoi(key_str);
      opts->SetMaxNumHspPerSequence(val);
    }
    else if (key == hashMap["qcovhsp_perc"])
    {
      double val = std::stod(key_str);
      opts->SetQueryCovHspPerc(val);
    }
    else if (key == hashMap["window_size"])
    {
      int val = std::stoi(key_str);
      opts->SetWindowSize(val);
    }
  }

  opts->Validate();

  this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(opts);

  blast_options = options;
  
  return opts;
}

// QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, Rcpp::List options, bool save_sequences)
// {
// #ifdef _OPENMP
//   this->num_threads = omp_get_num_threads();
// #else
//   this->num_threads = 1;
// #endif
//   arrow_wrapper = std::make_shared<ArrowWrapper>();
//   this->save_sequences = save_sequences;
//   this->program = program;
//   this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<Rcpp::List>(program, options));
//   // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
//   this->strand = strand;
//   this->seq_type = seq_type;
//   ok_promise.set_value(arrow::Status::OK());
// #ifdef _OPENMP
//   omp_init_lock(&hit_countLock);
// #endif
// }

// QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences)
QuickBLAST::Impl::Impl(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences, bool save_hsp_sequences)
    : seq_type(seq_type), strand(strand), program(program), blast_options(options), save_sequences(save_sequences), save_hsp_sequences(save_hsp_sequences)
{
  try
  {
    // Rprintf("DBG0.1 QB \n");
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    this->num_threads = omp_get_num_threads();
#else
    this->num_threads = 1;
#endif
    // SetThreadCount(omp_get_max_threads());
    // SetThreadCount(std::thread::hardware_concurrency());
    // Rprintf("DBG1 QB \n");
    arrow_wrapper = std::make_shared<ArrowWrapper>();
    const auto empty_rb_ = arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema());
    if(!empty_rb_.ok()){
      throw std::runtime_error("QuickBLAST::Impl(): Error creating empty record batch.");
    }
    empty_rb = empty_rb_.ValueOrDie();
    // arrow_wrapper = ArrowWrapper();
    // ArrowWrapper *arrow_ptr = new ArrowWrapper();
    // Rcpp::XPtr<ArrowWrapper> arrow_ptr_(arrow_ptr, true);
    // arrow_wrapper.reset(new ArrowWrapper(), true);
    this->save_sequences = save_sequences;
    this->save_hsp_sequences = save_hsp_sequences;
    this->program = program;
    // Rprintf("DBG2 QB \n");
    // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<std::string>(program, options));
    this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
    // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
    this->strand = strand;
    this->seq_type = seq_type;
    // ok_promise.set_value(arrow::Status::OK());
    // Rprintf("DBG3 QB \n");
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_init_lock(&hit_countLock);
#endif
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("QuickBLAST::Impl(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("QuickBLAST::Impl(): C++ Runtime Error : ") + e.what());
  }catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("QuickBLAST::Impl(): Rcpp::exception : ") + e.what());
  }
  catch(...){
    Rcpp::stop("QuickBLAST::Impl(): Unknown Exception");
  }
}

// QuickBLAST::~QuickBLAST()
QuickBLAST::Impl::~Impl()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_destroy_lock(&hit_countLock);
#endif

  // DO NOT DELETE NCBI C++ OBJECTs or PTRs or face Corruption
  //  delete self;
  //  opts->ReleaseReference();
  // delete opts;

  // delete arrow_wrapper;

  // Rprintf("~QuickBLAST::Impl ");
}

// // Function to process a single FASTA block
// void QuickBLAST::PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream)
// {
//   if (outputStream != nullptr)
//   {

//     // Print FastaSequenceData object
//     (*outputStream) << "No: " << data->rec_no << std::endl;
//     (*outputStream) << "Header: " << data->header << std::endl;
//     (*outputStream) << "Sequence: " << data->seq << std::endl;
//     (*outputStream) << std::endl;
//     outputStream->flush();
//   }
// }

// int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
int QuickBLAST::Impl::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
{
  int frame = 0;
  if (strand == eNa_strand_plus)
  {
    frame = (start % 3) + 1;
  }
  else if (strand == eNa_strand_minus)
  {
    frame = -(((int)length - start - 1) % 3 + 1);
  }
  return frame;
}

// std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data)
std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractFASTA(const FastaSequenceData &fasta_data)
{
  Rcpp::checkUserInterrupt();
  std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr;
  std::shared_ptr<arrow::Int64Builder> rec_no_builder;
  std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
  rec_no_builder = std::make_shared<arrow::Int64Builder>();
  fasta_seq_builder = std::make_shared<arrow::StringBuilder>();
  fasta_h_builder = std::make_shared<arrow::StringBuilder>();
  static_cast<void>(rec_no_builder->Append(fasta_data.rec_no));
  static_cast<void>(fasta_h_builder->Append(fasta_data.header));
  static_cast<void>(fasta_seq_builder->Append(fasta_data.seq));
  static_cast<void>(fasta_seq_builder->Finish(&seqArr));
  static_cast<void>(fasta_h_builder->Finish(&hArr));
  static_cast<void>(rec_no_builder->Finish(&recnoArr));
  return arrow::RecordBatch::Make(arrow_wrapper->GetFASTASchema(), 1, {recnoArr, hArr, seqArr});
}

static std::string sanitize_protein(const std::string &s_in, bool &had_bad) {
  static const std::string allowed = "ACDEFGHIKLMNPQRSTVWYBXZJUO"; // include common ambiguity letters
  std::string out;
  out.reserve(s_in.size());
  had_bad = false;
  for (unsigned char c0 : s_in) {
    if (c0 == '\0') { had_bad = true; continue; }
    unsigned char c = static_cast<unsigned char>(std::toupper(c0));
    if (c >= 'A' && c <= 'Z') {
      if (allowed.find(c) != std::string::npos) out.push_back(char(c));
      else { had_bad = true; /* optional: map unknown to 'X' */ out.push_back('X'); }
    } else {
      // skip digits, punctuation, whitespace
      had_bad = true;
    }
  }
  return out;
}

static std::string sanitize_nucleotide(const std::string &s_in, bool &had_bad) {
  static const std::string allowed = "ACGTUN*-"; // IUPAC nucleotides
  std::string out; out.reserve(s_in.size()); had_bad = false;
  for (unsigned char c0 : s_in) {
    if (c0 == '\0') { had_bad = true; continue; }
    unsigned char c = static_cast<unsigned char>(std::toupper(c0));
    if (c >= 'A' && c <= 'Z') {
      if (allowed.find(c) != std::string::npos) out.push_back(char(c));
      else { had_bad = true; /* map unknown to N */ out.push_back('N'); }
    } else {
      had_bad = true;
    }
  }
  return out;
}

// SSeqLoc *QuickBLAST::Impl::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope)
std::pair<std::shared_ptr<SSeqLoc>, CSeq_entry_Handle> QuickBLAST::Impl::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope)
{
  //https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/ident
  
  int rec_no = fasta_data.rec_no;
  std::string fastaID(fasta_data.header.data());
  std::string fastaSequence(fasta_data.seq.data());
  
  // Rcpp::Rcout << "Header : "  << rec_no << " : " << fastaID << std::endl << "Sequence : " << fastaSequence << std::endl << std::flush; //DEBUG
  
  // std::string cleanID = fastaID;
  // std::replace(cleanID.begin(), cleanID.end(), ' ', '_');
  // if(cleanID.length() > 50)
  //   cleanID.resize(50); //respecting the 50 char limit for FASTA headers by BLAST
  
  // create a seqdesc containing a title string (cleanID)
  // CRef<ncbi::objects::CSeqdesc> title_desc(new ncbi::objects::CSeqdesc());
  // title_desc->SetTitle(fastaID); //cleanID
  // 
  // // put it into a CSeq_descr and attach to the bioseq
  // CRef<ncbi::objects::CSeq_descr> descr(new ncbi::objects::CSeq_descr());
  // descr->Set().push_back(title_desc);

  const TSeqPos seqlen = fastaSequence.length();

  _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
  // ncbi::CRef<ncbi::objects::CSeq_interval> interval(new ncbi::objects::CSeq_interval());
  // interval->SetFrom(0);
  // interval->SetTo(seqlen - 1);
  // 
  // CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
  // id->Select(CSeq_id_Base::E_Choice::e_Local);
  // id->SetLocal().SetId(rec_no);
  // id->SetLocal().SetStr(fastaID);

  // std::string cleanID = fastaID;
  //   // std::replace(cleanID.begin(), cleanID.end(), ' ', '_');
  //   if(cleanID.length() > 50)
  //     cleanID.resize(50); //respecting the 50 char limit for FASTA headers by BLAST
  CRef<CSeq_id> id(new CSeq_id(fastaID, CSeq_id::fParse_RawText | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal));
  // // id->SetLocal().SetId(rec_no);
  // // id->SetLocal().SetStr(cleanID);
  // id->SetLocal().SetStr(fastaID);
  
  CRef<CSeq_loc>
      cseq_loc_obj(new CSeq_loc());
  cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
  cseq_loc_obj->SetId(*id);
  // cseq_loc_obj->SetWhole()
  //     .SetLocal()
  //     .SetStr(fastaID); //cleanID
  // cseq_loc_obj->SetWhole()
  //   .SetLocal()
  //   .SetId(rec_no);
  if (seq_type == ESeqType::eProtein)
  {
    cseq_loc_obj->SetStrand(eNa_strand_unknown);
  }
  else
  {
    switch (strand)
    {
    case EStrand::ePlus:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_plus);
      break;
    case EStrand::eMinus:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_minus);
      break;
    case EStrand::eUnknown:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_unknown);
      break;
    case EStrand::eBoth:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both);
      break;
    case EStrand::eBoth_rev:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both_rev);
      break;
    case EStrand::eOther:
      cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_other);
      break;
    }
  }

  CRef<CSeq_data> seq_data(new CSeq_data());
  seq_data->Select(seq_type == ESeqType::eProtein ? CSeq_data_Base::E_Choice::e_Iupacaa : CSeq_data_Base::E_Choice::e_Iupacna);
  switch (seq_type)
  {
  case ESeqType::eProtein:
  {
    seq_data->SetIupacaa(CIUPACaa(fastaSequence));
    break;
  }
  case ESeqType::eNucleotide:
  {
    seq_data->SetIupacna(CIUPACna(fastaSequence));
    break;
  }
  }

  CRef<CSeq_inst> seq_inst(new CSeq_inst());
  seq_inst->SetSeq_data(*seq_data);
  seq_inst->SetLength(fastaSequence.length());
  seq_inst->SetRepr(CSeq_inst::eRepr_raw);
  // other useful defaults:
  seq_inst->SetTopology(CSeq_inst::eTopology_linear);
  
  if (seq_type == ESeqType::eProtein)
  {
    // seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
    seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_aa);
  }
  else
  {
    // seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
    seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_na);
  }
  
  // 
  // // enum EStrand {
  // //   eStrand_not_set =   0,
  // //   eStrand_ss      =   1,  ///< single strand
  // //   eStrand_ds      =   2,  ///< double strand
  // //   eStrand_mixed   =   3,
  // //   eStrand_other   = 255  ///< default ds for DNA, ss for RNA, pept
  // // };
  // 
  // seq_inst->SetRepr(CSeq_inst_Base::ERepr::eRepr_raw);
  // seq_inst->SetTopology(CSeq_inst_Base::ETopology::eTopology_linear);
  seq_inst->SetLength(seqlen);
  // seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
  // 
  // // { CSeq_inst::eMol_not_set, " " } ,
  // // { CSeq_inst::eMol_dna,     "DNA" } ,
  // // { CSeq_inst::eMol_rna,     "RNA" } ,
  // // { CSeq_inst::eMol_aa,      "protein" } ,
  // // { CSeq_inst::eMol_na,      "nucleotide" } ,
  // // { CSeq_inst::eMol_other,   "other" }
  
  // CRef<ncbi::objects::CBioseq> bioseq(new CBioseq(*cseq_loc_obj, cleanID)); //fastaID
  // bioseq->SetInst(*seq_inst);
  // 
  // CRef<CSeq_entry>
  //     ret_entry(new CSeq_entry());
  // ret_entry->SetSeq(*bioseq);
  // 
  // parent_scope->AddTopLevelSeqEntry(*ret_entry);

  CRef<CBioseq> bioseq(new CBioseq());
  bioseq->SetId().push_back(id);      // ID that matches cseq_loc above
  bioseq->SetInst(*seq_inst);
  // bioseq->SetDescr(*descr);
  
  // optionally add description/title:
  // CRef<GBSeq> descr(new GBSeq()); ... bioseq->SetDescr(...);
  
  // -- put into an entry and add to scope (so the scope can resolve lcl|rec_no)
  CRef<CSeq_entry> ret_entry(new CSeq_entry());
  ret_entry->SetSeq(*bioseq);
  
  // Add the sequence entry to the scope so Seq-id <-> SeqMap lookups succeed.
  CSeq_entry_Handle tse_handle = parent_scope->AddTopLevelSeqEntry(*ret_entry);
  // parent_scope->AddBioseq(*bioseq);

  // //DEBUG
  // for (auto &id0 : ret_entry->GetSeq().GetId()) {
  //   Rcpp::Rcout << "Added seq id: " << id0->AsFastaString() << std::endl;
  // }
  // try {
  //   CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(*ret_entry->GetSeq().GetId().front());
  //   CBioseq_Handle bh = parent_scope->GetBioseqHandle(idh);
  //   Rcpp::Rcout << "Scope resolves id: " << idh.AsString() << std::endl;
  // } catch (const CException &e) {
  //   Rcpp::Rcerr << "Scope resolution failed after add: " << e.GetMsg() << std::endl;
  // } //DEBUG END
  
  // return new SSeqLoc(cseq_loc_obj.GetObject(), parent_scope.GetObject());
  
  // return new SSeqLoc(cseq_loc_obj, parent_scope);
  
  auto sseq = std::make_shared<ncbi::blast::SSeqLoc>(
       cseq_loc_obj,
       parent_scope
    );
  
  return std::make_pair(std::move(sseq), tse_handle);
}

// std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc)
std::string QuickBLAST::Impl::GetSSeqLocSequence(const SSeqLoc &seq_loc)
{
  const CSeq_id &id = *(seq_loc.seqloc->GetId());

  // Get the Bioseq using the Seq-id.
  CBioseq_Handle bioseq_handle = seq_loc.scope->GetBioseqHandle(id);

  // Terminate the program if the GI cannot be resolved to a Bioseq.
  if (!bioseq_handle)
  {
    ERR_POST(Fatal << "Bioseq not found");
  }

  // Get the sequence using CSeqVector.
  // Use Iupac encoding: CSeq_data::e_Iupacna or CSeq_data::e_Iupacaa.
  // const auto &length = bioseq_handle.GetBioseqLength();
  const auto &seq_vect_begin = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).begin();
  const auto &seq_vect_end = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).end();

  std::string str(seq_vect_begin, seq_vect_end);

  return NStr::PrintableString(str);
}

// SEXP QuickBLAST::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
SEXP QuickBLAST::Impl::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
{

  // Dispatch based on the data type of the array
  if (type->id() == arrow::Type::STRUCT)
  {

    auto struct_array = std::static_pointer_cast<arrow::StructArray>(array);
    int num_fields = struct_array->num_fields();

    // Create an Rcpp list to hold the data frames representing each field of the struct
    Rcpp::List struct_list(num_fields);
    Rcpp::CharacterVector names(num_fields);

    for (int i = 0; i < num_fields; i++)
    {

      auto field_array = struct_array->field(i);
      auto field_type = type->field(i)->type();
      auto field_name = type->field(i)->name();
      names[i] = field_name;
      struct_list[i] = Hits2RList_internal(field_array, field_type, field_name);
    }

    struct_list.names() = names;

    return struct_list;
  }
  else if (type->id() == arrow::Type::LIST)
  {

    auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
    auto value_type = type->field(0)->type();

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
      // else
      // {
      //   strings[i] = NA_STRING;
      // }
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
      // else
      // {
      //   ints[i] = NA_INTEGER;
      // }
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

// SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
SEXP QuickBLAST::Impl::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
{
  // Assuming the schema of the RecordBatch is accessible here

  auto rb_schema = rb->schema();

  // Convert each column of the RecordBatch to R objects and store in a list
  Rcpp::List result_list(rb_schema->num_fields());

  for (int i = 0; i < rb_schema->num_fields(); ++i)
  {

    auto array = rb->column(i);
    auto field_type = rb_schema->field(i)->type();
    auto field_name = rb_schema->field(i)->name();
    result_list[i] = Hits2RList_internal(array, field_type, field_name);
  }

  return result_list;
}

// SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
SEXP QuickBLAST::Impl::Hits2RList(const arrow::RecordBatchVector &rb_vector)
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

// std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input)
std::vector<std::pair<std::string, std::string>> QuickBLAST::Impl::BLASTOptionsFromString(const std::string &input)
{
  std::vector<std::pair<std::string, std::string>> keyValuePairs;
  std::istringstream iss(input);
  std::string token;

  while (iss >> token)
  {
    if (token[0] == '-')
    {
      // Extract key-value pair
      std::string key = token.substr(1);
      std::string value;

      if (iss >> value)
      {
        keyValuePairs.emplace_back(key, value);
      }
      else
      {
        // Handle error: Missing value for key
        cerr << "Error: Missing value for key '" << key << "'." << std::endl;
        break;
      }
    }
    else
    {
      // Handle error: Invalid token (not starting with '-')
      cerr << "Error: Invalid token '" << token << "'." << std::endl;
      break;
    }
  }

  return keyValuePairs;
}

// template <typename T1>
// std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return arrow_wrapper->SplitFilesIntoEntries(filename, delim, num_threads, Entry_callback, return_values);
}
// template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return pImpl->StreamFile(filename, delim, num_threads, Entry_callback, return_values);
}

// // template <typename OptionsType>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
// template <typename OptionsType>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
// {
//   this->blast_options = options;
//   this->program = program_name;
//   if constexpr (std::is_same_v<OptionsType, std::string> || std::is_same_v<OptionsType, Rcpp::String>)
//   {
//     // this->blast_options_str = options;
//     return SetQuickBLASTOptions<std::string>(program_name, options);
//   }
//   else if constexpr (std::is_same_v<OptionsType, Rcpp::List>)
//   {
//     // this->blast_options_list = options;
//     return SetQuickBLASTOptions<Rcpp::List>(program_name, options);
//   }
//   else
//   {
//     static_assert(std::is_same_v<OptionsType, OptionsType>, "Unsupported type");
//   }
// }

// // template <>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
// template <>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
// {
//   assert(!program_name.empty());

//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options.empty())
//   {
//     cout << "Using " << program_name << " Defaults..." << std::endl;
//     return opts;
//   }

//   std::vector<std::pair<std::string, std::string>> keyValuePairs = BLASTOptionsFromString(options);

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (const auto &pair : keyValuePairs)
//   {
//     std::string key_str = pair.second;
//     std::size_t key = std::hash<std::string>{}(pair.first);

//     if (key == hashMap["evalue"])
//     {
//       double val = std::stod(key_str);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = std::stod(key_str);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = (key_str == "TRUE" || key_str == "True" || key_str == "true" || key_str == "1");
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = key_str;
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = std::stod(key_str);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = std::stod(key_str);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// }

// // template <>
// // ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
// template <>
// ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
// {
//   assert(!program_name.empty());
//   Rcpp::List options_(options);
//   for (int i = 0; i < options_.size(); ++i)
//   {
//     // Check if the element is empty
//     if (Rf_isNull(options_[i]))
//     {
//       // Remove the empty element
//       options_.erase(i);
//       --i; // Decrement the index since the list size has changed
//     }
//   }
//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options_.size() == 0 || options_.isNULL())
//   {
//     // Rcpp::Rcout << "Using " << program_name << " Defaults..." << std::endl;
//     Rprintf("Using %s Defaults...\n", program_name);
//     return opts;
//   }

//   Rcpp::CharacterVector keys = options_.names();

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (int i = 0; i < keys.size(); ++i)
//   {
//     std::string key_str = Rcpp::as<std::string>(keys[i]);
//     std::size_t key = std::hash<std::string>{}(key_str);

//     if (key == hashMap["evalue"])
//     {
//       double val = Rcpp::as<double>(options_["evalue"]);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = Rcpp::as<double>(options_["pident"]);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = Rcpp::as<bool>(options_["gapped_mode"]);
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = Rcpp::as<std::string>(options_["filter_string"]);
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = Rcpp::as<int>(options_["effective_search_space"]);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = Rcpp::as<int>(options_["cutoff_score"]);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = Rcpp::as<double>(options_["gap_trigger"]);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff"]);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff_final"]);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = Rcpp::as<int>(options_["hit_list_size"]);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = Rcpp::as<double>(options_["low_score_percentage"]);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_subject"]);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_sequence"]);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = Rcpp::as<double>(options_["qcovhsp_perc"]);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = Rcpp::as<int>(options_["window_size"]);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// }

// std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
{
  // assert(num_threads > 0);
  /*   if (!arrow_wrapper || arrow_wrapper.get() == nullptr)
    {
      arrow_wrapper = std::make_shared<ArrowWrapper>();
    } */
  /* if (this->opts.Empty() || this->opts.IsNull())
  {
    // if (this->blast_options_list.size() > 0)
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_list);
    // }
    // else if (!this->blast_options_str.empty())
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_str);
    // }
    // else
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, "");
    // }
    this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, this->blast_options));
  } */

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  unsigned int n_threads = num_threads > omp_get_num_threads() ? omp_get_num_threads() : num_threads;
#else
  unsigned int n_threads = 1;
#endif

  n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
  
  SetThreadCount(n_threads);

  unsigned int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);

  unsigned int s_seq_count = arrow_wrapper->CountCharacter(subjectFile, '>', n_threads);

  if((s_seq_count > 10000 || q_seq_count > 10000) && outFile.empty())
    Rcpp::Rcerr << "Warning: Queries/Subjects > 10000, large inputs can crash during return" << std::endl << std::flush; 
  
  arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);
  if (!outfile_sts.ok())
  {
    /* Rcpp::Rcerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
                   << outfile_sts.message() << std::endl; */
    // cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
    //      << outfile_sts.message() << std::endl;
    // REprintf("ERROR : Could not create output file stream : %s \n %s \n", outfile_sts.detail()->ToString().c_str(), outfile_sts.message().c_str());
    throw std::runtime_error(std::string("ERROR : Could not create output file stream ") + outfile_sts.detail()->ToString() + std::string("\n") + outfile_sts.message());
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  
  const unsigned int totalIterations = q_seq_count * s_seq_count;
  if (blast_sequence_limit > 0)
  {
    blast_sequence_limit = blast_sequence_limit > totalIterations ? totalIterations : blast_sequence_limit;
    blast_sequence_limit = blast_sequence_limit > s_seq_count ? s_seq_count - 1 : blast_sequence_limit; // - 1;
  }
  else
  {
    blast_sequence_limit = s_seq_count - 1;
  }
  assert(totalIterations > 0);
  if(totalIterations <= 0)
    Rcpp::stop("Improperly formatted FASTA file.");

  //// int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
  ////  batch_size = 32 * n_threads; // batch_size > 0 ? batch_size : 1024;
  arrow_wrapper->SetBatchSize(batch_size);

  Progress progress_bar(totalIterations, show_progress);

  std::shared_ptr<arrow::RecordBatchVector> final_ret = StreamFile(
      queryFile, ">", n_threads, [this, &progress_bar, n_threads, subjectFile, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_q)
      {
        // Rcpp::Rcout << data_q->rec_no << "-Q->" << data_q->header << std::endl << std::flush; //DEBUG
        // Rcpp::Rcout << data_q->seq << std::endl << std::flush; //DEBUG
        // std::ostringstream osq; //DEBUG
        // for (size_t i = 0; i < 32; ++i) {
        //   osq << std::hex << std::setw(2) << std::setfill('0')
        //       << (static_cast<int>(static_cast<unsigned char>(data_q->seq[i]))) << (i+1==32 ? "" : " ");
        // }
        // Rcpp::Rcout << "Hex bytes: " << osq.str() << std::endl << std::flush; //DEBUG
        
        bool had_bad_q = false;
        switch (seq_type)
        {
        case ESeqType::eProtein:
        {
          data_q->seq = sanitize_protein(data_q->seq, had_bad_q);
          break;
        }
        case ESeqType::eNucleotide:
        {
          data_q->seq = sanitize_nucleotide(data_q->seq, had_bad_q);
          break;
        }
        }
        if (had_bad_q) {
          Rcpp::Rcout << "Warning: removed/normalized invalid bytes from Query sequence " << data_q->header << std::endl;
        }
        if(data_q->seq.empty()){
          return std::make_shared<arrow::RecordBatchVector>();
        }
            CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
            // const std::shared_ptr<SSeqLoc> query_loc(std::move(CreateSSeqLocFromType<FastaSequenceData>(*data_q, scope)));
            // const std::shared_ptr<SSeqLoc> query_loc(std::move(this->CreateSSeqLocFromType(*data_q, scope)));
            auto [ query_loc, query_seq_entry ] = CreateSSeqLocFromType(*data_q, scope);
            if (!query_loc) {
              Rcpp::stop("BLAST_files: CreateSSeqLocFromType(query) returned NULL.");
            }
            _ASSERT(query_loc->seqloc.NotEmpty());
            
            assert(!Progress::check_abort());
            Rcpp::checkUserInterrupt();
            // progress_bar.increment();
            
            std::shared_ptr<TSeqLocVector> subjects_loc_vec(new TSeqLocVector());
            std::shared_ptr<vector<CSeq_entry_Handle>> subjects_seqent_vec =  std::make_shared<vector<CSeq_entry_Handle>>();

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
            omp_lock_t query_locLock;
            omp_lock_t subjects_loc_vecLock;
            omp_lock_t subjects_seqent_vecLock;
            omp_init_lock(&query_locLock);
            omp_init_lock(&subjects_loc_vecLock);
            omp_init_lock(&subjects_seqent_vecLock);
#endif

            std::shared_ptr<arrow::RecordBatchVector> ret_results = StreamFile(
                subjectFile, ">", n_threads, [this, query_loc,
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                &query_locLock, &subjects_loc_vecLock, &subjects_seqent_vecLock,
#endif
                & scope, &subjects_loc_vec, &subjects_seqent_vec, &progress_bar, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_s)
                {
                    // const std::unique_ptr<SSeqLoc> subject_loc(CreateSSeqLocFromType<FastaSequenceData>(*data_s, scope));
                  bool had_bad_s = false;
                  switch (seq_type)
                  {
                  case ESeqType::eProtein:
                  {
                    data_s->seq = sanitize_protein(data_s->seq, had_bad_s);
                    break;
                  }
                  case ESeqType::eNucleotide:
                  {
                    data_s->seq = sanitize_nucleotide(data_s->seq, had_bad_s);
                    break;
                  }
                  }
                  if (had_bad_s) {
                    Rcpp::Rcout << "Warning: removed/normalized invalid bytes from Query sequence " << data_s->header << std::endl;
                  }
                    if(data_s->seq.empty()){
                      return std::make_shared<arrow::RecordBatchVector>();
                    }
                    
                    // const std::unique_ptr<SSeqLoc> subject_loc(std::move(this->CreateSSeqLocFromType(*data_s, scope))); 
                    auto [ subject_loc, subject_seq_entry] = CreateSSeqLocFromType(*data_s, scope);
                    if (!subject_loc) {
                      Rcpp::stop("BLAST_files: CreateSSeqLocFromType(subject) returned NULL.");
                    }
                    _ASSERT(subject_loc->seqloc.NotEmpty());
     
     
                   //DEBUG
                   // Rcpp::Rcout << "Query id: " << query_loc->seqloc->GetId()->GetSeqIdString(true) << std::endl;
                   // Rcpp::Rcout << "Subject id: " << subject_loc->seqloc->GetId()->GetSeqIdString(true) << std::endl;
                  //DEBUG
                    // Rcpp::Rcout << "strcmp(): " << strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) << std::endl << std::flush; //DEBUG
                    // 
                    // const CSeq_id* id1 = query_loc->seqloc->GetId();
                    // const CSeq_id* id2 = subject_loc->seqloc->GetId();
                    // if (!id1 || !id2) { /* handle multi-part locs or error */ }
                    // 
                    // CSeq_id_Handle h1 = CSeq_id_Handle::GetHandle(*id1);
                    // CSeq_id_Handle h2 = CSeq_id_Handle::GetHandle(*id2);
                    
                    if (strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) != 0)
                    {
                        assert(!Progress::check_abort());
                        Rcpp::checkUserInterrupt();
                        // progress_bar.increment();
                        
                        CBl2Seq* blaster;

                        try
                        {

                            switch (blast_sequence_limit)
                            {
                            case 0:
                            {
                                blaster = new CBl2Seq(*query_loc, *subject_loc, this->GetQuickBLASTOptions());
                                // arrow::RecordBatchVector tmp_rbv = { ExtractHits<SSeqLoc>(blaster->Run(), *query_loc, *subject_loc, *scope) };
                                TSeqAlignVector alignments = blaster->Run();
                                AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);
                                
                                arrow::RecordBatchVector tmp_rbv = { ExtractHits(alignments, *query_loc, *subject_loc, *scope, progress_bar, return_values) }; 
                                
                                scope->RemoveTopLevelSeqEntry(subject_seq_entry);
                                // subject_seq_entry.Reset();
                                subject_loc.reset();
                                
                                if (return_values)
                                {
                                    return std::make_shared<arrow::RecordBatchVector>(tmp_rbv);
                                }
                                else
                                {
                                    tmp_rbv.clear();
                                    return std::make_shared<arrow::RecordBatchVector>();
                                }
                            }
                            break;
                            default:
                            {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                omp_set_lock(&subjects_loc_vecLock);
                                omp_set_lock(&subjects_seqent_vecLock);
#endif
                                subjects_loc_vec->emplace_back(*subject_loc);
                                subjects_seqent_vec->emplace_back(subject_seq_entry);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                omp_unset_lock(&subjects_loc_vecLock);
                                omp_unset_lock(&subjects_seqent_vecLock);
#endif

                                if (subjects_loc_vec->size() >= blast_sequence_limit)
                                {
                                    TSeqLocVector subjects_buffer_vec;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                    omp_set_lock(&subjects_loc_vecLock);

#endif
                                    subjects_buffer_vec.reserve(subjects_loc_vec->size());
                                    subjects_buffer_vec.swap(*subjects_loc_vec);
                                    subjects_loc_vec->clear();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                    omp_unset_lock(&subjects_loc_vecLock);
#endif
                                    blaster = new CBl2Seq(*query_loc, subjects_buffer_vec, this->GetQuickBLASTOptions(), true);

                                    AddHitCount(subjects_buffer_vec.size());
                                    TSeqAlignVector alignments = blaster->Run();
                                    AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);
                                    
                                    std::shared_ptr<arrow::RecordBatchVector> tmp_rbv = ExtractHits(alignments, *query_loc, subjects_buffer_vec, *scope, progress_bar, return_values); 
                                    //return ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                    omp_set_lock(&subjects_seqent_vecLock);
#endif
                                    for(auto s_ent: *subjects_seqent_vec){
                                      scope->RemoveTopLevelSeqEntry(s_ent);
                                      // s_ent.Reset();
                                    }
                                    subjects_seqent_vec->clear();
                                    subjects_seqent_vec->shrink_to_fit();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                    omp_unset_lock(&subjects_seqent_vecLock);
#endif
                                    // for(auto s_loc: subjects_buffer_vec){
                                    //   s_loc.reset();
                                    // }
                                    
                                    if (return_values)
                                    {
                                        return tmp_rbv;
                                    }
                                    else
                                    {
                                        tmp_rbv->clear();
                                        return std::make_shared<arrow::RecordBatchVector>();
                                    }
                                }
                                
                            }
                            break;
                            }
                        }
                        catch (const CException& e)
                        {
                          // Rcpp::Rcout << data_s->rec_no << "-S->" << data_s->header << std::endl << std::flush; //DEBUG
                          // Rcpp::Rcout << data_s->seq << std::endl << std::flush; //DEBUG
                          // std::ostringstream oss; //DEBUG
                          // for (size_t i = 0; i < 32; ++i) {
                          //   oss << std::hex << std::setw(2) << std::setfill('0')
                          //       << (static_cast<int>(static_cast<unsigned char>(data_s->seq[i]))) << (i+1==32 ? "" : " ");
                          // }
                          // Rcpp::Rcout << "Hex bytes: " << oss.str() << std::endl << std::flush; //DEBUG
                          
                          // Handle exception ...
                          Rcpp::Rcerr << e.GetFunction() << std::endl;
                          Rcpp::Rcerr << e.GetErrCodeString() << std::endl;
                            Rcpp::Rcerr << e.GetErrCode() << std::endl;
                            Rcpp::Rcerr << e.GetModule() << std::endl;
                            Rcpp::Rcerr << e.GetPredecessor() << std::endl;
                            Rcpp::Rcerr << e.GetFile() << std::endl;
                            Rcpp::Rcerr << e.GetLine() << std::endl;
                            Rcpp::Rcerr << e.GetMsg() << std::endl;
                            Rcpp::Rcerr << e.GetStackTrace() << std::endl;
                            Rcpp::Rcerr << e.GetStackTraceLevel() << std::endl;
                            Rcpp::Rcerr << e.GetClass() << std::endl
                                << std::flush;
                            Rcpp::stop("BLAST_files(): NCBI Exception: Stopping run");
                            // throw;
                        }
                        catch (const std::exception& e)
                        {
                            Rcpp::stop(std::string("BLAST_files(): C++ Exception: ") + e.what());
                        }
                        catch(const Rcpp::exception &e){
                          Rcpp::stop(std::string("BLAST_files() - Rcpp Exception : ") + e.what());
                        }catch(...){
                          Rcpp::stop("BLAST_files() - Unknown Exception");
                        }
                        
                        // Rcpp::Rcout << "Processing Subject: " << data_s->rec_no << std::endl;
                        
                    }

                    return std::make_shared<arrow::RecordBatchVector>(); // EMPTY ERROR Return
                },
                return_values);

            if (subjects_loc_vec->size() > 0)
            {
                CBl2Seq blaster(*query_loc, *subjects_loc_vec, this->GetQuickBLASTOptions(), true);
                AddHitCount(subjects_loc_vec->size());
                TSeqAlignVector alignments = blaster.Run();
                AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);
                
                std::shared_ptr<arrow::RecordBatchVector> ret_vec = ExtractHits(alignments, *query_loc, *subjects_loc_vec, *scope, progress_bar, return_values);

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                omp_set_lock(&subjects_seqent_vecLock);
#endif
                for(auto s_ent: *subjects_seqent_vec){
                  scope->RemoveTopLevelSeqEntry(s_ent);
                  // s_ent.Reset();
                }
                subjects_seqent_vec->clear();
                subjects_seqent_vec->shrink_to_fit();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                omp_unset_lock(&subjects_seqent_vecLock);
#endif
                
                // for(auto s_loc: subjects_loc_vec.get()){
                //   s_loc.reset();
                // }
                
                if (return_values)
                {
                    ret_results->insert(ret_results->end(), ret_vec->begin(), ret_vec->end());
                }
                else
                {
                    ret_vec->clear();
                }
            }

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
            omp_destroy_lock(&query_locLock);
            omp_destroy_lock(&subjects_loc_vecLock);
#endif
            scope->RemoveTopLevelSeqEntry(query_seq_entry);
            // query_seq_entry.Reset();
            query_loc.reset();
            
            if (return_values) {
              return ret_results;
            }
            else
            {
              ret_results->clear();
              return std::make_shared<arrow::RecordBatchVector>();
            } },
      return_values);

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp barrier
#endif


  Rcpp::Rcout << "Final Batch Size: " << arrow_wrapper->GetBatchSize() << std::endl << std::flush; 

  arrow_wrapper->FinishOutputStream();
  
  if (return_values)
  {
    return final_ret;
  }
  else
  {
    final_ret->clear();
    final_ret->shrink_to_fit();
    return std::make_shared<arrow::RecordBatchVector>();
  }
}

// std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject)
std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::BLAST_seqs(const std::string &query, const std::string &subject)
{
  try {
    // 0) quick user interrupt check if this runs under R
    Rcpp::checkUserInterrupt();
    // 1) ensure arrow_wrapper exists
    if (!this->arrow_wrapper) {
      Rcpp::stop("BLAST_seqs: arrow_wrapper is null (not initialized).");
    }
    // 2) convert inputs via arrow wrapper and validate
    auto q_type = this->arrow_wrapper->CastToType(query);
    auto s_type = this->arrow_wrapper->CastToType(subject);
    if (q_type.header.empty() || q_type.seq.empty()) {
      Rcpp::stop("BLAST_seqs: query header/sequence is empty.");
    }
    if (s_type.header.empty() || s_type.seq.empty()) {
      Rcpp::stop("BLAST_seqs: subject header/sequence is empty.");
    }

    // 3) get CObjectManager instance and create scope
    CRef<CObjectManager> objMgr(CObjectManager::GetInstance());
    if (!objMgr) {
      Rcpp::stop("BLAST_seqs: CObjectManager::GetInstance() returned NULL.");
    }

    CRef<ncbi::CScope> scope(new ncbi::CScope(*objMgr));
    // optional: add defaults or other initialization your app normally does:
    scope->AddDefaults();
    // 4) build SSeqLoc objects and check for null
    
    // const std::unique_ptr<SSeqLoc> query_seqloc(std::move(this->CreateSSeqLocFromType(q_type, scope)));
    auto [ query_seqloc, query_seq_entry ] = CreateSSeqLocFromType(q_type, scope);
    if (!query_seqloc) {
      Rcpp::stop("BLAST_seqs: CreateSSeqLocFromType(query) returned NULL.");
    }
    _ASSERT(query_seqloc->seqloc.NotEmpty());
    // const std::unique_ptr<SSeqLoc> subject_seqloc(std::move(this->CreateSSeqLocFromType(s_type, scope)));
    auto [ subject_seqloc, subject_seq_entry ] = CreateSSeqLocFromType(s_type, scope);
    if (!subject_seqloc) {
      Rcpp::stop("BLAST_seqs: CreateSSeqLocFromType(subject) returned NULL.");
    }
    _ASSERT(subject_seqloc->seqloc.NotEmpty());
    // 5) create the blaster and run
    CBl2Seq blaster(*query_seqloc, *subject_seqloc, this->GetQuickBLASTOptions());
    // // run and extract hits
    // CRef<ncbi::CSeq_align_set> result = blaster.Run();
    // if ( ! result ) {
    //   Rcpp::stop("BLAST_seqs: blaster.Run() returned NULL.");
    // }
    TSeqAlignVector alignments = blaster.Run();
    // CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);
    
    Progress progress_bar(1, true);
    
    return this->ExtractHits(alignments, *query_seqloc, *subject_seqloc, *scope, progress_bar, true); 
  }
  catch (const ncbi::CException &e) {
    // NCBI toolkit exceptions
    std::string msg = "NCBI Toolkit exception in BLAST_seqs: ";
    try { msg += e.GetMsg(); } catch(...) { msg += "(failed to read message)"; }
    Rcpp::stop(msg);
  }
  catch (const std::exception &e) {
    Rcpp::stop(std::string("BLAST_seqs() - C++ exception : ") + e.what());
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("BLAST_seqs() - Rcpp Exception : ") + e.what());
  }
  catch (...) {
    Rcpp::stop("BLAST_seqs(): Unknown exception");
  }
  
  // // Rcpp::checkUserInterrupt();
  // 
  // CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
  // 
  // // std::unique_ptr<SSeqLoc>
  // //     query_seqloc(CreateSSeqLocFromType<std::string>(query, scope));
  // // std::unique_ptr<SSeqLoc> subject_seqloc(CreateSSeqLocFromType<std::string>(subject, scope));
  // 
  // std::unique_ptr<SSeqLoc>
  //     query_seqloc(CreateSSeqLocFromType(arrow_wrapper->CastToType(query), scope));
  // std::unique_ptr<SSeqLoc> subject_seqloc(CreateSSeqLocFromType(arrow_wrapper->CastToType(subject), scope));
  // 
  // CBl2Seq blaster(*query_seqloc, *subject_seqloc, GetQuickBLASTOptions());
  // 
  // // return ExtractHits<SSeqLoc>(blaster.Run(), *query_seqloc, *subject_seqloc, *scope);
  // return ExtractHits(blaster.Run(), *query_seqloc, *subject_seqloc, *scope);
}

//' @name BLAST C++ Call
//' @title BLAST C++ Call
//'
//' @description BLAST 2 Files/Seqs. This is for the QuickBLAST C++ object exposed in R
//'
//' @param query (string) Query FASTA File/Seq
//' @param subject (string) Subject FASTA File/Seq
//' @param outputFile (string) Output Filename (Arrow Feather/IPC Format)  - Not used for Sequence BLAST
//' @param input_type - (QuickBLAST::EInputType) 0 - eFile, 1 - eSequenceString
//' @param blast_sequence_limit (int) Batch Size to BLAST at a time  - Not used for Sequence BLAST
//' @param show_progress (bool) TRUE - Show progress, Set FALSE for multiple instances - Not used for Sequence BLAST
//' @return Nested List of BLAST Hits
// auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
auto QuickBLAST::Impl::BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
{

  assert(std::filesystem::exists(query));
  assert(std::filesystem::exists(subject));

  switch (input_type)
  {
  case QuickBLAST::EInputType::eFile:
  {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    unsigned int n_threads = omp_get_num_threads();
#else
    unsigned int n_threads = 1;
#endif
    int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
    batch_size = 32 * n_threads;
    if(blast_sequence_limit > 0)
      batch_size = batch_size > 0 ? batch_size : blast_sequence_limit;
    else
      batch_size = batch_size > 0 ? batch_size : 0;
    return BLAST_files(query, subject, outputFile, outFormat, blast_sequence_limit, n_threads, show_progress, true, batch_size);
    // return Hits2RList(*ret_val);
  }
  break;
  case QuickBLAST::EInputType::eSequenceString:
  {
    arrow::RecordBatchVector ret_val;
    ret_val.emplace_back(BLAST_seqs(query, subject));
    return std::make_shared<arrow::RecordBatchVector>(ret_val);
    // return Hits2RList(ret_val);
  }
  break;
  default:
  {
    // Rcpp::Rcerr << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
    // cout << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
    REprintf("input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !");
    // return false; //Rcpp::wrap(false);
    return std::make_shared<arrow::RecordBatchVector>();
  }
  break;
  }
}

std::string QuickBLAST::GetQuickBLASTOptionString(){
  return pImpl->GetQuickBLASTOptionString(); 
}

std::string QuickBLAST::Impl::GetQuickBLASTOptionString(){
  return blast_options;
}

ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality = CBlastOptions::eLocal)
{
  return pImpl->SetQuickBLASTOptions(program_name, options, locality);
}

unsigned int QuickBLAST::GetObjectID()
{
  return pImpl->GetObjectID();
}

void QuickBLAST::SetObjectID(unsigned int id)
{
  pImpl->SetObjectID(id);
}

std::shared_ptr<arrow::Schema> QuickBLAST::GetSchema() { return pImpl->GetSchema(); };

std::string QuickBLAST::GetProgram(){ return pImpl->GetProgram(); };

std::string QuickBLAST::Impl::GetProgram(){
  return program;
}

void QuickBLAST::SetThreadCount(unsigned int num_threads)
{
  pImpl->SetThreadCount(num_threads);
}
unsigned int QuickBLAST::GetThreadCount()
{
  return pImpl->GetThreadCount();
}
int QuickBLAST::GetHitCount()
{
  return pImpl->GetHitCount();
}
void QuickBLAST::AddHitCount(int val)
{
  pImpl->AddHitCount(val);
}
ncbi::blast::CBlastOptionsHandle &QuickBLAST::GetQuickBLASTOptions()
{
  return pImpl->GetQuickBLASTOptions();
}
void QuickBLAST::ResetHitCount() { pImpl->ResetHitCount(); }

// template <>
// SSeqLoc *QuickBLAST::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope){
//   return pImpl->CreateSSeqLocFromType(fasta_data, parent_scope);
// }

// int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand){
//   return pImpl->GetFrame(start, length, strand);
// }

// template <>
// std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, const CScope &scope){
//   return pImpl->ExtractHits(alignments, qloc, sloc, scope);
// }

// std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data){
//   return pImpl->ExtractFASTA(fasta_data);
// }

// std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc){
//   return pImpl->GetSSeqLocSequence(seq_loc);
// }

SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
{
  return pImpl->Hits2RList(rb);
}

SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
{
  return pImpl->Hits2RList(rb_vector);
}

// std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input){
//   return pImpl->BLASTOptionsFromString(input);
// }

// template <typename OptionsType>
// ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options){
//   return pImpl->SetQuickBLASTOptions(program_name, options);
// }

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const TSeqLocVector &sloc, CScope &scope, Progress &progress_bar, const bool &return_values) 
{
  try{
    Rcpp::checkUserInterrupt();
    assert(!alignments.empty());
    
    std::shared_ptr<arrow::RecordBatchVector> ret_rbv;
    
    if (!qloc.seqloc) {
      Rprintf("ERROR: ExtractHits: qloc.seqloc is NULL\n");
      //ret_rbv->emplace_back(arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie());
      ret_rbv->emplace_back(empty_rb);
      return ret_rbv;
    }
    // if (!sloc.seqloc) {
    //   Rprintf("ERROR: ExtractHits: sloc.seqloc is NULL\n");
    //   ret_rbv->emplace_back(arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie());
    //   return ret_rbv;
    // }
    
    // if constexpr (std::is_same_v<T, TSeqLocVector>)
    // {
    std::shared_ptr<arrow::RecordBatchVector> recBth_vec = std::make_shared<arrow::RecordBatchVector>();
    
    for (const auto &s_it : sloc)
    {
      Rcpp::checkUserInterrupt();
  
      std::shared_ptr<arrow::RecordBatch> rb = ExtractHits(alignments, qloc, s_it, scope, progress_bar, return_values);
  
      if(return_values)
        if (rb)
        {
          recBth_vec->emplace_back(std::move(rb));
        }
    }
  
    // if(return_values){
      return recBth_vec;
    // }
    
    //Already adding RBs to Batch in the overloaded ExtractHits()
    // else{
    //   const auto &wrt_sts = arrow_wrapper->AddRBV2Batch(*recBth_vec);
    //   // recBth_vec.clear();
    //   // recBth_vec.shrink_to_fit();
    //   if (!wrt_sts.ok())
    //   {
    //     throw std::runtime_error(std::string("Warn : Invalid Alignment RBV (Returning Empty) : ") + wrt_sts.detail()->ToString() + std::string(" : ") + wrt_sts.message());
    //   }
    //   return std::make_shared<arrow::RecordBatchVector>();
    // }
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("ExtractHits(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("ExtractHits(): C++ Runtime Error : ") + e.what());
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("ExtractHits() - Rcpp Exception : ") + e.what());
  }
  catch(...){
    Rcpp::stop("ExtractHits(): Unknown Exception");
  }
  // }
  // // For SSeqLoc
  // else if constexpr (std::is_same_v<T, SSeqLoc>)
  // {
  //     std::shared_ptr<arrow::RecordBatch> rb = ExtractHits(alignments, qloc, sloc, scope);

  //     const auto &wrt_sts = arrow_wrapper->AddRB2Batch(rb);
  //     if (wrt_sts.ok())
  //     {
  //         return wrt_sts.ValueOrDie();
  //     }
  //     else
  //     {
  //         cerr << "Warn : Invalid Alignment RB (Returning Empty) : " << wrt_sts.status().detail() << std::endl
  //              << wrt_sts.status().message() << std::endl;
  //     }

  //     return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
  // }
  // else
  // {
  //     static_assert(std::is_same_v<T, T>, "Unsupported type, only ncbi::blast::TSeqLocVector & ncbi::blast::SSeqLoc are supported");
  // }
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, CScope &scope, Progress &progress_bar, const bool &return_values) 
{
  try{
    Rcpp::checkUserInterrupt();
    // assert(!alignments.empty());
    if (alignments.empty()) {
      // return an empty but typed record batch
      progress_bar.increment();
      return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
    }
  
    if (!qloc.seqloc) {
      Rprintf("ERROR: ExtractHits: qloc.seqloc is NULL\n");
      progress_bar.increment();
      return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
    }
    if (!sloc.seqloc) {
      Rprintf("ERROR: ExtractHits: sloc.seqloc is NULL\n");
      progress_bar.increment();
      return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
    }
    
    // CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
    
    std::string strand;
    arrow::StringBuilder strand_builder;
    
    auto query_strand = qloc.seqloc->GetStrand();
    auto subject_strand = sloc.seqloc->GetStrand();
  
    switch (query_strand)
    {
    case eNa_strand_minus:
      strand = strand + "-";
      break;
    case eNa_strand_plus:
      strand = strand + "+";
      break;
    default:
      strand = strand + "*";
      break;
    }
  
    switch (subject_strand)
    {
    case eNa_strand_minus:
      strand = strand + "/-";
      break;
    case eNa_strand_plus:
      strand = strand + "/+";
      break;
    default:
      strand = strand + "/*";
      break;
    }
  
    std::string qseq = "", sseq = "", frame = "*/*", qseq_id, sseq_id;
    arrow::StringBuilder qseqid_builder, sseqid_builder; // qseq_title_builder, sseq_title_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder;
    arrow::LargeStringBuilder qhsp_builder, shsp_builder;
    arrow::Int64Builder qlen_builder, slen_builder, num_alignments_builder;
  
    qseq_id = qloc.seqloc->GetId()->GetSeqIdString(true);
    sseq_id = sloc.seqloc->GetId()->GetSeqIdString(true);
    
    // Rcpp::Rcout << "here2.0.0" << std::endl << std::flush;
    // CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(*qloc.seqloc->GetId());
    // Rcpp::Rcout << "here2.0.1" << std::endl << std::flush;
    // CBioseq_Handle q_bh = scope->GetBioseqHandle(q_idh);
    // Rcpp::Rcout << "here2.0.2" << std::endl << std::flush;
    // // const auto q_b = q_bh.GetCompleteObject();
    // std::string qseq_title = qseq_id;
    // CSeqdesc_CI qdesc(q_bh, CSeqdesc::e_Title);
    // if (qdesc) {
    //   const CSeqdesc& qd = *qdesc;
    //   if (qd.IsTitle()) {
    //     // const auto &titles = d.GetTitle().Get();
    //     // if (!titles.empty()) return titles.front();
    //     if (qd.IsTitle() && !qd.GetTitle().empty()) 
    //       qseq_title = qd.GetTitle();
    //   }
    // }
    // Rcpp::Rcout << "here2.0.3" << std::endl << std::flush;
    // // const auto qdesc = q_b->GetDescr().Get();
    // // Rcpp::Rcout << "here2.0.4" << std::endl << std::flush;
    // // std::string qseq_title = qseq_id;
    // // for (auto &d : qdesc) {
    // //   if (d->IsTitle() && !d->GetTitle().empty()) 
    // //     qseq_title = d->GetTitle();
    // // }
    // Rcpp::Rcout << "here2.0.5" << std::endl << qseq_title << std::endl << std::flush;
    //   
    // CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(*sloc.seqloc->GetId());
    // CBioseq_Handle s_bh = scope->GetBioseqHandle(s_idh);
    // // const auto s_b = s_bh.GetCompleteObject();
    // // const auto sdesc = s_b->GetDescr().Get();
    // std::string sseq_title = sseq_id;
    // CSeqdesc_CI sdesc(s_bh, CSeqdesc::e_Title);
    // if (sdesc) {
    //   const CSeqdesc& sd = *sdesc;
    //   if (sd.IsTitle()) {
    //     // const auto &titles = d.GetTitle().Get();
    //     // if (!titles.empty()) return titles.front();
    //     if (sd.IsTitle() && !sd.GetTitle().empty()) 
    //       sseq_title = sd.GetTitle();
    //   }
    // }
    // // for (auto &d : sdesc) {
    // //   if (d->IsTitle() && !d->GetTitle().empty()) 
    // //     sseq_title = d->GetTitle();
    // // }
    // Rcpp::Rcout << "here2.0.6" << std::endl << sseq_title << std::endl << std::flush;
    
    // qseq = "";
    // sseq = "";
    // switch (save_sequences)
    // {
    // case true:
    //   qseq = GetSSeqLocSequence(qloc);
    //   sseq = GetSSeqLocSequence(sloc);
    //   break;
    // }
  
    arrow::Int64Builder hsp_offset_builder;
  
    int num_rows = 0;
  
    arrow::Int64Builder length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    arrow::StringBuilder frames_builder;
    
    for (const auto &seq_align_set : alignments)
    {
  
      if (seq_align_set->IsEmpty())
      {
        break;
      }
      assert(seq_align_set->IsSet());
      assert(seq_align_set->CanGet());
      const auto &seq_aligns = seq_align_set->Get();
      assert(!seq_aligns.empty());
  
      if (seq_aligns.size() > 0) // FILL UP THE ARRAYS
      {
  
        try
        {
  
          for (const auto &it : seq_aligns)
          {
            if (!it)
            {
              continue;
            }
            assert(!it.IsNull());
            if (!it.NotEmpty())
            {
              continue;
            }
            
            // const auto& seq_titles = GetTitlesFromSeqAlign(it, scope);
            // std::string qseq_title = seq_titles.first();
            // std::string sseq_title = seq_titles.second();
            
            assert(it->CanGetScore());
            // int score, n_splices, num_ident, aln_len, gaps, mismatches, positive, qstart, qend, sstart, send, negative_count;
            // double bits, evalue, blast_score, pident, aln_len01, pident_gap, qcovhsp, sum_evalue, product_coverage, overall_identity, high_quality_percent_coverage, exon_identity, consensus_splices, comp_adj_method, matches;
            // std::string frames;
            // 
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices);
            // 
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices);
            // it->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method);
            // 
            // aln_len01 = it->AlignLengthRatio();
            // 
            // qstart = it->GetSeqStart(0);
            // qend = it->GetSeqStop(0);
            // sstart = it->GetSeqStart(1);
            // send = it->GetSeqStop(1);
        
            std::string q_full = "", s_full = "";
            std::string q_hsp = "", s_hsp = "", q_aligned = "", s_aligned = "";
            
            if (it->GetSegs().IsDenseg()) {
              const CDense_seg& dseg = it->GetSegs().GetDenseg();
              
              // // Get sequence ids (rows)
              // if (dseg.CanGetIds()) {
              //   const auto &ids = dseg.GetIds();
              //   // print/inspect id strings:
              //   for (size_t r = 0; r < ids.size(); ++r) {
              //     if (ids[r]) {
              //       NcbiCout << "Row " << r << " id: " << ids[r]->GetSeqIdString(true) << NcbiEndl;
              //     }
              //   }
              // }
              
              // Full sequences for the two first rows (query, subject)
              if (dseg.CanGetIds()) {
                // try to fetch full sequences for rows 0 and 1
                if (dseg.GetIds().size() > 0) {
                  GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[0]), scope, q_full);
                }
                if (dseg.GetIds().size() > 1) {
                  GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[1]), scope, s_full);
                }
              }
              
              if(save_hsp_sequences){
                // HSP sequences
                bool ok = GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned);
              }
              // NcbiCout << "Full query length: " << q_full.size() << " HSP ungapped length: " << q_hsp.size() << NcbiEndl;
              // NcbiCout << "Full subject length: " << s_full.size() << " HSP ungapped length: " << s_hsp.size() << NcbiEndl;
              // NcbiCout << "Aligned strings length: " << q_aligned.size() << " / " << s_aligned.size() << NcbiEndl;
              // NcbiCout << "Query FASTA: " << q_full.substr(0, 200) << NcbiEndl;   // print only prefix
              // NcbiCout << "Subject FASTA: " << s_full.substr(0, 200) << NcbiEndl;
              // NcbiCout << "Query HSP: " << q_hsp.substr(0, 200) << NcbiEndl;   // print only prefix
              // NcbiCout << "Subject HSP: " << s_hsp.substr(0, 200) << NcbiEndl;
            }
           
           switch (save_sequences)
           {
           case true:
             qseq = q_full;
             sseq = s_full;
             break;
           }
            
            // // handle Std-seg (a sequence of local 'loc' entries)
            // else if (it->GetSegs().IsStd()) {
            //   // const CStd_seg &stdseg = it->GetSegs().GetStd();
            //   // stdseg has a list of segments; each segment has a list of locs for each row
            //   // iterate and extract using the loc's intervals
            //   // For brevity, here's a simple approach that attempts to extract by using GetSeqStart/GetSeqStop
            //   int qstart = it->GetSeqStart(0);
            //   int qstop  = it->GetSeqStop(0);
            //   int sstart = it->GetSeqStart(1);
            //   int sstop  = it->GetSeqStop(1);
            //   // fetch sequences by slicing the bioseq handles (if available)
            //   // (You may prefer to iterate stdseg.Get() entries to get exact block-level offsets)
            //   NcbiCout << "Std-seg: q[" << qstart << "," << qstop << "] s[" << sstart << "," << sstop << "]" << NcbiEndl;
            //   // you can reuse GetFullSequenceString + substringing with CSeqVector for exact subrange
            // }
            // else {
            //   // Other seg types: disc, spliced, packed-int, etc.
            //   NcbiCout << "Unhandled seg type; implement specialized extraction if needed" << NcbiEndl;
            // }
            
            double aln_len01 = 0;
            double aln_len = 0;
            double bits = 0.0;
            double blast_score = 0.0;
            double pident = 0.0;
            double pident_gap = 0.0;
            double gaps = 0;
            double evalue = 0.0;
            double num_ident = 0;
            double mismatches = 0;
            double qcovhsp = 0.0;
            double score = 0;
            double positive = 0;
            double n_splices = 0;
            double sum_evalue = 0.0;
            double product_coverage = 0.0;
            double overall_identity = 0.0;
            double negative_count = 0;
            double matches = 0.0;
            double high_quality_percent_coverage = 0.0;
            double exon_identity = 0.0;
            double consensus_splices = 0.0;
            double comp_adj_method = 0.0;
            std::string frames = "*/*";
            
            // For each requested score, call GetNamedScore and check result
            bool ok;
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
            // Rcpp::Rcout << "AlignLength present: " << ok << " value: " << aln_len << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits);
            // Rcpp::Rcout << "BitScore present: " << ok << " value: " << bits << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
            // Rcpp::Rcout << "Blast score present: " << ok << " value: " << blast_score << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
            // Rcpp::Rcout << "IdentityCount present: " << ok << " value: " << num_ident << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident);
            // Rcpp::Rcout << "PercentIdentity_Ungapped present: " << ok << " value: " << pident << std::endl;
            
            if(!ok)
              if (aln_len > 0){
                pident = 100.0 * (double)num_ident / (double)aln_len;
                // Rcpp::Rcout << "computed pident = " << pident << std::endl;
              }
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap);
            // Rcpp::Rcout << "PercentIdentity (gapped) present: " << ok << " value: " << pident_gap << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps);
            // Rcpp::Rcout << "GapCount present: " << ok << " value: " << gaps << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue);
            // Rcpp::Rcout << "EValue present: " << ok << " value: " << evalue << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
            // Rcpp::Rcout << "MismatchCount present: " << ok << " value: " << mismatches << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp);
            // Rcpp::Rcout << "PercentCoverage present: " << ok << " value: " << qcovhsp << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score);
            // Rcpp::Rcout << "Score present: " << ok << " value: " << score << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
            // Rcpp::Rcout << "PositiveCount present: " << ok << " value: " << positive << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices);
            // Rcpp::Rcout << "Splices present: " << ok << " value: " << n_splices << std::endl;
            
            // extended ones
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
            // Rcpp::Rcout << "SumEValue present: " << ok << " value: " << sum_evalue << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
            // Rcpp::Rcout << "ProductCoverage present: " << ok << " value: " << product_coverage << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
            // Rcpp::Rcout << "OverallIdentity present: " << ok << " value: " << overall_identity << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count);
            // Rcpp::Rcout << "NegativeCount present: " << ok << " value: " << negative_count << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches);
            // Rcpp::Rcout << "Matches present: " << ok << " value: " << matches << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage);
            // Rcpp::Rcout << "HQPercentCoverage present: " << ok << " value: " << high_quality_percent_coverage << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity);
            // Rcpp::Rcout << "ExonIdentity present: " << ok << " value: " << exon_identity << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices);
            // Rcpp::Rcout << "ConsensusSplices present: " << ok << " value: " << consensus_splices << std::endl;
            
            ok = it->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method);
            // Rcpp::Rcout << "CompAdjMethod present: " << ok << " value: " << comp_adj_method << std::endl;
            
            // Seq positions
            int qstart = it->GetSeqStart(0);
            int qend   = it->GetSeqStop(0);
            int sstart = it->GetSeqStart(1);
            int send   = it->GetSeqStop(1);
            aln_len01 = it->AlignLengthRatio();
            
            // Rcpp::Rcout << "qstart=" << qstart << " qend=" << qend << " sstart=" << sstart << " send=" << send << std::endl;
            // Rcpp::Rcout << "aln_len (reported): " << aln_len << " nident=" << num_ident << " pident=" << pident << " mismatches=" << mismatches << " gaps=" << gaps << std::endl;
  
            frames = std::to_string(GetFrame(qstart, aln_len, query_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, subject_strand));
            
            static_cast<void>(qhsp_builder.Append(q_hsp));
            static_cast<void>(shsp_builder.Append(s_hsp));
            static_cast<void>(frames_builder.Append(frames));
            static_cast<void>(qstart_builder.Append(qstart));
            static_cast<void>(qend_builder.Append(qend));
            static_cast<void>(sstart_builder.Append(sstart));
            static_cast<void>(send_builder.Append(send));
            static_cast<void>(pident_builder.Append(pident));
            static_cast<void>(evalue_builder.Append(evalue));
            static_cast<void>(length_builder.Append(aln_len));
            static_cast<void>(aln_len01_builder.Append(aln_len01));
            static_cast<void>(bitscore_builder.Append(bits));
            static_cast<void>(score_builder.Append(score));
            static_cast<void>(qcovhsp_builder.Append(qcovhsp));
            static_cast<void>(blast_score_builder.Append(blast_score));
            static_cast<void>(pident_gap_builder.Append(pident_gap));
            static_cast<void>(gaps_builder.Append(gaps));
            static_cast<void>(nident_builder.Append(num_ident));
            static_cast<void>(mismatch_builder.Append(mismatches));
            static_cast<void>(positive_builder.Append(positive));
            static_cast<void>(n_splices_builder.Append(n_splices));
            static_cast<void>(hsp_cnt_builder.Append(num_rows + 1));
            static_cast<void>(sum_evalue_builder.Append(sum_evalue));
            static_cast<void>(product_coverage_builder.Append(product_coverage));
            static_cast<void>(overall_identity_builder.Append(overall_identity));
            static_cast<void>(negative_count_builder.Append(negative_count));
            static_cast<void>(matches_builder.Append(matches));
            static_cast<void>(high_quality_percent_coverage_builder.Append(high_quality_percent_coverage));
            static_cast<void>(exon_identity_builder.Append(exon_identity));
            static_cast<void>(consensus_splices_builder.Append(consensus_splices));
            static_cast<void>(comp_adj_method_builder.Append(comp_adj_method));
  
            /// SEQ INFO
            static_cast<void>(qseqid_builder.Append(qseq_id));
            static_cast<void>(sseqid_builder.Append(sseq_id));
            static_cast<void>(qseq_builder.Append(qseq));
            static_cast<void>(sseq_builder.Append(sseq));
            static_cast<void>(qlen_builder.Append(q_full.length()));
            static_cast<void>(slen_builder.Append(s_full.length()));
            static_cast<void>(num_alignments_builder.Append(seq_aligns.size()));
  
            static_cast<void>(strand_builder.Append(strand));
            static_cast<void>(hsp_offset_builder.Append(1));
  
            // static_cast<void>(qseq_title_builder.Append(qseq_title));
            // static_cast<void>(sseq_title_builder.Append(sseq_title));
            
            num_rows++;
          }
        }
        catch(const std::exception &e){
          Rcpp::stop(std::string("ExtractHits(): C++ Exception : ") + e.what());
        }
        catch(const std::runtime_error &e){
          Rcpp::stop(std::string("ExtractHits(): C++ Runtime Error : ") + e.what());
        }
        catch(const Rcpp::exception &e){
          Rcpp::stop(std::string("ExtractHits() - Rcpp Exception : ") + e.what());
        }
        catch(...){
          Rcpp::stop("ExtractHits(): Unknown Exception");
        }
      }
      else
      {
        return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
      }
    }
    
    if (num_rows == 0)
    {
      return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
    }
  
    std::shared_ptr<arrow::Array> qhsp_array;
    static_cast<void>(qhsp_builder.Finish(&qhsp_array));
    std::shared_ptr<arrow::Array> shsp_array;
    static_cast<void>(shsp_builder.Finish(&shsp_array));
    std::shared_ptr<arrow::Array> frames_array;
    static_cast<void>(frames_builder.Finish(&frames_array));
    std::shared_ptr<arrow::Array> pident_array;
    static_cast<void>(pident_builder.Finish(&pident_array));
    std::shared_ptr<arrow::Array> pident_gap_array;
    static_cast<void>(pident_gap_builder.Finish(&pident_gap_array));
    std::shared_ptr<arrow::Array> evalue_array;
    static_cast<void>(evalue_builder.Finish(&evalue_array));
    std::shared_ptr<arrow::Array> length_array;
    static_cast<void>(length_builder.Finish(&length_array));
    std::shared_ptr<arrow::Array> qstart_array;
    static_cast<void>(qstart_builder.Finish(&qstart_array));
    std::shared_ptr<arrow::Array> qend_array;
    static_cast<void>(qend_builder.Finish(&qend_array));
    std::shared_ptr<arrow::Array> sstart_array;
    static_cast<void>(sstart_builder.Finish(&sstart_array));
    std::shared_ptr<arrow::Array> send_array;
    static_cast<void>(send_builder.Finish(&send_array));
    std::shared_ptr<arrow::Array> aln_len01_array;
    static_cast<void>(aln_len01_builder.Finish(&aln_len01_array));
    std::shared_ptr<arrow::Array> bitscore_array;
    static_cast<void>(bitscore_builder.Finish(&bitscore_array));
    std::shared_ptr<arrow::Array> score_array;
    static_cast<void>(score_builder.Finish(&score_array));
    std::shared_ptr<arrow::Array> qcovhsp_array;
    static_cast<void>(qcovhsp_builder.Finish(&qcovhsp_array));
    std::shared_ptr<arrow::Array> blast_score_array;
    static_cast<void>(blast_score_builder.Finish(&blast_score_array));
    std::shared_ptr<arrow::Array> gaps_array;
    static_cast<void>(gaps_builder.Finish(&gaps_array));
    std::shared_ptr<arrow::Array> nident_array;
    static_cast<void>(nident_builder.Finish(&nident_array));
    std::shared_ptr<arrow::Array> mismatch_array;
    static_cast<void>(mismatch_builder.Finish(&mismatch_array));
    std::shared_ptr<arrow::Array> positive_array;
    static_cast<void>(positive_builder.Finish(&positive_array));
    std::shared_ptr<arrow::Array> n_splices_array;
    static_cast<void>(n_splices_builder.Finish(&n_splices_array));
    std::shared_ptr<arrow::Array> hsp_cnt_array;
    static_cast<void>(hsp_cnt_builder.Finish(&hsp_cnt_array));
    std::shared_ptr<arrow::Array> sum_evalue_array;
    static_cast<void>(sum_evalue_builder.Finish(&sum_evalue_array));
    std::shared_ptr<arrow::Array> product_coverage_array;
    static_cast<void>(product_coverage_builder.Finish(&product_coverage_array));
    std::shared_ptr<arrow::Array> overall_identity_array;
    static_cast<void>(overall_identity_builder.Finish(&overall_identity_array));
    std::shared_ptr<arrow::Array> negative_count_array;
    static_cast<void>(negative_count_builder.Finish(&negative_count_array));
    std::shared_ptr<arrow::Array> matches_array;
    static_cast<void>(matches_builder.Finish(&matches_array));
    std::shared_ptr<arrow::Array> high_quality_percent_coverage_array;
    static_cast<void>(high_quality_percent_coverage_builder.Finish(&high_quality_percent_coverage_array));
    std::shared_ptr<arrow::Array> exon_identity_array;
    static_cast<void>(exon_identity_builder.Finish(&exon_identity_array));
    std::shared_ptr<arrow::Array> consensus_splices_array;
    static_cast<void>(consensus_splices_builder.Finish(&consensus_splices_array));
    std::shared_ptr<arrow::Array> comp_adj_method_array;
    static_cast<void>(comp_adj_method_builder.Finish(&comp_adj_method_array));
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({
  qhsp_array,
  shsp_array,
  pident_array,
  pident_gap_array,
  frames_array,
  evalue_array,
  length_array,
  aln_len01_array,
  qstart_array,
  qend_array,
  sstart_array,
  send_array,
  bitscore_array,
  score_array,
  qcovhsp_array,
  blast_score_array,
  gaps_array,
  nident_array,
  mismatch_array,
  positive_array,
  n_splices_array,
  hsp_cnt_array,
  sum_evalue_array,
  product_coverage_array,
  overall_identity_array,
  negative_count_array,
  matches_array,
  high_quality_percent_coverage_array,
  exon_identity_array,
  consensus_splices_array,
  comp_adj_method_array},
                                                                                                   {"qhsp", "shsp", "pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});
  
    assert(aln_struct_array.ok());
  
    std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
  
    // std::shared_ptr<arrow::Array> qseq_title_array;
    // static_cast<void>(qseq_title_builder.Finish(&qseq_title_array));
    // 
    // std::shared_ptr<arrow::Array> sseq_title_array;
    // static_cast<void>(sseq_title_builder.Finish(&sseq_title_array));
    
    std::shared_ptr<arrow::Array> qseqid_array;
    static_cast<void>(qseqid_builder.Finish(&qseqid_array));
  
    std::shared_ptr<arrow::Array> sseqid_array;
    static_cast<void>(sseqid_builder.Finish(&sseqid_array));
  
    std::shared_ptr<arrow::Array> qseq_array;
    static_cast<void>(qseq_builder.Finish(&qseq_array));
  
    std::shared_ptr<arrow::Array> sseq_array;
    static_cast<void>(sseq_builder.Finish(&sseq_array));
    
    std::shared_ptr<arrow::Array> qlen_array;
    static_cast<void>(qlen_builder.Finish(&qlen_array));
  
    std::shared_ptr<arrow::Array> slen_array;
    static_cast<void>(slen_builder.Finish(&slen_array));
  
    std::shared_ptr<arrow::Array> strand_array;
    static_cast<void>(strand_builder.Finish(&strand_array));
  
    std::shared_ptr<arrow::Array> num_alignment_array;
    static_cast<void>(num_alignments_builder.Finish(&num_alignment_array));
    
    // Create the seq_info struct array and populate with the arrays
    // std::shared_ptr<arrow::StructArray> seqtitle_struct_array = *arrow::StructArray::Make({qseq_title_array, sseq_title_array}, {arrow::field("qseq_title", arrow::utf8()), arrow::field("sseq_title", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
  
    arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
                                                                                                seqids_struct_array,
                                                                                                  seqs_struct_array,
                                                                                                  strand_array,
                                                                                                  lengths_struct_array},
                                                                                                 {"num_alignments", "seqids", "seqs", "strands", "lengths"});
  
    assert(seq_info_array.ok());
  
    std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
  
    // Rprintf("\n%d\n", num_rows); //DEBUG
    std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                                num_rows,
                                                                                {seq_info_array_, aln_struct_array_});
                   
     if(alignment_rb->num_rows() <= 0){
       Rcpp::stop("ExtractHits() - arrow::RecordBatch() - No alignments could be computed.");
     }
     
     arrow::Status align_sts = alignment_rb->ValidateFull();
     if(!align_sts.ok()){
       // Rcpp::Rcout << align_sts.message()  << std::endl << align_sts.ToString() << std::endl << "rows:" << alignment_rb->num_rows() << "\ncols:" << alignment_rb->num_columns()  << std::endl << std::flush; //DEBUG
       throw std::runtime_error("ExtractHits() - arrow::RecordBatch() - Alignments failed validation.");
     }                                                             
    
    if (alignment_rb)
    {
      // Rcpp::Rcout << alignment_rb->ToString() << std::endl << std::flush; //DEBUG
      // tmp_extracted++;
      // Rcpp::Rcout << "Extracted: " << tmp_extracted << std::endl << std::flush; //DEBUG
      
      const auto &wrt_sts = arrow_wrapper->AddRB2Batch(alignment_rb);
      if (!wrt_sts.ok())
      {
        throw runtime_error("ExtractHits() - Error adding RecordBatch to write buffer...");
      }
      
      progress_bar.increment();
      
      if(return_values){
        return alignment_rb;
      }else{
        return empty_rb;
      }
    }
    Rcpp::Rcerr << "ExtractHits() - Empty alignment_rb..." << std::endl << std::flush; //DEBUG
    return empty_rb; //arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // //ERROR RETURN, END
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("ExtractHits(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("ExtractHits(): C++ Runtime Error : ") + e.what());
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("ExtractHits() - Rcpp Exception : ") + e.what());
  }
  catch(...){
    Rcpp::stop("ExtractHits(): Unknown Exception");
  }
}


// // helper to get the molecular type for a Bioseq: returns true==protein, false==nucleotide
// bool QuickBLAST::Impl::IsProteinBioseq(const CBioseq_Handle &bh)
// {
//   if (!bh) return false;
//   try {
//     auto mol = bh.GetInst_Mol();
//     return (mol == CSeq_inst::eMol_aa);
//   } catch (...) {
//     return false;
//   }
// }

// // Extract full sequence string for a given Seq-id (requires scope to resolve)
// bool QuickBLAST::Impl::GetFullSequenceString(CRef<CSeq_id> id, CScope &scope, std::string &out_seq)
// {
//   if (!id) return false; //|| !scope
//   CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(*id);
//   if (!idh) return false;
//   CBioseq_Handle bh = scope.GetBioseqHandle(idh);
//   if (!bh) return false;
//   // determine coding (aa vs na)
//   CSeqVector sv(bh);
//   // Reserve to avoid repeated reallocations
//   out_seq.clear();
//   out_seq.reserve(bh.GetInst().GetLength() ? (size_t)bh.GetInst().GetLength() : 1024);
//   
//   // Get whole sequence
//   TSeqPos len = bh.GetInst().GetLength();
//   if (len == 0) return false;
//   
//   // CSeqVector::GetSeqData overloads differ by toolkit version; commonly:
//   // sv.GetSeqData(start, length, out_string, CSeqVector::eCoding_Iupac...)
//   // We'll use a robust loop to fetch in blocks to avoid relying on a single signature.
//   const TSeqPos block = 65536; // fetch in chunks if very large
//   for (TSeqPos pos = 0; pos < len; pos += block) {
//     TSeqPos fetch_len = std::min(block, len - pos);
//     std::string chunk;
//     // prefer IUPAC coding for DNA, IUPACaa for proteins
//     sv.GetSeqData(pos, fetch_len, chunk);
//     out_seq.append(chunk);
//   }
//   return true;
// }

// Return true on success, writes full sequence text to out_seq
bool QuickBLAST::Impl::GetFullSequenceString(CRef<CSeq_id> id, CScope &scope, std::string &out_seq)
{
  if (!id) return false;
  CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(*id);
  if (!idh) return false;
  
  CBioseq_Handle bh = scope.GetBioseqHandle(idh);
  if (!bh) return false;
  
  // If the CBioseq has explicit Seq-data in IUPAC, read it directly.
  if (bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_inst::TSeq_data &sd = bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    // If it's IUPACaa or IUPACna, extract directly:
    if (sd.IsIupacna()) {
      // sd.GetIupacna() returns a CSeq_data::TIupacna (string-like)
      out_seq.assign(sd.GetIupacna().Get());
      return true;
    } else if (sd.IsIupacaa()) {
      out_seq.assign(sd.GetIupacaa().Get());
      return true;
    }
    // otherwise fall back to SeqVector (it will do conversions)
  }
  
  // Fallback: use CSeqVector and force IUPAC coding so we get readable letters.
  CSeqVector sv(bh);
  sv.SetIupacCoding();                 // <-- critical fix
  out_seq.clear();
  TSeqPos len = bh.GetInst().GetLength();
  if (len == 0) return false;
  
  // Fetch in chunks (use pos, stop overload — stop is inclusive)
  const TSeqPos block = 65536;
  for (TSeqPos pos = 0; pos < len; pos += block) {
    TSeqPos fetch_len = std::min(block, len - pos);
    TSeqPos stop = pos + fetch_len - 1;
    std::string chunk;
    sv.GetSeqData(pos, stop, chunk); // returns ASCII IUPAC letters now
    out_seq.append(chunk);
  }
  return true;
}


// // // Extract HSP sequences (ungapped concatenation) and optionally aligned-with-gaps strings
// // // returns true on success
// // bool QuickBLAST::Impl::GetHSPSequencesFromDenseg(const CDense_seg& dseg, CScope &scope,
// //                                       std::string &q_hsp_ungapped, std::string &s_hsp_ungapped,
// //                                       std::string *q_aligned_with_gaps, std::string *s_aligned_with_gaps)
// // {
// //   // if (!scope) { Rcpp::Rcout << "here 0.1" << std::endl << std::flush; return false; }
// //   q_hsp_ungapped.clear();
// //   s_hsp_ungapped.clear();
// //   if (q_aligned_with_gaps) q_aligned_with_gaps->clear();
// //   if (s_aligned_with_gaps) s_aligned_with_gaps->clear();
// //   
// //   // rows and segments
// //   const size_t num_rows = dseg.CheckNumRows();
// //   const size_t num_segs = dseg.GetNumseg();            
// //   if (num_rows < 2 || num_segs == 0) { Rcpp::Rcout << "here 0.2" << std::endl << std::flush; return false; }
// //   
// //   // we assume row 0 = query, row 1 = subject (common for pairwise BLAST)
// //   const size_t q_row = 0;
// //   const size_t s_row = 1;
// //   
// //   // Seq ids
// //   if (!dseg.CanGetIds()) { Rcpp::Rcout << "here 0.3" << std::endl << std::flush; return false; }
// //   const auto &ids = dseg.GetIds();
// //   if (ids.size() <= q_row || ids.size() <= s_row) { Rcpp::Rcout << "here 0.4" << std::endl << std::flush; return false; }
// //   
// //   // prepare SeqVectors for each row
// //   CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(*ids[q_row]);
// //   CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(*ids[s_row]);
// //   if (!q_idh || !s_idh) { Rcpp::Rcout << "here 0.5" << std::endl << std::flush; return false; }
// //   Rcpp::Rcout << "here 0.5.1 : q_idh : " << q_idh.GetSeqId()->GetSeqIdString(true) << std::endl << std::flush;
// //   Rcpp::Rcout << "here 0.5.2 : s_idh : " << s_idh.GetSeqId()->GetSeqIdString(true) << std::endl << std::flush;
// //   CRef<CSeq_id> qid_copy(new CSeq_id(q_idh.GetSeqId()->AsFastaString(), CSeq_id::fParse_AnyRaw | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal)); 
// //   CRef<CSeq_id> sid_copy(new CSeq_id(s_idh.GetSeqId()->AsFastaString(), CSeq_id::fParse_AnyRaw | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal)); 
// //   // CRef<CSeq_id> qid_copy(new CSeq_id(*ids[q_row])); // copy ctor
// //   // CRef<CSeq_id> sid_copy(new CSeq_id(*ids[s_row])); // copy ctor
// //   qid_copy->SetLocal().SetStr( q_idh.GetSeqId()->GetSeqIdString(true) );
// //   sid_copy->SetLocal().SetStr( s_idh.GetSeqId()->GetSeqIdString(true) );
// //   Rcpp::Rcout << "here 0.5.3 : qid_copy(LOCAL) : " << qid_copy->GetSeqIdString(true) << std::endl << std::flush;
// //   Rcpp::Rcout << "here 0.5.4 : sid_copy(LOCAL) : " << sid_copy->GetSeqIdString(true) << std::endl << std::flush;
// //   Rcpp::Rcout << "here 0.5.5 : qid_copy(LOCAL) : " << qid_copy->AsFastaString() << std::endl << std::flush;
// //   Rcpp::Rcout << "here 0.5.6 : sid_copy(LOCAL) : " << sid_copy->AsFastaString() << std::endl << std::flush;
// //   CBioseq_Handle q_bh = scope.GetBioseqHandle(CSeq_id_Handle::GetHandle(*qid_copy)); //*ids[q_row] //q_idh
// //   CBioseq_Handle s_bh = scope.GetBioseqHandle(CSeq_id_Handle::GetHandle(*sid_copy)); //*ids[s_row] //s_idh
// //   if (!q_bh || !s_bh) {
// //     Rcpp::Rcout << "here 0.6.1 : " << q_bh.GetSeqId()->AsFastaString() << std::endl << std::flush;
// //     Rcpp::Rcout << "here 0.6.2 : " << s_bh.GetSeqId()->AsFastaString() << std::endl << std::flush;
// //     
// //     vector<CSeq_entry_Handle> tse_handles;
// //     scope.GetAllTSEs(tse_handles);        // fills vector with all top-level seq-entry handles
// //     Rcpp::Rcout << tse_handles.size() << std::endl << std::flush;
// //     for (const CSeq_entry_Handle &seh : tse_handles) {
// //       // CBioseq_CI enumerates all Bioseqs inside the seq-entry (one level or recursive
// //       // depending on flags; default is to enumerate the bioseqs in the entry)
// //       for (CBioseq_CI bs_it(seh); bs_it; ++bs_it) {
// //         CBioseq_Handle bh = *bs_it;                    // get handle to this Bioseq
// //         CConstRef<CBioseq> seq = bh.GetCompleteBioseq(); // load full object
// //         
// //         // print all Seq-ids for this Bioseq
// //         for (const CRef<CSeq_id> &id : seq->GetId()) {
// //           NcbiCout << "Bioseq: " << id->GetSeqIdString(true) << NcbiEndl;
// //         }
// //       }
// //     }
// //     return false;
// //   }
// //   
// //   Rcpp::Rcout << q_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data() << std::endl << std::flush; //DEBUG
// //   Rcpp::Rcout << s_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data() << std::endl << std::flush; //DEBUG
// //   
// //   CSeqVector qsv(q_bh);
// //   CSeqVector ssv(s_bh);
// //   
// //   // denseg stores starts as a flat vector: starts[seg_index * num_rows + row]
// //   // and lengths as lengths[seg_index]
// //   const auto &starts = dseg.GetStarts();
// //   const auto &lens = dseg.GetLens();
// //   
// //   // iterate segments
// //   for (size_t seg = 0; seg < num_segs; ++seg) {
// //     TSeqPos seg_len = lens[seg];
// //     // compute indexes into starts vector
// //     ssize_t q_start = starts[seg * num_rows + q_row];
// //     ssize_t s_start = starts[seg * num_rows + s_row];
// //     
// //     // ungapped concatenation: append only real residues (skip gaps)
// //     if (q_start >= 0 && seg_len > 0) {
// //       TSeqPos q_stop = (TSeqPos)q_start + (TSeqPos)seg_len - 1;
// //       std::string qchunk;
// //       qsv.GetSeqData((TSeqPos)q_start, q_stop, qchunk);
// //       q_hsp_ungapped.append(qchunk);
// //       Rcpp::Rcout << "Ungapped qchunk size:" << qchunk.size() << std::endl << std::flush; //DEBUG
// //     }
// //     if (s_start >= 0 && seg_len > 0) {
// //       TSeqPos s_stop = (TSeqPos)s_start + (TSeqPos)seg_len - 1;
// //       std::string schunk;
// //       ssv.GetSeqData((TSeqPos)s_start, s_stop, schunk);
// //       s_hsp_ungapped.append(schunk);
// //       Rcpp::Rcout << "Ungapped schunk size:" << schunk.size() << std::endl << std::flush; //DEBUG
// //     }
// //     
// //     // aligned-with-gaps (optional): insert '-' for gaps
// //     if (q_aligned_with_gaps || s_aligned_with_gaps) {
// //       if (q_aligned_with_gaps) {
// //         if (q_start >= 0) {
// //           TSeqPos q_stop = (TSeqPos)q_start + (TSeqPos)seg_len - 1;
// //           std::string qchunk;
// //           qsv.GetSeqData((TSeqPos)q_start, q_stop, qchunk);
// //           q_aligned_with_gaps->append(qchunk);
// //           Rcpp::Rcout << "Gapped qchunk size:" << qchunk.size() << std::endl << std::flush; //DEBUG
// //         } else {
// //           q_aligned_with_gaps->append(seg_len, '-');
// //         }
// //       }
// //       if (s_aligned_with_gaps) {
// //         if (s_start >= 0) {
// //           TSeqPos s_stop = (TSeqPos)s_start + (TSeqPos)seg_len - 1;
// //           std::string schunk;
// //           ssv.GetSeqData((TSeqPos)s_start, s_stop, schunk);
// //           s_aligned_with_gaps->append(schunk);
// //           Rcpp::Rcout << "Gapped schunk size:" << schunk.size() << std::endl << std::flush; //DEBUG
// //         } else {
// //           s_aligned_with_gaps->append(seg_len, '-');
// //         }
// //       }
// //     }
// //   }
// //   
// //   NcbiCout << "dseg rows=" << num_rows << " segs=" << num_segs
// //            << " starts=" << starts.size() << " lens=" << lens.size() << NcbiEndl;
// //   for (size_t seg=0; seg < std::min(num_segs, (size_t)10); ++seg) {
// //     NcbiCout << "seg " << seg << " len=" << lens[seg]
// //              << " qstart=" << starts[seg * num_rows + 0]
// //              << " sstart=" << starts[seg * num_rows + 1] << NcbiEndl;
// //   }
// //   return true;
// // }
// 
// bool QuickBLAST::Impl::GetHSPSequencesFromDenseg(const CDense_seg& dseg, CScope &scope,
//                                                  std::string &q_hsp_ungapped, std::string &s_hsp_ungapped,
//                                                  std::string *q_aligned_with_gaps, std::string *s_aligned_with_gaps)
// {
//   q_hsp_ungapped.clear();
//   s_hsp_ungapped.clear();
//   if (q_aligned_with_gaps) q_aligned_with_gaps->clear();
//   if (s_aligned_with_gaps) s_aligned_with_gaps->clear();
//   
//   const size_t num_rows = dseg.CheckNumRows();
//   const size_t num_segs = dseg.GetNumseg();
//   if (num_rows < 2 || num_segs == 0) { Rcpp::Rcout << "here 0.2" << std::endl; return false; }
//   
//   const size_t q_row = 0;
//   const size_t s_row = 1;
//   
//   if (!dseg.CanGetIds()) { Rcpp::Rcout << "here 0.3" << std::endl; return false; }
//   const auto &ids = dseg.GetIds();
//   if (ids.size() <= q_row || ids.size() <= s_row) { Rcpp::Rcout << "here 0.4" << std::endl; return false; }
//   
//   CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(*ids[q_row]);
//   CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(*ids[s_row]);
//   if (!q_idh || !s_idh) { Rcpp::Rcout << "here 0.5" << std::endl; return false; }
//   
//   // Try to find bioseq handles in scope (use handle created consistently)
//   CBioseq_Handle q_bh = scope.GetBioseqHandle(q_idh);
//   CBioseq_Handle s_bh = scope.GetBioseqHandle(s_idh);
//   if (!q_bh || !s_bh) {
//     // If direct handle fails, try tolerant fallback (token parse, local id) -- see earlier ResolveBioseqHandleFromAlignId
//     Rcpp::Rcout << "Fallback: direct handle failed for ids: "
//                 << q_idh.AsString() << " / " << s_idh.AsString() << std::endl;
//     return false;
//   }
//   
//   // Debug: extract stored sequence string
//   std::string stored_q, stored_s;
//   if (q_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
//     const CSeq_data &sd = q_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
//     if (sd.IsIupacaa()) stored_q = sd.GetIupacaa().Get();
//     else if (sd.IsIupacna()) stored_q = sd.GetIupacna().Get();
//   }
//   if (s_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
//     const CSeq_data &sd = s_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
//     if (sd.IsIupacaa()) stored_s = sd.GetIupacaa().Get();
//     else if (sd.IsIupacna()) stored_s = sd.GetIupacna().Get();
//   }
//   Rcpp::Rcout << "Stored q len: " << stored_q.size() << " s len: " << stored_s.size() << std::endl;
//   
//   CSeqVector qsv(q_bh);
//   CSeqVector ssv(s_bh);
//   
//   const auto &starts = dseg.GetStarts();
//   const auto &lens = dseg.GetLens();
//   
//   for (size_t seg = 0; seg < num_segs; ++seg) {
//     TSeqPos seg_len = lens[seg];
//     ssize_t q_start = starts[seg * num_rows + q_row];
//     ssize_t s_start = starts[seg * num_rows + s_row];
//     
//     if (seg_len == 0) continue;
//     
//     // Preferred: use GetSeqData(from, length, out)
//     if (q_start >= 0) {
//       std::string qchunk;
//       qsv.GetSeqData((TSeqPos)q_start, (TSeqPos)seg_len, qchunk);
//       Rcpp::Rcout << "Ungapped qchunk size:" << qchunk.size() << std::endl;
//       q_hsp_ungapped.append(qchunk);
//       
//       // debug compare with stored string (if available)
//       if (stored_q.size() >= (size_t)q_start + seg_len) {
//         std::string expected = stored_q.substr(q_start, seg_len);
//         if (expected != qchunk) {
//           // print tiny hex-dump for quick diagnosis
//           auto hd = [](const std::string &s, size_t N = 32) {
//             std::ostringstream o; o<<std::hex<<std::setfill('0');
//             for (size_t i=0;i<std::min(N,s.size());++i) o << std::setw(2) << (int)(unsigned char)s[i] << " ";
//             return o.str();
//           };
//           Rcpp::Rcout << "Mismatch in qchunk bytes. expected hex: " << hd(expected,40)
//                       << "\n got hex: " << hd(qchunk,40) << std::endl;
//         }
//       }
//     } else {
//       if (q_aligned_with_gaps) q_aligned_with_gaps->append(seg_len, '-');
//       if (q_start < 0) Rcpp::Rcout << "q_start negative => gap for " << seg_len << std::endl;
//     }
//     
//     if (s_start >= 0) {
//       std::string schunk;
//       ssv.GetSeqData((TSeqPos)s_start, (TSeqPos)seg_len, schunk);
//       Rcpp::Rcout << "Ungapped schunk size:" << schunk.size() << std::endl;
//       s_hsp_ungapped.append(schunk);
//       
//       if (stored_s.size() >= (size_t)s_start + seg_len) {
//         std::string expected_s = stored_s.substr(s_start, seg_len);
//         if (expected_s != schunk) {
//           Rcpp::Rcout << "Mismatch in schunk bytes (hex)..." << std::endl;
//         }
//       }
//     } else {
//       if (s_aligned_with_gaps) s_aligned_with_gaps->append(seg_len, '-');
//     }
//     
//     // aligned-with-gaps building (same pattern)
//     if (q_aligned_with_gaps) {
//       if (q_start >= 0) {
//         std::string qchunk;
//         qsv.GetSeqData((TSeqPos)q_start, (TSeqPos)seg_len, qchunk);
//         q_aligned_with_gaps->append(qchunk);
//         Rcpp::Rcout << "Gapped qchunk size:" << qchunk.size() << std::endl;
//       } else q_aligned_with_gaps->append(seg_len, '-');
//     }
//     if (s_aligned_with_gaps) {
//       if (s_start >= 0) {
//         std::string schunk;
//         ssv.GetSeqData((TSeqPos)s_start, (TSeqPos)seg_len, schunk);
//         s_aligned_with_gaps->append(schunk);
//         Rcpp::Rcout << "Gapped schunk size:" << schunk.size() << std::endl;
//       } else s_aligned_with_gaps->append(seg_len, '-');
//     }
//   }
//   
//   return true;
// }

bool QuickBLAST::Impl::GetHSPSequencesFromDenseg(const CDense_seg& dseg, CScope &scope,
                                                 std::string &q_hsp_ungapped, std::string &s_hsp_ungapped,
                                                 std::string *q_aligned_with_gaps, std::string *s_aligned_with_gaps)
{
  q_hsp_ungapped.clear();
  s_hsp_ungapped.clear();
  if (q_aligned_with_gaps) q_aligned_with_gaps->clear();
  if (s_aligned_with_gaps) s_aligned_with_gaps->clear();
  
  const size_t num_rows = dseg.CheckNumRows();
  const size_t num_segs = dseg.GetNumseg();
  if (num_rows < 2 || num_segs == 0) { return false; }
  
  const size_t q_row = 0, s_row = 1;
  if (!dseg.CanGetIds()) { return false; }
  const auto &ids = dseg.GetIds();
  if (ids.size() <= q_row || ids.size() <= s_row) { return false; }
  
  CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(*ids[q_row]);
  CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(*ids[s_row]);
  if (!q_idh || !s_idh) { return false; }
  
  CBioseq_Handle q_bh = scope.GetBioseqHandle(q_idh);
  CBioseq_Handle s_bh = scope.GetBioseqHandle(s_idh);
  if (!q_bh || !s_bh) {
    // Rcpp::Rcout << "Fallback: direct handle failed for ids: "
                // << q_idh.AsString() << " / " << s_idh.AsString() << std::endl;
    return false;
  }
  
  // Extract the stored IUPAC strings (if present) — prefer these
  std::string stored_q, stored_s;
  if (q_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_data &sd = q_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    if (sd.IsIupacaa()) stored_q = sd.GetIupacaa().Get();
    else if (sd.IsIupacna()) stored_q = sd.GetIupacna().Get();
  }
  if (s_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_data &sd = s_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    if (sd.IsIupacaa()) stored_s = sd.GetIupacaa().Get();
    else if (sd.IsIupacna()) stored_s = sd.GetIupacna().Get();
  }
  // Rcpp::Rcout << "Stored q len: " << stored_q.size() << " s len: " << stored_s.size() << std::endl << std::flush;
  
  CSeqVector qsv(q_bh);
  CSeqVector ssv(s_bh);
  
  qsv.SetIupacCoding(); 
  ssv.SetIupacCoding(); 
  
  const auto &starts = dseg.GetStarts();
  const auto &lens = dseg.GetLens();
  
  for (size_t seg = 0; seg < num_segs; ++seg) {
    TSeqPos seg_len = lens[seg];
    if (seg_len == 0) continue;
    
    ssize_t q_start = starts[seg * num_rows + q_row];
    ssize_t s_start = starts[seg * num_rows + s_row];
    
    // First try: use stored string for exact ASCII residues (recommended)
    if (q_start >= 0 && stored_q.size() >= (size_t)q_start + seg_len) {
      std::string qchunk = stored_q.substr((size_t)q_start, (size_t)seg_len);
      q_hsp_ungapped.append(qchunk);
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(qchunk);
      // Rcpp::Rcout << "Ungapped qchunk size:" << qchunk.size() << std::endl << std::flush;
    } else if (q_start >= 0) {
      // fallback: pull from CSeqVector (might be encoded; handle with care)
      std::string qchunk;
      qsv.GetSeqData((TSeqPos)q_start, (TSeqPos)seg_len, qchunk); // from,length
      q_hsp_ungapped.append(qchunk);
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(qchunk);
      // Rcpp::Rcout << "Ungapped qchunk (from CSeqVector) size:" << qchunk.size() << std::endl << std::flush;
    } else {
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(seg_len, '-');
    }
    
    if (s_start >= 0 && stored_s.size() >= (size_t)s_start + seg_len) {
      std::string schunk = stored_s.substr((size_t)s_start, (size_t)seg_len);
      s_hsp_ungapped.append(schunk);
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(schunk);
      // Rcpp::Rcout << "Ungapped schunk size:" << schunk.size() << std::endl << std::flush;
    } else if (s_start >= 0) {
      std::string schunk;
      ssv.GetSeqData((TSeqPos)s_start, (TSeqPos)seg_len, schunk);
      s_hsp_ungapped.append(schunk);
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(schunk);
      // Rcpp::Rcout << "Ungapped schunk (from CSeqVector) size:" << schunk.size() << std::endl << std::flush;
    } else {
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(seg_len, '-');
    }
  }
  
  return true;
}


/*
 * Utility: AddAllAvailableScoresToAlign
 * Computes and attaches a large set of named scores for a single CSeq_align.
 *
 * Parameters:
 *   align  - CRef<CSeq_align> (the alignment to augment)
 *   scope  - CRef<CScope> containing the Bioseqs referenced by the align (must be same scope used elsewhere)
 *   effective_search_space - optional: set to non-zero to configure evalue/bit-score computations
 *
 * After calling this, many calls to align->GetNamedScore(...) will return true.
 */
void QuickBLAST::Impl::AddAllAvailableScoresToAlign(CRef<CSeq_align> align,
                                         CRef<CScope> scope,
                                         double effective_search_space = 0.0)
{
  if (!align || !scope) {
    return;
  }
  
  try {
    CScoreBuilder scorer;
    
    if (effective_search_space > 0.0) {
      // configuring effective search space can produce correct e-values in some contexts
      scorer.SetEffectiveSearchSpace(effective_search_space);
    }
    
    // Common/core scores
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_AlignLength);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_BitScore);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_Blast);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_PercentIdentity);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_GapCount);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_EValue);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_IdentityCount);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_MismatchCount);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_PercentCoverage);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_Score);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_PositiveCount);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_Splices);
    
    // Extended scores
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_SumEValue);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_ProductCoverage);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_OverallIdentity);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_NegativeCount);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_Matches);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_HighQualityPercentCoverage);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_ExonIdentity);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_ConsensusSplices);
    scorer.AddScore(*scope, *align, CSeq_align::EScoreType::eScore_CompAdjMethod);
    
    // For spliced aligners, add Splign-specific scores if relevant
    try {
      scorer.AddSplignScores(*align);
    } catch (const CException& e) {
      // Not fatal — only some alignments need splign scores
      LOG_POST(Warning << "AddSplignScores() failed or not applicable: " << e.GetMsg());
    }
    
    // Fallback: ensure PercentIdentity is present; if not compute from IdentityCount / AlignLength
    double pident = 0.0;
    int aln_len = 0;
    int nident = 0;
    bool has_pident = align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident);
    bool has_alnlen = align->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
    bool has_nident = align->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, nident);
    
    if (!has_pident && has_alnlen && has_nident && aln_len > 0) {
      double computed = 100.0 * double(nident) / double(aln_len);
      // set the named score explicitly on the CSeq_align:
      align->SetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, computed);
    }
    
  } catch (const CException& e) {
    // don't let toolkit exceptions propagate uncontrolled into R; log and continue
    ERR_POST(Warning << "AddAllAvailableScoresToAlign(): NCBI Exception: " << e.GetMsg());
  } catch (const std::exception& e) {
    ERR_POST(Warning << "AddAllAvailableScoresToAlign(): std::exception: " << e.what());
  }catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("AddAllAvailableScoresToAlign(): Rcpp Exception : ") + e.what());
  }
  catch (...) {
    ERR_POST(Warning << "AddAllAvailableScoresToAlign(): unknown exception");
  }
}


/*
 * Utility: AddAllAvailableScoresToAlignList
 * Adds scores to every alignment in a list/vector.
 */
void QuickBLAST::Impl::AddAllAvailableScoresToAlignList(std::list< CRef<CSeq_align> >& aligns, CRef<CScope> scope, double effective_search_space = 0.0)
{
  if (!scope) return;
  try {
    CScoreBuilder scorer;
    if (effective_search_space > 0.0) scorer.SetEffectiveSearchSpace(effective_search_space);
    
    // Compute batch scores (AddScore has an overload for list)
    // We'll ask for a set of scores in a loop to leverage internal batching
    std::vector<CSeq_align::EScoreType> score_types = {
      CSeq_align::EScoreType::eScore_AlignLength,
      CSeq_align::EScoreType::eScore_BitScore,
      CSeq_align::EScoreType::eScore_Blast,
      CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped,
      CSeq_align::EScoreType::eScore_PercentIdentity,
      CSeq_align::EScoreType::eScore_GapCount,
      CSeq_align::EScoreType::eScore_EValue,
      CSeq_align::EScoreType::eScore_IdentityCount,
      CSeq_align::EScoreType::eScore_MismatchCount,
      CSeq_align::EScoreType::eScore_PercentCoverage,
      CSeq_align::EScoreType::eScore_Score,
      CSeq_align::EScoreType::eScore_PositiveCount,
      CSeq_align::EScoreType::eScore_Splices,
      CSeq_align::EScoreType::eScore_SumEValue,
      CSeq_align::EScoreType::eScore_ProductCoverage,
      CSeq_align::EScoreType::eScore_OverallIdentity,
      CSeq_align::EScoreType::eScore_NegativeCount,
      CSeq_align::EScoreType::eScore_Matches,
      CSeq_align::EScoreType::eScore_HighQualityPercentCoverage,
      CSeq_align::EScoreType::eScore_ExonIdentity,
      CSeq_align::EScoreType::eScore_ConsensusSplices,
      CSeq_align::EScoreType::eScore_CompAdjMethod
    };
    
    for (auto st : score_types) {
      try {
        scorer.AddScore(*scope, aligns, st);
      } catch (const CException& e) {
        // non-fatal; continue with others
        ERR_POST(Warning << "AddScore for type " << static_cast<int>(st) << " failed: " << e.GetMsg());
      }
    }
    
    // Splign scores for each alignment if needed
    for (CRef<CSeq_align> &a : aligns) {
      try {
        scorer.AddSplignScores(*a);
      } catch (...) { /* ignore */ }
    }
    
    // compute percent identity fallback per alignment if missing
    for (CRef<CSeq_align> &a : aligns) {
      if (!a) continue;
      int aln_len = 0;
      int nident = 0;
      double pident = 0.0;
      bool hasp = a->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident);
      bool haslen = a->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
      bool hasid = a->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, nident);
      if (!hasp && haslen && hasid && aln_len > 0) {
        double computed = 100.0 * double(nident) / double(aln_len);
        a->SetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, computed);
      }
    }
    
  } catch (const CException& e) {
    ERR_POST(Error << "AddAllAvailableScoresToAlignList: NCBI Exception: " << e.GetMsg());
  }
}


/*
 * Convenience overloads for common container types
 */

// For CRef<CSeq_align_set>
void QuickBLAST::Impl::AddAllAvailableScoresToAlignSet(CRef<CSeq_align_set> alnset, CRef<CScope> scope, double effective_search_space = 0.0)
{
  if (!alnset || !scope) return;
  
  // Convert align set to list
  std::list< CRef<CSeq_align> > aligns;
  if (alnset->IsSet()) {
    for (auto &r : alnset->Get()) {
      aligns.push_back(r);
    }
  }
  
  AddAllAvailableScoresToAlignList(aligns, scope, effective_search_space);
  
  // The CSeq_align objects inside alnset are updated (we modified the same CRef objects)
}

// For TSeqAlignVector (vector<CRef<CSeq_align>>)
void QuickBLAST::Impl::AddAllAvailableScoresToSeqAlignVector(TSeqAlignVector &alnvec, CRef<CScope> scope, double effective_search_space = 0.0)
{
  std::list< CRef<CSeq_align> > aligns;
  for (const auto& r : alnvec) {
    if (!r) continue;
    for (const auto& aln : r->Get()) {
      aligns.push_back(aln);
    }
  }
  AddAllAvailableScoresToAlignList(aligns, scope, effective_search_space);
}



// Returns empty string on failure.
// std::string QuickBLAST::Impl::GetTitleForSeqId(const CSeq_id &seq_id, CRef<CScope> scope) {
//   try {
//     CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(seq_id);
//     CBioseq_Handle bh = scope->GetBioseqHandle(idh);     // may throw if not present
//     // Prefer title descriptor
//     CSeqdesc_CI tit(bh, CSeqdesc::e_Title);
//     if (tit) {
//       const CSeqdesc &d = *tit;
//       if (d.IsTitle() && !d.GetTitle().Get().empty()) {
//         return d.GetTitle().Get().front();
//       }
//     }
//     // Fallback: use the first Seq-id on the Bioseq as display string
//     const CBioseq &b = bh.GetCompleteObject();
//     if (b.IsSetId() && !b.GetId().empty()) {
//       return b.GetId().front()->AsFastaString();
//     }
//   } catch (const CException &e) {
//     // scope couldn't resolve seq_id or other ObjMgr error
//     // return empty or log error
//   }
//   // As last resort, return the Seq-id itself as a string
//   return seq_id.AsFastaString();
// }
// 
// // Returns pair<query_title, subject_title>
// std::pair<std::string, std::string> QuickBLAST::Impl::GetTitlesFromSeqAlign(const CSeq_align &align, CRef<CScope> scope) {
//   std::string qtitle, stitle;
//   
//   // some align types store ids per row (GetDim() rows)
//   int dim = align.GetDim();  // number of rows
//   if (dim >= 1) {
//     const CSeq_id& qid = align.GetSeq_id(0);
//     qtitle = GetTitleForSeqId(qid, scope);
//   }
//   if (dim >= 2) {
//     const CSeq_id& sid = align.GetSeq_id(1);
//     stitle = GetTitleForSeqId(sid, scope);
//   }
//   // If dim < 2, you may need to inspect the segs or other fields
//   return {qtitle, stitle};
// }


// template <>
// SSeqLoc *QuickBLAST::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope)
// {
//   int rec_no = fasta_data.rec_no;
//   std::string fastaID(fasta_data.header.data());
//   std::string fastaSequence(fasta_data.seq.data());

//   const TSeqPos seqlen = fastaSequence.length();

//   _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
//   ncbi::CRef<ncbi::objects::CSeq_interval> interval(new ncbi::objects::CSeq_interval());
//   interval->SetFrom(0);
//   interval->SetTo(seqlen - 1);

//   CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
//   id->Select(CSeq_id_Base::E_Choice::e_Local);
//   id->SetLocal().SetId(rec_no);
//   id->SetLocal().SetStr(fastaID);

//   CRef<CSeq_loc>
//       cseq_loc_obj(new CSeq_loc());
//   cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
//   cseq_loc_obj->SetWhole()
//       .SetLocal()
//       .SetStr(fastaID);
//   if (seq_type == ESeqType::eProtein)
//   {
//     cseq_loc_obj->SetStrand(eNa_strand_unknown);
//   }
//   else
//   {
//     switch (strand)
//     {
//     case EStrand::ePlus:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_plus);
//       break;
//     case EStrand::eMinus:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_minus);
//       break;
//     case EStrand::eUnknown:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_unknown);
//       break;
//     case EStrand::eBoth:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both);
//       break;
//     case EStrand::eBoth_rev:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both_rev);
//       break;
//     case EStrand::eOther:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_other);
//       break;
//     }
//   }

//   CRef<CSeq_data> seq_data(new CSeq_data());
//   seq_data->Select(seq_type == ESeqType::eProtein ? CSeq_data_Base::E_Choice::e_Iupacaa : CSeq_data_Base::E_Choice::e_Iupacna);
//   switch (seq_type)
//   {
//   case ESeqType::eProtein:
//   {
//     seq_data->SetIupacaa(CIUPACaa(fastaSequence));
//     break;
//   }
//   case ESeqType::eNucleotide:
//   {
//     seq_data->SetIupacna(CIUPACna(fastaSequence));
//     break;
//   }
//   }

//   CRef<CSeq_inst> seq_inst(new CSeq_inst());
//   seq_inst->SetSeq_data(*seq_data);
//   seq_inst->SetLength(fastaSequence.length());
//   seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
//   seq_inst->SetTopology(CSeq_inst::eTopology_linear);
//   seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
//   seq_inst->SetRepr(CSeq_inst_Base::ERepr::eRepr_raw);
//   seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
//   seq_inst->SetTopology(CSeq_inst_Base::ETopology::eTopology_linear);
//   seq_inst->SetStrand(CSeq_inst_Base::EStrand::eStrand_ss);
//   seq_inst->SetLength(seqlen);

//   CRef<ncbi::objects::CBioseq> bioseq(new CBioseq(*cseq_loc_obj, fastaID));
//   bioseq->SetInst(*seq_inst);

//   CRef<CSeq_entry>
//       ret_entry(new CSeq_entry());
//   ret_entry->SetSeq(*bioseq);

//   parent_scope->AddTopLevelSeqEntry(*ret_entry);

//   return new SSeqLoc(cseq_loc_obj.GetObject(), parent_scope.GetObject());
// }

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
{
  return pImpl->BLAST_files(queryFile, subjectFile, outFile, outFormat, blast_sequence_limit, num_threads, show_progress, return_values, batch_size);
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject)
{
  return pImpl->BLAST_seqs(query, subject);
}

auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
{
  return pImpl->BLAST(query, subject, outputFile, outFormat, input_type, blast_sequence_limit, show_progress);
}

QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences, bool save_hsp_sequences)
    : pImpl(std::make_unique<Impl>(seq_type, strand, program, options, save_sequences, save_hsp_sequences)) {}

// Destructor: default implementation (pImpl will automatically clean up)
QuickBLAST::~QuickBLAST() = default;



