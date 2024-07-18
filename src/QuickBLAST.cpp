// [[Rcpp::plugins(openmp)]]
#define CALLDEF(name, n)       \
  {                            \
    #name, (DL_FUNC) & name, n \
  }

// #include <cassert>
// #include <chrono>
#include <stdlib.h>
// #include <experimental/filesystem>
#include <filesystem>
// #include <functional>

#include <R_ext/Rdynload.h>
// #include <R_ext/Visibility.h>

#include <RcppCommon.h>
#include <Rcpp.h>

#include <cassert>
#include <chrono>
#include <functional>

#include <map>
#include <tuple>

// #include <algo\blast\QuickBLAST\QuickBLAST.hpp>
using namespace Rcpp;

std::map<unsigned int, std::tuple<int, int, int, std::string, std::string, bool>> obj_list;

#if defined(MINGW32) || defined(WIN32)
#ifdef QBLIBRARY_IMPORTS
#define QBLIBRARY_API __declspec(dllimport)
#endif
#else
#define QBLIBRARY_API //__declspec(dllimport)
#endif

RCPP_EXPOSED_CLASS(QuickBLAST);
RCPP_EXPOSED_CLASS(ArrowWrapper);

// RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::ESeqType);
// RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::EStrand);
// RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::EInputType);

enum ESeqType
{
  eNucleotide = 0,
  eProtein = 1
};
enum EStrand
{
  ePlus = 0,
  eMinus = 1,
  eBoth = 2,
  eBoth_rev = 3,
  eOther = 4,
  eUnknown = 5
};
enum EInputType
{
  eFile = 0,
  eSequenceString = 1,
  eFolder = 2
};

extern "C"
{
  /*QuickBLASTHandle GetQuickBLASTInstance(int id);
   int GetInstanceCount();
   int GetInstanceID(std::shared_ptr<QuickBLAST> ptr);
   QuickBLASTHandle CreateQuickBLASTInstance(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences = false);
   void DeleteQuickBLASTInstance(QuickBLAST* ptr);
   //ArrowWrapper* CreateArrowWrapperInstance();
   std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
   std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
   void SetQuickBLASTOptions(std::shared_ptr<QuickBLAST> ptr, std::string program_name, std::string options);*/

  // QBLIBRARY_API QuickBLASTHandle GetQuickBLASTInstance(int id);
  // QBLIBRARY_API int GetInstanceCount();
  // QBLIBRARY_API int GetInstanceID(QuickBLASTHandle ptr);
  // QBLIBRARY_API QuickBLASTHandle CreateQuickBLASTInstance(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences = false);
  // QBLIBRARY_API void DeleteQuickBLASTInstance(QuickBLASTHandle ptr);
  // // ArrowWrapper* CreateArrowWrapperInstance();
  // // std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
  // // std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
  // QBLIBRARY_API void SetQuickBLASTOptions(QuickBLASTHandle ptr, std::string program_name, std::string options);
  // QBLIBRARY_API ArrowRBVHandle BLAST2Seqs(QuickBLASTHandle ptr, std::string query, std::string subject);
  // QBLIBRARY_API ArrowRBVHandle BLAST2Files(QuickBLASTHandle ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
  // QBLIBRARY_API void BLAST2Folders(QuickBLASTHandle ptr, std::string query, std::string subject, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
  // QBLIBRARY_API void BLAST1Folder(QuickBLASTHandle ptr, std::string input_folder, std::string extension, std::string out_folder, int num_threads, bool reciprocal_hits, int min_batch_size = 1024);
  // QBLIBRARY_API std::string getFilenameWithoutExtension(const std::string &filename);
  // // QBLIBRARY_API List rm_null(List x);
  // QBLIBRARY_API std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension);
  // // QBLIBRARY_API void R_init_QuickBLAST(DllInfo *dll);

  QBLIBRARY_API SEXP libQB_GetInstanceCount();
  QBLIBRARY_API SEXP libQB_CreateQuickBLASTInstance(SEXP seq_type, SEXP strand, SEXP program, SEXP options, SEXP save_sequences);
  QBLIBRARY_API SEXP libQB_DeleteQuickBLASTInstance(SEXP ptr_id);
  //     //QBLIBRARY_API ArrowWrapper* CreateArrowWrapperInstance();
  //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
  //     //QBLIBRARY_API std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
  QBLIBRARY_API SEXP libQB_SetQuickBLASTOptions(SEXP ptr_id, SEXP program_name, SEXP options);
  QBLIBRARY_API SEXP BLAST2Seqs(SEXP ptr_id, SEXP query, SEXP subject);
  QBLIBRARY_API SEXP BLAST2Files(SEXP ptr_id, SEXP query, SEXP subject, SEXP out_file, SEXP seq_limit, SEXP num_threads, SEXP show_progress, SEXP return_values, SEXP min_batch_size);

  // __declspec(dllexport) SEXP test_QBR();
  // __declspec(dllexport) SEXP test_QBR_cpp();
  QBLIBRARY_API SEXP libQB_isQuickBLASTLoaded();

  /*SEXP GetQuickBLASTInstance(INTSXP id);
  INTSXP GetInstanceCount();
  INTSXP GetInstanceID(SEXP ptr);
  SEXP CreateQuickBLASTInstance(INTSXP seq_type, INTSXP strand, STRSXP program, STRSXP options, LGLSXP save_sequences);
  void DeleteQuickBLASTInstance(SEXP ptr);
  //ArrowWrapper* CreateArrowWrapperInstance();
  //std::shared_ptr<arrow::RecordBatch> cpp_BLAST2Seqs(std::shared_ptr<QuickBLAST> ptr, std::string query, std::string subject);
  //std::shared_ptr<arrow::RecordBatchVector> cpp_BLAST2Files(std::shared_ptr<QuickBLAST> ptr, std::string queryFile, std::string subjectFile, std::string outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
  void SetQuickBLASTOptions(SEXP ptr, STRSXP program_name, STRSXP options);
  SEXP BLAST2Seqs(SEXP ptr, STRSXP query, std::string subject);
  SEXP BLAST2Files(SEXP ptr, STRSXP queryFile, STRSXP subjectFile, STRSXP outFile, INTSXP blast_sequence_limit, INTSXP num_threads, LGLSXP show_progress, LGLSXP return_values, INTSXP batch_size = 1024);
  LGLSXP BLAST2Folders(SEXP ptr, STRSXP query, STRSXP subject, STRSXP extension, STRSXP out_folder, INTSXP num_threads, LGLSXP reciprocal_hits, INTSXP min_batch_size = 1024);
  LGLSXP BLAST1Folder(SEXP ptr, STRSXP input_folder, STRSXP extension, STRSXP out_folder, INTSXP num_threads, LGLSXP reciprocal_hits, INTSXP min_batch_size = 1024);*/
}

//' Stub function that always returns true. Only to test the connection of DLLs and C function calls
//' @examples
//' \dontrun{
//' QuickBLAST::isQuickBLASTLoaded()
//' }
//'
//' @return Always TRUE
//' @md
//' @export
// [[Rcpp::export]]
SEXP isQuickBLASTLoaded()
{
  std::string ret_str = Rcpp::as<std::string>(libQB_isQuickBLASTLoaded());
  Rcpp::Rcout << ret_str << std::endl;
  Rcpp::Rcout << "Rcpp - QuickBLAST dependencies Loaded!" << std::endl;
  return Rcpp::wrap(true);
}

//' @export
// [[Rcpp::export]]
SEXP GetInstanceCount()
{
  Rprintf("testing R side");
  int ret_val = Rcpp::as<int>(libQB_GetInstanceCount());
  Rprintf("testing R side");
  if (ret_val == obj_list.size())
  {
    Rcpp::Rcout << "Instance Count : " << ret_val << std::endl;
    return Rcpp::wrap(true);
  }
  else
  {
    Rcpp::Rcerr << "Instance Count : " << ret_val << " does not match with obj_list size : " << obj_list.size() << std::endl;
    return R_NilValue;
  }
}

// // [[Rcpp::export]]
// RcppExport SEXP test_QBR()
// {
//   Rcpp::Rcout << "Hello from QuickBLAST R side" << std::endl;
//   return Rcpp::wrap(true);
// }

// // [[Rcpp::export]]
// RcppExport SEXP test_QBR_cpp()
// {
//   Rcpp::Rcout << "Hello from QuickBLAST R side connected with cpp" << std::endl;
//   std::string ret_str = Rcpp::as<std::string>(test_QBcpp());
//   Rcpp::Rcout << ret_str << std::endl;
//   return Rcpp::wrap(true);
// }

// // [[Rcpp::export]]
// RcppExport SEXP test_Rnil()
// {
//   return R_NilValue;
// }

// bool QueryOOBESupport()
// {
//   return true;
// }

/*extern "C" attribute_visible void R_unload_QuickBLAST(DllInfo *dll)
{

}*/

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

      options_ += "-" + name + " " + Rcpp::as<std::string>(value) + " ";
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
  Rcpp::Rcout << options_ << std::endl;
  return options_;
}

//' @name CreateNewBLASTInstance
//'
//' @title Create new BLAST instance with seq_type, strand, program and BLAST options.
//'
//' @description Create a new QuickBLAST C++ object with seq_type, strand, program and BLAST options, which can be used in QuickBLAST::BLAST2Files() and QuickBLAST::BLAST2Seqs()
//'
//' @param seq_info (R List) list((int) "seq_type" : 0 - (eNucleotide) (OR) 1 - (eProtein), (int) "strand" = 0 - (ePlus) (OR) 1 - (eMinus), (bool) "save_sequences" = Save Sequences in the arrow::RecordBatch?). Named or unnamed List is accepted but in order.
//' @param program (string) Name of the BLAST program
//' @param options (string (or) Named List) List of BLAST options - check QuickBLAST::GetAvailableBLASTOptions()
//' @return C++ Pointer of QuickBLAST Object (Cannot be used in R)
//' @export
//' @useDynLib QuickBLAST, .registration = TRUE, .fixes="QB_"
// [[Rcpp::export]]
SEXP CreateNewBLASTInstance(SEXP seq_info, SEXP program, SEXP options)
{
  assert(TYPEOF(seq_info) == LISTSXP);
  assert(TYPEOF(program) == STRSXP);
  assert(TYPEOF(options) == LISTSXP || TYPEOF(options) == VECSXP || TYPEOF(options) == STRSXP);
  Rcpp::List seq_info_ = Rcpp::as<Rcpp::List>(seq_info);
  assert(seq_info_.size() == 3);
  int seq_type_ = Rcpp::as<int>(seq_info_[0]);
  int strand_ = Rcpp::as<int>(seq_info_[1]);
  bool save_sequences_ = Rcpp::as<bool>(seq_info_[2]);
  std::string options_ = ConvertBLASTOptions2String(options);
  std::string program_ = Rcpp::as<std::string>(program);

  SEXP seq_type_SEXP = Rcpp::wrap(seq_type_);
  SEXP strand_SEXP = Rcpp::wrap(strand_);
  SEXP options_SEXP = Rcpp::wrap(options_);
  SEXP save_sequences_SEXP = Rcpp::wrap(save_sequences_);

  SEXP obj_id = libQB_CreateQuickBLASTInstance(seq_type_SEXP, strand_SEXP, program, options_SEXP, save_sequences_SEXP);
  unsigned int obj_id_ = Rcpp::as<unsigned int>(obj_id);

  obj_list.insert(std::make_pair(obj_id_, std::make_tuple(obj_id_, seq_type_, strand_, program_, options_, save_sequences_)));
  // obj_list[obj_id_] = obj_tup_;
  return Rcpp::wrap(obj_list.size());
}

/*//' @name BLAST2Folders
//' @title BLAST 2 Folders
//'
//' @description BLAST files of query folder with files of subject folder.
//'
//' @note Output file is named as $query-$subject.hits in the out_folder.
//'
//' @param ptr (QuickBLAST ptr) Pointer to QuickBLAST object created with CreateNewBLASTInstance().
//' @param query (std::string) Query FASTA Folder
//' @param subject (std::string) Subject FASTA Folder
//' @param extension (std::string) Extension of FASTA files.
//' @param out_folder (string) Output Folder
//' @param num_threads (int) Number of Threads
//' @param reciprocal_hits (bool) BLAST bi-directionally?
//' @param min_batch_size (int) Minimum Batch Size to start writing to output file. (Default - 1024)
//' @return List of output filenames
// [[Rcpp::export]]
RcppExport SEXP QB_BLAST2Folders(SEXP ptr, SEXP query, SEXP subject, SEXP extension, SEXP out_folder, SEXP num_threads, SEXP reciprocal_hits, int min_batch_size = 1024)
{
  // int typ1 = TYPEOF(query);
  // int typ2 = TYPEOF(subject);

  // assert(typ1 == LISTSXP || typ1 == VECSXP);
  // assert(typ2 == LISTSXP || typ2 == VECSXP);
  assert(ptr != nullptr);
  int seq_limit = -1;
  auto start = std::chrono::high_resolution_clock::now();

  Rcpp::XPtr<QuickBLAST> qb_ptr = static_cast<Rcpp::XPtr<QuickBLAST>>(ptr);
  // convert inputs to appropriate C++ types
  // Rcpp::XPtr<QuickBLAST> ptr_(qb_ptr, false);
  std::shared_ptr<QuickBLAST> ptr_(qb_ptr.get());
  std::string query_ = as<std::string>(query);
  std::string subject_ = as<std::string>(subject);
  std::string extension_ = as<std::string>(extension);
  std::string outFolder_ = as<std::string>(out_folder);
  int num_threads_ = as<int>(num_threads);
  int min_batch_size_ = min_batch_size;
  bool reciprocal_hits_ = as<bool>(reciprocal_hits);
  std::filesystem::path outPath(outFolder_);
  std::filesystem::create_directory(outPath);
  if (!std::filesystem::is_empty(outPath))
  {
    Rcpp::Rcerr << "out_folder : Folder must be empty.";
    return Rcpp::wrap(false);
  }
  assert(min_batch_size_ > 0);
  assert(!query_.empty());
  assert(!subject_.empty());

  std::filesystem::path query_path = query_;
  std::filesystem::path subject_path = subject_;
  std::vector<std::string> queryFiles = getFilesInDir(query_, extension_);
  std::vector<std::string> subjectFiles = getFilesInDir(subject_, extension_);

  assert(!queryFiles.empty());
  assert(!subjectFiles.empty());

  int iterations = queryFiles.size() * subjectFiles.size();

  assert(iterations > 0);

  Rcpp::List ret_lst(iterations);
  Rcpp::CharacterVector names(iterations);
  // Progress progress_bar(iterations, true);

  for (int i = 0; i < (int)queryFiles.size(); i++)
  {
    for (int j = 0; j < (int)subjectFiles.size(); j++)
    {
      Rcpp::checkUserInterrupt();
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
        Rcpp::Rcerr << "Warn : File not found : " << query_input.string() << " or " << subject_input.string() << std::endl;
        continue;
      }

      if (!reciprocal_hits_)
      {
        std::string base_name_rbh = getFilenameWithoutExtension(subjectFiles[j]) + "-" + getFilenameWithoutExtension(queryFiles[i]);
        std::filesystem::path outFile_rbh = outPath / (base_name_rbh + ".hits");

        if (std::filesystem::exists(outFile_rbh.string()))
        {
          continue;
        }
      }
      int one_dim_index = i + queryFiles.size() * j;
      Rcpp::Rcout << "BLASTing : " << base_name << std::endl;
      QuickBLASTHandle in_val;
      in_val.ptr = ptr_;
      // static_cast<void>(ptr_->BLAST_files(query_input.string(), subject_input.string(), outFile_.string(), seq_limit, num_threads_, true, false, min_batch_size_));
      // std::shared_ptr<QuickBLAST> sh_ptr(ptr_.get());
      static_cast<void>(BLAST2Files(in_val, query_input.string(), subject_input.string(), outFile_.string(), seq_limit, num_threads_, true, false, min_batch_size_));
      ret_lst[one_dim_index] = outFile_.string();
      names[one_dim_index] = base_name;
    }
  }
  ret_lst.names() = names;
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  // Print the time in seconds
  Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;
  return rm_null(ret_lst);
}*/

/*//' @name BLAST1Folder
//' @title BLAST 1 Folder
//'
//' @description BLAST the contents of a folder.
//'
//' @note Output file is named as $query-$subject.hits in the out_folder.
//'
//' @param ptr (QuickBLAST ptr) Pointer to QuickBLAST object created with CreateNewBLASTInstance().
//' @param input_folder (std::string) Input Folder with FASTA (Query and Subjects)
//' @param extension (std::string) Extension of FASTA files.
//' @param out_folder (string) Output Folder
//' @param num_threads (int) Number of Threads
//' @param reciprocal_hits (bool) BLAST bi-directionally?
//' @param min_batch_size (int) Minimum Batch Size to start writing to output file. (Default - 1024)
//' @return List of output filenames
// [[Rcpp::export]]
RcppExport SEXP QB_BLAST1Folder(SEXP ptr, SEXP input_folder, SEXP extension, SEXP out_folder, SEXP num_threads, SEXP reciprocal_hits, int min_batch_size = 1024)
{

  // int typ1 = TYPEOF(input_folder);

  // assert(typ1 == LISTSXP || typ1 == VECSXP);

  assert(ptr != nullptr);
  int seq_limit = -1;
  auto start = std::chrono::high_resolution_clock::now();

  Rcpp::XPtr<QuickBLAST> qb_ptr = static_cast<Rcpp::XPtr<QuickBLAST>>(ptr);
  // convert inputs to appropriate C++ types
  // Rcpp::XPtr<QuickBLAST> ptr_(qb_ptr, false);
  std::shared_ptr<QuickBLAST> ptr_(qb_ptr.get());
  std::string input_folder_ = as<std::string>(input_folder);
  std::string extension_ = as<std::string>(extension);
  std::string outFolder_ = as<std::string>(out_folder);
  int num_threads_ = as<int>(num_threads);
  int min_batch_size_ = min_batch_size;
  bool reciprocal_hits_ = as<bool>(reciprocal_hits);
  std::filesystem::path outPath(outFolder_);
  std::filesystem::create_directory(outPath);
  if (!std::filesystem::is_empty(outPath))
  {
    Rcpp::Rcerr << "out_folder : Folder must be empty.";
    return Rcpp::wrap(false);
  }
  assert(min_batch_size_ > 0);
  assert(!input_folder_.empty());

  std::vector<std::string> inputFiles = getFilesInDir(input_folder_, extension_);

  assert(!inputFiles.empty());

  std::filesystem::path folder = input_folder_;

  int iterations = inputFiles.size() * inputFiles.size();

  assert(iterations > 0);

  Rcpp::List ret_lst(iterations);
  Rcpp::CharacterVector names(iterations);
  for (int i = 0; i < (int)inputFiles.size(); i++)
  {
    for (int j = 0; j < (int)inputFiles.size(); j++)
    {
      Rcpp::checkUserInterrupt();

      if (inputFiles[i] == inputFiles[j])
      {
        continue;
      }

      std::string base_name = getFilenameWithoutExtension(inputFiles[i]) + "-" + getFilenameWithoutExtension(inputFiles[j]);
      std::filesystem::path outFile_ = outPath / (base_name + ".hits");

      if (!reciprocal_hits_)
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
      Rcpp::Rcout << "BLASTing : " << base_name << std::endl;
      QuickBLASTHandle in_val;
      in_val.ptr = ptr_;
      // static_cast<void>(ptr_->BLAST_files(qry_filePath.string(), subj_filePath.string(), outFile_.string(), seq_limit, num_threads_, true, false, min_batch_size_));
      // std::shared_ptr<QuickBLAST> sh_ptr(ptr_.get());
      static_cast<void>(BLAST2Files(in_val, qry_filePath.string(), subj_filePath.string(), outFile_.string(), seq_limit, num_threads_, true, false, min_batch_size_));
      ret_lst[one_dim_index] = outFile_.string();
      names[one_dim_index] = base_name;
    }
  }
  ret_lst.names() = names;
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  // Print the time in seconds
  Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;

  return rm_null(ret_lst);
}*/

static const R_CallMethodDef callMethods[] = {
    {"isQuickBLASTLoaded", (DL_FUNC)&isQuickBLASTLoaded, 0},
    {"GetInstanceCount", (DL_FUNC)&GetInstanceCount, 0},
    {"BLAST2Files", (DL_FUNC)&BLAST2Files, 9},
    // {"BLAST1Folder", (DL_FUNC)&QB_BLAST1Folder, 7},
    // {"BLAST2Folders", (DL_FUNC)&QB_BLAST2Folders, 8},
    {"BLAST2Seqs", (DL_FUNC)&BLAST2Seqs, 3},
    {"CreateNewBLASTInstance", (DL_FUNC)&CreateNewBLASTInstance, 3},
    {NULL, NULL, 0}};

// static const R_CMethodDef cMethods[] = {
//     {"isQuickBLASTLoaded", (DL_FUNC)&QB_isQuickBLASTLoaded, 0},
//     {"GetInstanceCount", (DL_FUNC)&QB_GetInstanceCount, 0},
//     {"BLAST2Files", (DL_FUNC)&BLAST2Files, 9},
//     // {"BLAST1Folder", (DL_FUNC)&QB_BLAST1Folder, 7},
//     // {"BLAST2Folders", (DL_FUNC)&QB_BLAST2Folders, 8},
//     {"BLAST2Seqs", (DL_FUNC)&BLAST2Seqs, 3},
//     {"CreateNewBLASTInstance", (DL_FUNC)&QB_CreateNewBLASTInstance, 3},
//     {NULL, NULL, 0}};

/*
// R_CMethodDef
// R_ExternalMethodDef
static const R_CMethodDef cMethods[] = {
    {"test_QBR", (DL_FUNC)&test_QBR, 0},
    {"test_QBR_cpp", (DL_FUNC)&test_QBR_cpp, 0},
    //{"test_QBcpp", (DL_FUNC)&test_QBcpp, 0},
    // {"GetQuickBLASTInstance", (DL_FUNC)&GetQuickBLASTInstance, 1},
    // {"GetInstanceCount", (DL_FUNC)&GetInstanceCount, 0},
    // {"GetInstanceID", (DL_FUNC)&GetInstanceID, 1},
    // {"CreateQuickBLASTInstance", (DL_FUNC)&QB_CreateNewBLASTInstance, 5},
    // {"SetQuickBLASTOptions", (DL_FUNC)&SetQuickBLASTOptions, 3},
    // {"BLAST2Seqs", (DL_FUNC)&QB_BLAST2Seqs, 3},
    // {"BLAST2Files", (DL_FUNC)&QB_BLAST2Files, 9},
    // {"BLAST2Folders", (DL_FUNC)&QB_BLAST2Folders, 8},
    // {"BLAST1Folder", (DL_FUNC)&QB_BLAST1Folder, 7},
    NULL};*/

void R_init_QuickBLAST(DllInfo *dll)
{
  // R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);
}

/*RCPP_MODULE(blast_module)
{
  class_<QuickBLAST>("QuickBLAST")
      .constructor<QuickBLAST::ESeqType, QuickBLAST::EStrand, std::string, std::string, bool>()
      .constructor<QuickBLAST::ESeqType, QuickBLAST::EStrand, std::string, Rcpp::List, bool>()
      .method("BLAST", &QuickBLAST::BLAST);
}*/
