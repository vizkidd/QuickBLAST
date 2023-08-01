
// [[Rcpp::plugins(openmp)]]

#include <RcppCommon.h>
#include <Rcpp.h>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <functional>
#include <progress.hpp>
#include "quick_blast-functions.hpp"
#include "arrow_wrapper-functions.hpp"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(QuickBLAST);
RCPP_EXPOSED_CLASS(ArrowWrapper);

RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::ESeqType);
RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::EStrand);
RCPP_EXPOSED_ENUM_NODECL(QuickBLAST::EInputType);

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
// [[Rcpp::export]]
RcppExport SEXP CreateNewBLASTInstance(SEXP seq_info, SEXP program, SEXP options)
{
  // convert inputs to appropriate C++ types
  Rcpp::List seq_info_ = as<Rcpp::List>(seq_info);
  assert(seq_info_.size() == 3);
  QuickBLAST::ESeqType seq_type = static_cast<QuickBLAST::ESeqType>(as<int>(seq_info_[0]));
  QuickBLAST::EStrand strand = static_cast<QuickBLAST::EStrand>(as<int>(seq_info_[1]));
  bool save_sequences = as<bool>(seq_info_[2]);
  std::string program_ = as<std::string>(program);

  ncbi::blast::CBlastOptionsHandle *opts;
  int typ = TYPEOF(options);
  switch (typ)
  {
  case LISTSXP:
  case VECSXP:
  {
    Rcpp::List options_list = Rcpp::as<Rcpp::List>(options);
    QuickBLAST *ret_ptr = new QuickBLAST(seq_type, strand, program_, options_list, save_sequences);
    Rcpp::XPtr<QuickBLAST>
        ptr(ret_ptr, false);
    // return the external pointer to the R side
    return ptr;
    break;
  }
  case STRSXP:
  {

    std::string options_str = Rcpp::as<std::string>(options);
    QuickBLAST *ret_ptr = new QuickBLAST(seq_type, strand, program_, options_str, save_sequences);
    Rcpp::XPtr<QuickBLAST>
        ptr(ret_ptr, false);
    // return the external pointer to the R side
    return ptr;
    break;
  }
  default:
    Rcpp::Rcerr << "Only named list or string allowed for BLAST options : Check QuickBLAST::GetAvailableBLASTOptions()";
    return Rcpp::wrap(false);
    break;
  }
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

//' @name BLAST2Seqs
//' @title BLAST 2 Sequences Strings
//'
//' @description BLAST 2 Sequence Strings and return the RecordBatchVector of BLAST Hits
//'
//' @param ptr (QuickBLAST ptr) Pointer to QuickBLAST object created with CreateNewBLASTInstance().
//' @param query (string) Query Sequence
//' @param subject (string) Subject Sequence
//' @return RecordBatchVector of BLAST Hits
// [[Rcpp::export]]
RcppExport SEXP BLAST2Seqs(SEXP ptr, SEXP query, SEXP subject)
{
  auto start = std::chrono::high_resolution_clock::now();

  Rcpp::XPtr<QuickBLAST> ptr_(ptr);
  // convert inputs to appropriate C++ types
  std::string query_ = as<std::string>(query);
  std::string subject_ = as<std::string>(subject);

  assert(!query_.empty());
  assert(!subject_.empty());

  std::shared_ptr<arrow::RecordBatch> ret_val_ = ptr_->BLAST_seqs(query_, subject_);

  if (ret_val_)
  {
    arrow::Status rb_sts = ret_val_->ValidateFull();
    if (!rb_sts.ok())
    {
      Rcpp::Rcerr << "ERR : Invalid RB : " << rb_sts.message() << rb_sts.detail() << std::endl;
      return Rcpp::wrap(false);
    }
  }
  else
  {
    return Rcpp::wrap(false);
  }
  Rcpp::List ret_val = as<List>(ptr_->Hits2RList(ret_val_));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  // Print the time in seconds
  Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;

  return rm_null(ret_val);
}

std::vector<std::string> getFilesInDir(const std::string &folderPath, const std::string &extension)
{
  std::vector<std::string> outFiles;

  if (!extension.empty())
  {
    for (const auto &entry : filesystem::directory_iterator(folderPath))
    {
      if (entry.is_regular_file() && entry.path().extension() == extension)
      {
        outFiles.emplace_back(entry.path().filename().string());
      }
    }
  }
  else
  {
    for (const auto &entry : filesystem::directory_iterator(folderPath))
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

//' @name BLAST2Folders
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
//' @return List of output filenames
// [[Rcpp::export]]
RcppExport SEXP BLAST2Folders(SEXP ptr, SEXP query, SEXP subject, SEXP extension, SEXP out_folder, SEXP num_threads, SEXP reciprocal_hits)
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
  Rcpp::XPtr<QuickBLAST> ptr_(qb_ptr, false);
  std::string query_ = as<std::string>(query);
  std::string subject_ = as<std::string>(subject);
  std::string extension_ = as<std::string>(extension);
  std::string outFolder_ = as<std::string>(out_folder);
  int num_threads_ = as<int>(num_threads);
  bool reciprocal_hits_ = as<bool>(reciprocal_hits);
  std::filesystem::path outPath(outFolder_);
  std::filesystem::create_directory(outPath);
  if (!std::filesystem::is_empty(outPath))
  {
    Rcpp::Rcerr << "out_folder : Folder must be empty.";
    return Rcpp::wrap(false);
  }
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

  for (int i = 0; i < queryFiles.size(); i++)
  {
    for (int j = 0; j < subjectFiles.size(); j++)
    {
      Rcpp::checkUserInterrupt();
      // progress_bar.increment();
      if (queryFiles[i] == subjectFiles[j])
      {
        continue;
      }

      std::string base_name = getFilenameWithoutExtension(queryFiles[i]) + "-" + getFilenameWithoutExtension(subjectFiles[j]);
      std::filesystem::path outFile_ = outPath / base_name + ".hits";
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
        std::filesystem::path outFile_rbh = outPath / base_name_rbh + ".hits";

        if (std::filesystem::exists(outFile_rbh.string()))
        {
          continue;
        }
      }
      int one_dim_index = i + queryFiles.size() * j;
      Rcpp::Rcout << "BLASTing : " << base_name << std::endl;
      static_cast<void>(ptr_->BLAST_files(query_input.string(), subject_input.string(), outFile_.string(), seq_limit, num_threads_, true, false));
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
}

//' @name BLAST1Folder
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
//' @return List of output filenames
// [[Rcpp::export]]
RcppExport SEXP BLAST1Folder(SEXP ptr, SEXP input_folder, SEXP extension, SEXP out_folder, SEXP num_threads, SEXP reciprocal_hits)
{

  // int typ1 = TYPEOF(input_folder);

  // assert(typ1 == LISTSXP || typ1 == VECSXP);

  assert(ptr != nullptr);
  int seq_limit = -1;
  auto start = std::chrono::high_resolution_clock::now();

  Rcpp::XPtr<QuickBLAST> qb_ptr = static_cast<Rcpp::XPtr<QuickBLAST>>(ptr);
  // convert inputs to appropriate C++ types
  Rcpp::XPtr<QuickBLAST> ptr_(qb_ptr, false);
  std::string input_folder_ = as<std::string>(input_folder);
  std::string extension_ = as<std::string>(extension);
  std::string outFolder_ = as<std::string>(out_folder);
  int num_threads_ = as<int>(num_threads);
  bool reciprocal_hits_ = as<bool>(reciprocal_hits);
  std::filesystem::path outPath(outFolder_);
  std::filesystem::create_directory(outPath);
  if (!std::filesystem::is_empty(outPath))
  {
    Rcpp::Rcerr << "out_folder : Folder must be empty.";
    return Rcpp::wrap(false);
  }
  assert(!input_folder_.empty());

  std::vector<std::string> inputFiles = getFilesInDir(input_folder_, extension_);

  assert(!inputFiles.empty());

  std::filesystem::path folder = input_folder_;

  int iterations = inputFiles.size() * inputFiles.size();

  assert(iterations > 0);

  Rcpp::List ret_lst(iterations);
  Rcpp::CharacterVector names(iterations);
  for (int i = 0; i < inputFiles.size(); i++)
  {
    for (int j = 0; j < inputFiles.size(); j++)
    {
      Rcpp::checkUserInterrupt();

      if (inputFiles[i] == inputFiles[j])
      {
        continue;
      }

      std::string base_name = getFilenameWithoutExtension(inputFiles[i]) + "-" + getFilenameWithoutExtension(inputFiles[j]);
      std::filesystem::path outFile_ = outPath / base_name + ".hits";

      if (!reciprocal_hits_)
      {
        std::string base_name_rbh = getFilenameWithoutExtension(inputFiles[j]) + "-" + getFilenameWithoutExtension(inputFiles[i]);
        std::filesystem::path outFile_rbh = outPath / base_name_rbh + ".hits";

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
      static_cast<void>(ptr_->BLAST_files(qry_filePath.string(), subj_filePath.string(), outFile_.string(), seq_limit, num_threads_, true, false));
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
}

//' @name BLAST2Files
//' @title BLAST 2 Files
//'
//' @description BLAST 2 Files and return the RecordBatchVector of BLAST Hits Asynchronously
//'
//' @note Keep return_values_ = FALSE for large files to avoid choking R.
//'
//' @param ptr (QuickBLAST ptr) Pointer to QuickBLAST object created with CreateNewBLASTInstance().
//' @param query (string) Query FASTA File
//' @param subject (string) Subject FASTA File
//' @param out_file (string) Output Filename (Arrow Feather/IPC Format)
//' @param seq_limit (int) Batch Size to BLAST at a time. { -1 = Whole File, 0 - One sequence at a time or > 0 }
//' @param num_threads (int) Number of Threads
//' @param show_progress (bool) TRUE - Show progress, Set FALSE for multiple instances
//' @param return_values (bool) TRUE - Returns values as list, Default - FALSE - Does not return values (Return true on completion)
//' @return RecordBatchVector of Asynchronous BLAST Hits
// [[Rcpp::export]]
RcppExport SEXP BLAST2Files(SEXP ptr, SEXP query, SEXP subject, SEXP out_file, SEXP seq_limit, SEXP num_threads, SEXP show_progress, SEXP return_values)
{
  auto start = std::chrono::high_resolution_clock::now();
  Rcpp::XPtr<QuickBLAST> ptr_(ptr);
  // convert inputs to appropriate C++ types
  int seq_limit_ = as<int>(seq_limit);
  std::string query_ = as<std::string>(query);
  std::string subject_ = as<std::string>(subject);
  bool show_progress_ = as<bool>(show_progress);
  bool return_values_ = as<bool>(return_values);
  std::string outFile_ = as<std::string>(out_file);
  int num_threads_ = as<int>(num_threads);

  assert(!query_.empty());
  assert(!subject_.empty());
  assert(!outFile_.empty());

  if (return_values_)
  {
    std::shared_ptr<arrow::RecordBatchVector> ret_vals = ptr_->BLAST_files(query_, subject_, outFile_, seq_limit_, num_threads_, show_progress_, return_values_);
    Rcpp::List ret_vals_ = as<List>(ptr_->Hits2RList(*ret_vals));

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    // Print the time in seconds
    Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;

    return rm_null(ret_vals_);
  }
  else
  {
    static_cast<void>(ptr_->BLAST_files(query_, subject_, outFile_, seq_limit_, num_threads_, show_progress_, return_values_));
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    // Print the time in seconds
    Rcpp::Rcout << "Clock : " << elapsed_seconds.count() << " seconds" << std::endl;
    return Rcpp::wrap(true);
  }
}

RCPP_MODULE(blast_module)
{
  class_<QuickBLAST>("QuickBLAST")
      .constructor<QuickBLAST::ESeqType, QuickBLAST::EStrand, std::string, std::string, bool>()
      .method("BLAST", &QuickBLAST::BLAST);
}
