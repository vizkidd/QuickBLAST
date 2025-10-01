#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#include "omp.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <cassert>
#include <thread>
#include <filesystem>
#include <tuple>
#include <unistd.h>
#include <sys/mman.h>
#include <regex>
#include <mutex>
#include <cstdio> //tmpfile

#include <string>
#include <pwd.h>      // getpwuid, struct passwd
#include <sys/types.h>
#include <cstdlib>    // getenv
#include <limits.h> 

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

static std::string get_username_safe() {
  // 1) try getlogin_r (thread-safe)
#if defined(_POSIX_LOGIN_NAME_MAX)
  size_t bufsize = _POSIX_LOGIN_NAME_MAX + 1;
#else
  size_t bufsize = 256;
#endif
  std::string name;
  std::vector<char> buf(bufsize);
  if (getlogin_r(buf.data(), buf.size()) == 0) {
    if (buf[0] != '\0') {
      return std::string(buf.data());
    }
  }
  
  // 2) fallback to passwd entry for the effective UID
  struct passwd *pw = getpwuid(geteuid());
  if (pw && pw->pw_name) {
    return std::string(pw->pw_name);
  }
  
  // 3) fallback to environment variable
  const char* envu = std::getenv("USER");
  if (envu && envu[0] != '\0') {
    return std::string(envu);
  }
  
  // 4) last resort
  return std::string("unknown");
}

ArrowWrapper::Impl::Impl()
{
  try{
    // // Rprintf("DBG1 AW \n");
    // Rcpp::Rcout <<"DBG1 AW " << std::endl;
  
    // ok_promise.set_value(arrow::Status::OK());
  
    // outputStream = std::make_shared<std::ostringstream>();
    // outputStream = Rcpp::XPtr<std::ostringstream>();
  
    arrow_LFS = arrow::fs::LocalFileSystem();
    std::string username = "";
  #if defined(linux) //|| defined(MINGW32)
    username = get_username_safe(); //getlogin();
    // Rcpp::Rcout <<"here2" << username << std::endl;
  // #elif defined(WIN32)
  //   char username_arr[UNLEN + 1];
  //   DWORD username_len = UNLEN + 1;
  //   GetUserNameA(username_arr, &username_len);
  //   username = username_arr;
  #endif
  
    if (username.empty())
    {
      username = "Unknown";
    }
  
    username += "(QuickBLAST)";
    // // Rprintf("DBG2 AW \n");
    // Rcpp::Rcout <<"DBG2 AW " << std::endl;
    // props_bldr.compression(arrow::Compression::LZ4_FRAME);
    // props_bldr.created_by(username);
    // props_bldr.data_page_version(parquet::ParquetDataPageVersion::V2);
    // props_bldr.write_batch_size(1024);
    // props_bldr.encoding(parquet::Encoding::RLE);
    // props_bldr.version(parquet::ParquetVersion::PARQUET_2_LATEST);
    // props_bldr.enable_write_page_index();
  
    // arrow_props_bldr.set_engine_version(parquet::ArrowWriterProperties::EngineVersion::V2);
    // arrow_props_bldr.set_use_threads(true);
  
    // parquet_writer_props = props_bldr.build();
    // arrow_writer_props = arrow_props_bldr.build();
  
    ipc_options.allow_64bit = true; //false;
    ipc_options.use_threads = true;
    ipc_options.metadata_version = arrow::ipc::MetadataVersion::V5;
    ipc_options.codec = arrow::util::Codec::Create(arrow::Compression::LZ4_FRAME).ValueOrDie();
  
    ipc_options.write_legacy_ipc_format = false;
  
    csv_options.include_header = true; //false;
    csv_options.batch_size = 1024;
    csv_options.delimiter = '\t';
    csv_options.null_string = "<NULL>";
    csv_options.quoting_style = arrow::csv::QuotingStyle::Needed;
    
    ipc_options.write_legacy_ipc_format = false;
  
    parquet_props = WriterProperties::Builder()
      .max_row_group_length(64 * 1024)
      ->created_by(username)
      ->version(ParquetVersion::PARQUET_2_6)
      ->data_page_version(ParquetDataPageVersion::V2)
      ->compression(Compression::LZ4)
      ->build();
  
    parquet_arrow_props = ArrowWriterProperties::Builder().store_schema()
                          ->set_use_threads(true)
                          ->set_engine_version(ArrowWriterProperties::V2)
                          ->build();
  
    rbv_batch = std::make_shared<arrow::RecordBatchVector>();
    // rbv_batch = std::make_shared<std::deque<std::shared_ptr<arrow::RecordBatch>>>();
    blast_metadata = std::make_shared<arrow::KeyValueMetadata>();
    AddFASTAMetadata("format", "Arrow IPC");
    AddFASTAMetadata("Created By", username);
    fasta_schema = arrow::schema({arrow::field("index", arrow::int64()), arrow::field("header", arrow::utf8()), arrow::field("sequence", arrow::utf8())});
  
    seq_info_type = arrow::struct_({
        arrow::field("num_alignments", arrow::int64()),
        // arrow::field("seq_titles", arrow::struct_({arrow::field("qseq_title", arrow::utf8()),
                                        // arrow::field("sseq_title", arrow::utf8())})),
        arrow::field("seqids", arrow::struct_({arrow::field("qseqid", arrow::utf8()),
                                               arrow::field("sseqid", arrow::utf8())})),
        arrow::field("seqs", arrow::struct_({arrow::field("qseq", arrow::large_utf8()),
                                             arrow::field("sseq", arrow::large_utf8())})),
        arrow::field("strands", arrow::utf8()),
        arrow::field("lengths", arrow::struct_({arrow::field("qlen", arrow::int64()),
                                                arrow::field("slen", arrow::int64())})),
    });
    this->hsp_type = arrow::struct_({arrow::field("qhsp", arrow::large_utf8()),
                                     arrow::field("shsp", arrow::large_utf8()),
                                     arrow::field("pident", arrow::float64()),
                                     arrow::field("pident_gap", arrow::float64()),
                                     arrow::field("frames", arrow::utf8()),
                                     arrow::field("evalue", arrow::float64()),
                                     arrow::field("length", arrow::int64()),
                                     arrow::field("length01", arrow::float64()),
                                     arrow::field("qstart", arrow::int64()),
                                     arrow::field("qend", arrow::int64()),
                                     arrow::field("sstart", arrow::int64()),
                                     arrow::field("send", arrow::int64()),
                                     arrow::field("bitscore", arrow::float64()),
                                     arrow::field("score", arrow::float64()),
                                     arrow::field("qcovhsp", arrow::float64()),
                                     arrow::field("blast_score", arrow::float64()),
                                     arrow::field("gaps", arrow::int64()),
                                     arrow::field("nident", arrow::int64()),
                                     arrow::field("mismatch", arrow::int64()),
                                     arrow::field("positive", arrow::int64()),
                                     arrow::field("n_splices", arrow::int64()),
                                     arrow::field("hsp_num", arrow::int64()),
                                     arrow::field("sum_evalue", arrow::float64()),
                                     arrow::field("product_coverage", arrow::float64()),
                                     arrow::field("overall_identity", arrow::float64()),
                                     arrow::field("negative_count", arrow::int64()),
                                     arrow::field("matches", arrow::float64()),
                                     arrow::field("high_quality_percent_coverage", arrow::float64()),
                                     arrow::field("exon_identity", arrow::float64()),
                                     arrow::field("consensus_splices", arrow::float64()),
                                     arrow::field("comp_adj_method", arrow::float64())});
    this->alignment_scores_type = arrow::list({hsp_type});
  
    blast_schema = arrow::schema({arrow::field("seq_info", seq_info_type),
                                  arrow::field("hsps", hsp_type)});
    // // Rprintf("DBG3 AW \n");
    // Rcpp::Rcout <<"DBG3 AW " << std::endl;
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_init_lock(&rec_countLock);
    omp_init_lock(&writer_threadsLock);
    omp_init_lock(&rbv_batchLock);
    omp_init_lock(&rec_writerLock);
  #endif
    // // Rprintf("DBG4 AW \n");
    // Rcpp::Rcout <<"DBG4 AW " << std::endl;
    Rcpp::Rcout << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("ArrowWrapper::Impl(): Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("ArrowWrapper::Impl(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("ArrowWrapper::Impl(): C++ Runtime Error : ") + e.what());
  }
  catch(...){
    Rcpp::stop( "ArrowWrapper::Impl(): Unknown Exception" );
  }
}

ArrowWrapper::Impl::~Impl()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_destroy_lock(&rec_countLock);
  omp_destroy_lock(&writer_threadsLock);
  omp_destroy_lock(&rbv_batchLock);
  omp_destroy_lock(&rec_writerLock);
#endif

  // std::cout << "~ArrowWrapper " << std::endl;
  // Rprintf("~ArrowWrapper::Impl");
}

ArrowWrapper::ArrowWrapper()
    : pImpl(std::make_unique<Impl>()) {}

// Destructor: default implementation (pImpl will automatically clean up)
ArrowWrapper::~ArrowWrapper() = default;

void ArrowWrapper::Impl::SetBatchSize(unsigned int batch_size)
{
  // assert(batch_size > 0);
  this->rb_batch_size.store(batch_size);
  if(batch_size > 0)
  {
    csv_options.batch_size = batch_size;
    parquet_batch_size = batch_size;
  }
  Rcpp::Rcout << "Batch size: " << batch_size << std::endl << std::flush; //DEBUG
}
void ArrowWrapper::SetBatchSize(unsigned int batch_size)
{
  pImpl->SetBatchSize(batch_size);
}

unsigned int ArrowWrapper::Impl::GetBatchSize()
{
  return this->rb_batch_size.load();
}
unsigned int ArrowWrapper::GetBatchSize()
{
  return pImpl->GetBatchSize();
}

int ArrowWrapper::Impl::GetRecordCount() { return rec_count; }
int ArrowWrapper::GetRecordCount() { return pImpl->GetRecordCount(); }

void ArrowWrapper::Impl::ResetRecordCount()
{
  rec_count = 1;
}
void ArrowWrapper::ResetRecordCount()
{
  pImpl->ResetRecordCount();
}

void ArrowWrapper::Impl::AddRecordCount()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_set_lock(&rec_countLock);
#endif
  rec_count++;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_unset_lock(&rec_countLock);
#endif
}
void ArrowWrapper::AddRecordCount()
{
  pImpl->AddRecordCount();
}

void ArrowWrapper::Impl::SetThreadCount(int num_threads)
{
  if (num_threads > 1)
  {
    ipc_options.use_threads = true;
    // #if defined(_OPENMP) || defined(WIN32)
    //     omp_set_num_threads(num_threads);
    // #endif
  }
  else
  {
    ipc_options.use_threads = false;
    // #if defined(_OPENMP) || defined(WIN32)
    //     omp_set_num_threads(1);
    // #endif
  }
  n_threads = num_threads;
}
void ArrowWrapper::SetThreadCount(int num_threads)
{
  pImpl->SetThreadCount(num_threads);
}

// std::string ArrowWrapper::Impl::CastToType(const std::string &full_entry)
// {
//   return full_entry;
// }
// std::string ArrowWrapper::CastToType(const std::string &full_entry)
// {
//   return full_entry;
// }

// FastaSequenceData ArrowWrapper::Impl::CastToType(const std::string &full_entry)
// {
//   // std::regex pattern("^>(\\w+)[\\r\\n]+([\\w\\W]+)$");
//   // std::regex punct_pattern("[[:punct:][:space:]\\n\\t]");
//   // std::smatch match;
//   // 
//   // std::string full_entry_str(static_cast<std::string>(full_entry));
//   // FastaSequenceData fasta_data;
//   // fasta_data.rec_no = GetRecordCount();
//   // if (std::regex_match(full_entry, match, pattern))
//   // {
//   //   fasta_data.header = match[1].str();
//   //   fasta_data.seq = std::regex_replace(match[2].str(), punct_pattern, "");
//   // }
//   // else
//   // {
//   // 
//   //   fasta_data.seq = std::regex_replace(full_entry, punct_pattern, "");
//   // }
//   // fasta_data.header = fasta_data.header.empty() ? std::to_string(fasta_data.rec_no) : fasta_data.header;
//   
//   static const std::regex pattern(R"(^>([^\r\n]+)\r?\n([\s\S]+)$)");
//   static const std::regex seq_clean(R"([^A-Za-z])"); // remove anything that's not a letter
//   
//   FastaSequenceData fasta_data;
//   fasta_data.rec_no = GetRecordCount();
//   
//   std::smatch match;
//   if (std::regex_match(full_entry, match, pattern) && full_entry[0] == '>') {
//     fasta_data.header = match[1].str();         // whole header line (keeps pipes/spaces)
//     std::string raw_seq = match[2].str();      // sequence possibly with newlines/spaces
//     fasta_data.seq = std::regex_replace(raw_seq, seq_clean, "");
//   } else {
//     // If pattern fails, fallback: take only letters from the whole input (conservative)
//     fasta_data.seq = std::regex_replace(full_entry, seq_clean, "");
//   }
//   
//   if (fasta_data.header.empty()) {
//     fasta_data.header = std::to_string(fasta_data.rec_no);
//   }
//   
//   AddRecordCount();
//   
//   // Rcpp::Rcout << ">" << fasta_data.header << std::endl << std::flush; //DEBUG
//   // Rcpp::Rcout << fasta_data.seq << std::endl << std::flush; //DEBUG
//   
//   return fasta_data;
// }
// FastaSequenceData ArrowWrapper::CastToType(const std::string &full_entry)
// {
//   return pImpl->CastToType(full_entry);
// }

FastaSequenceData ArrowWrapper::Impl::CastToType(const std::string_view &full_entry_sv)
{
  // full_entry_sv points inside the mmap buffer; DO NOT store
  // string_view beyond the lifetime of that buffer.
  FastaSequenceData fasta_data;
  fasta_data.rec_no = GetRecordCount(); // you already increment record count elsewhere
  
  // Expect format:
  // >header-line\r?\n
  // sequence lines (maybe many, with whitespace)
  // We avoid regex; do simple manual parsing.
  
  // Rcpp::Rcout << full_entry_sv << std::endl << std::flush; //DEBUG
  
  if (full_entry_sv.empty()) {
    fasta_data.header = std::to_string(fasta_data.rec_no);
    fasta_data.seq.clear();
    return fasta_data;
  }
  
  // If entry starts with '>' we parse header; otherwise treat entire input as sequence.
  size_t pos = 0;
  if (full_entry_sv[0] == '>') {
    // find first newline
    size_t newline_pos = full_entry_sv.find_first_of("\r\n");
    if (newline_pos == std::string_view::npos) {
      // header only
      std::string header(full_entry_sv);
      fasta_data.header = header;
      fasta_data.seq.clear();
      return fasta_data;
    }
    // header is between 1 and newline_pos
    std::string_view header_sv = full_entry_sv.substr(1, newline_pos - 1);
    // trim trailing spaces (cheap)
    auto trim_end = header_sv.size();
    while (trim_end > 0 && std::isspace(static_cast<unsigned char>(header_sv[trim_end-1]))) --trim_end;
    fasta_data.header.assign(header_sv.data(), trim_end);
    
    // the sequence region starts after the newline(s)
    size_t seq_start = newline_pos;
    // skip possible \r or \n sequences
    while (seq_start < full_entry_sv.size() && (full_entry_sv[seq_start] == '\r' || full_entry_sv[seq_start] == '\n')) ++seq_start;
    pos = seq_start;
  } else {
    // no header: entire entry is sequence
    fasta_data.header = std::to_string(fasta_data.rec_no);
    pos = 0;
  }
  
  // Now parse sequence region from pos to end; keep only letters.
  size_t remaining = full_entry_sv.size() - pos;
  if (remaining == 0) {
    fasta_data.seq.clear();
    return fasta_data;
  }
  
  // Reserve roughly remaining size to avoid repeated reallocation
  std::string seq;
  seq.reserve(remaining);
  
  const char* data = full_entry_sv.data();
  for (size_t i = pos; i < full_entry_sv.size(); ++i) {
    unsigned char ch = static_cast<unsigned char>(data[i]);
    // ASCII letters only; faster than locale-aware
    if ((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) {
      seq.push_back(static_cast<char>(ch));
    }
    // else skip digits/punctuation/whitespace
  }
  
  // shrink_to_fit if you want to reduce memory after large sequences (optional)
  // seq.shrink_to_fit();
  
  
  
  fasta_data.seq.swap(seq); // move efficiently
  
  // Rcpp::Rcout << fasta_data.rec_no << ">" << fasta_data.header << std::endl << std::flush; //DEBUG
  // Rcpp::Rcout << fasta_data.seq << std::endl << std::flush; //DEBUG
  
  return fasta_data;
}
FastaSequenceData ArrowWrapper::CastToType(const std::string_view &full_entry_sv)
{
  return pImpl->CastToType(full_entry_sv);
}

ArrowWrapper::EOutputFormat ArrowWrapper::Impl::OutputFormat2Enum(const std::string& str){
  if (str == "ipc") return ArrowWrapper::EOutputFormat::eIPC;
  if (str == "csv")  return ArrowWrapper::EOutputFormat::eCSV;
  if (str == "parquet")  return ArrowWrapper::EOutputFormat::eParquet;
  return ArrowWrapper::EOutputFormat::unknown;
}

// arrow::Result<std::shared_ptr<arrow::RecordBatch>> ArrowWrapper::Impl::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
arrow::Status ArrowWrapper::Impl::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
{
  try{
    if(!save2file)
    {
      return arrow::Status::OK(); //arrow::Result<std::shared_ptr<arrow::RecordBatch>>(rb_);
    }
    
    arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
    if (rb_)
    {
 
     if(rb_->num_rows() == 0)
       return arrow::Status::OK();
     
  //     //Wait until current batch is written
  //     {
  //       std::unique_lock<std::mutex> lk(rbv_mutex);
  //       rbv_not_full.wait(lk, [this]() {
  //         Rcpp::checkUserInterrupt();
  // #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  //         omp_set_lock(&rbv_batchLock);
  // #endif
  //         auto rbv_batch_size = rbv_batch->size();
  // #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  //         omp_unset_lock(&rbv_batchLock);
  // #endif
  //         return !writer_writing.load() && (rbv_batch_size < rb_batch_size || rb_batch_size == 0);
  //       });
  //       lk.unlock();
  //     }

  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_set_lock(&rbv_batchLock);
  #endif
      rbv_batch->emplace_back(std::move(rb_));
      unsigned int ret_size = rbv_batch->size();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_unset_lock(&rbv_batchLock);
  #endif

      if(writer_writing.load() || (ret_size < rb_batch_size.load() || rb_batch_size.load() == 0)){
        return arrow::Status::OK(); //arrow::Result<std::shared_ptr<arrow::RecordBatch>>(rb_);
      }


      {
        std::thread write_thread([this]()
        {
            try{
              // if(this->writer_threads.size() > 1){
              //   if(writer_threads.front().joinable())
              //     writer_threads.front().join();
              //   writer_threads.erase(writer_threads.begin());
              // }
              auto tmp_batch_size = rb_batch_size.load();
              // if(writer_threads.size() <= 1){
              if(writer_threads.empty()){
                tmp_batch_size++;
              }
              else{
                if(tmp_batch_size > 0)
                  tmp_batch_size--;
              }
              rb_batch_size.store(tmp_batch_size);
              // this->SetBatchSize(tmp_batch_size); // causes stack-smashing?
              static_cast<void>(this->WriteBatch2File());
              writer_writing.store(false);
              // rbv_not_full.notify_one();
              // rbv_not_full.notify_all();
            }catch(const Rcpp::exception &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): Rcpp Exception : ") + e.what());
            }
            catch(const std::exception &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): C++ Exception : ") + e.what());
            }
            catch(const std::runtime_error &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): C++ Runtime Exception : ") + e.what());
            }
            catch(...){
              Rcpp::stop("thread::WriteBatch2File(): Unknown Exception");
            }
          });
        // write_thread.detach(); // when detached, threads cannot join
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_set_lock(&writer_threadsLock);
#endif
        static_cast<void>(writer_threads.emplace_back(std::move(write_thread)));
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&writer_threadsLock);
#endif
      }

      // Rcpp::Rcout << "AddRB2Batch(): Writer threads : " << writer_threads.size() << " : " << ret_size << std::endl << std::flush; //DEBUG
  
      return arrow::Status::OK();
    }
    return error_sts; //arrow::Result<std::shared_ptr<arrow::RecordBatch>>(error_sts);
  }catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("AddRB2Batch(): Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("AddRB2Batch(): C++ Exception : ") + e.what() );
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("AddRB2Batch(): C++ Runtime Exception : ") + e.what() );
  }
  catch(...){
    Rcpp::stop("AddRB2Batch(): Unknown Exception");
  }
}
// arrow::Result<std::shared_ptr<arrow::RecordBatch>> ArrowWrapper::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
arrow::Status ArrowWrapper::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
{
  return pImpl->AddRB2Batch(rb_);
}

// arrow::Result<std::shared_ptr<arrow::RecordBatchVector>> ArrowWrapper::Impl::AddRBV2Batch(arrow::RecordBatchVector &rbv_)
arrow::Status ArrowWrapper::Impl::AddRBV2Batch(arrow::RecordBatchVector &rbv_)
{
  try{
    if(!save2file)
    {
      return arrow::Status::OK(); //arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(std::make_shared<arrow::RecordBatchVector>(rbv_));
    }
    
    arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
    if (rbv_.size() > 0)
    {

  //     {
  //       std::unique_lock<std::mutex> lk(rbv_mutex);
  //       rbv_not_full.wait(lk, [this]() {
  //         Rcpp::checkUserInterrupt();
  // #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  //         omp_set_lock(&rbv_batchLock);
  // #endif
  //         auto rbv_batch_size = rbv_batch->size();
  // #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  //         omp_unset_lock(&rbv_batchLock);
  // #endif
  //         return (!writer_writing.load() && (rbv_batch_size < rb_batch_size || rb_batch_size == 0)); 
  //       });
  //       lk.unlock();
  //     }
  
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_set_lock(&rbv_batchLock);
  #endif
      rbv_batch->insert(rbv_batch->end(), rbv_.begin(), rbv_.end());
      unsigned int ret_size = rbv_batch->size();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_unset_lock(&rbv_batchLock);
  #endif

      rbv_.clear();
      rbv_.shrink_to_fit();

      if(writer_writing.load() || (ret_size < rb_batch_size.load() || rb_batch_size.load() == 0)){
        return arrow::Status::OK(); //arrow::Result<std::shared_ptr<arrow::RecordBatch>>(rb_);
      }


      {
        std::thread write_thread([this]()
        {
            try{
              // if(this->writer_threads.size() > 1){
              //     if(writer_threads.front().joinable())
              //       writer_threads.front().join();
              //     writer_threads.erase(writer_threads.begin());
              // }
              auto tmp_batch_size = rb_batch_size.load();
              // if(writer_threads.size() <= 1){
              if(writer_threads.empty()){
                tmp_batch_size++;
              }
              else{
                if(tmp_batch_size > 0)
                  tmp_batch_size--;
              }
              rb_batch_size.store(tmp_batch_size);
              // this->SetBatchSize(tmp_batch_size); // causes stack-smashing?
              static_cast<void>(this->WriteBatch2File());
              writer_writing.store(false);
              // rbv_not_full.notify_one();
              // rbv_not_full.notify_all();
            }catch(const Rcpp::exception &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): Rcpp Exception : ") + e.what());
            }
            catch(const std::exception &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): C++ Exception : ") + e.what());
            }
            catch(const std::runtime_error &e){
              Rcpp::stop(std::string("thread::WriteBatch2File(): C++ Runtime Exception : ") + e.what());
            }
            catch(...){
              Rcpp::stop("thread::WriteBatch2File(): Unknown Exception");
            }
          });
        // write_thread.detach();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_set_lock(&writer_threadsLock);
#endif
        static_cast<void>(writer_threads.emplace_back(std::move(write_thread)));
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&writer_threadsLock);
#endif
      }

      // Rcpp::Rcout << "AddRBV2Batch(): Writer threads : " << writer_threads.size() << std::endl << std::flush; //DEBUG
      
      return arrow::Status::OK();
    }
    else{
      return arrow::Status::OK();
    }
    //}
    return error_sts; //arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(error_sts);
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("AddRBV2Batch(): Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("AddRBV2Batch(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("AddRBV2Batch(): C++ Runtime Exception : ") + e.what());
  }
  catch(...){
    Rcpp::stop("AddRBV2Batch(): Unknown Exception");
  }
}
// arrow::Result<std::shared_ptr<arrow::RecordBatchVector>> ArrowWrapper::AddRBV2Batch(arrow::RecordBatchVector &rbv_)
arrow::Status ArrowWrapper::AddRBV2Batch(arrow::RecordBatchVector &rbv_)
{
  return pImpl->AddRBV2Batch(rbv_);
}

std::string ArrowWrapper::Impl::GetOutputFormat(){
  return output_format;
}
std::string ArrowWrapper::GetOutputFormat(){
  return pImpl->GetOutputFormat();
}

arrow::Status ArrowWrapper::Impl::CreateOutputStream(std::string &outFile, const std::string& outputFormat)
{
  try{
    output_format = outputFormat;
    
    if(outFile.empty()){
      // outFile = std::tmpnam(nullptr); 
      save2file = false;
    }else{
      save2file = true;
      output_filename = outFile;
      Rcpp::Rcout << "Writing to : " << output_filename << std::endl << std::flush; //DEBUG
      Rcpp::Rcout << "Output Format : " << output_format << std::endl << std::flush; //DEBUG
    
      // // writer_loop() CODE
      writer_running.store(true);
      writer_writing.store(false);
      writer_finishing.store(false);
      writer_failed.store(false);
      // writer_thread = std::thread([this](){ 
      //     rbv_batch->clear();
      //     rbv_batch->shrink_to_fit();
      //     this->writer_loop(); 
      //   });
      
      rbv_batch->clear();
      rbv_batch->shrink_to_fit();
      
      switch(OutputFormat2Enum(output_format)){
      case ArrowWrapper::EOutputFormat::eIPC: {
          auto outFileStream_res = arrow_LFS.OpenAppendStream(output_filename, blast_metadata);
          if(!outFileStream_res.ok()){
            throw std::runtime_error(std::string("CreateOutputStream() - Could not open append stream to ") + output_filename);
          }
          outFileStream = outFileStream_res.ValueOrDie();
          auto writer_ = arrow::ipc::MakeFileWriter(outFileStream.get(), GetBLASTSchema(), GetArrowIPCOptions(), GetBLASTMetadata());
          if (!writer_.ok())
          {
            throw std::runtime_error(std::string("WriteBatch2File() - Error initiating IPC file writer: ") + writer_.status().message());
            return writer_.status();
          }
          rec_writer = writer_.ValueOrDie();
          break;
        }
      case ArrowWrapper::EOutputFormat::unknown:
      default: {
        throw std::runtime_error(std::string("CreateOutputStream() - Unsupported output format. Supported values ipc/csv/parquet"));
      }
      case ArrowWrapper::EOutputFormat::eCSV: {
          auto outFileStream_res = arrow_LFS.OpenAppendStream(output_filename, blast_metadata);
          if(!outFileStream_res.ok()){
            throw std::runtime_error(std::string("CreateOutputStream() - Could not open append stream to ") + output_filename);
          }
          outFileStream = outFileStream_res.ValueOrDie();
          auto writer_ = arrow::csv::MakeCSVWriter(outFileStream.get(), GetBLASTSchema(), GetArrowCSVOptions());
          if (!writer_.ok())
          {
            throw std::runtime_error(std::string("CreateOutputStream() - Error initiating CSV file writer: ") + writer_.status().message());
            return writer_.status();
          }
          rec_writer = writer_.ValueOrDie();
          break;
        }
      case ArrowWrapper::EOutputFormat::eParquet: {
          ARROW_ASSIGN_OR_RAISE(parquetFileStream,  arrow::io::FileOutputStream::Open(output_filename, /*append=*/ true));
          break;
        }
      }

      // finisher_running.store(true);
      finisher_thread = std::thread([this]() {
        try {
          while ( writer_running.load() || !writer_threads.empty() ) {
            // assert(!Progress::check_abort()); // R API calls inside threads will crash c++ runtime
            // Rcpp::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
            std::thread thread_to_join;
            {
              // scoped lock for safe access to writer_threads
              std::unique_lock<std::mutex> lk(writer_threads_mutex);
              if (writer_threads.empty()) { //writer_threads.size() >= 1
                // nothing to join right now — wait a bit (or better: use condition_variable)
                lk.unlock();
                // std::this_thread::sleep_for(std::chrono::milliseconds(50));
                continue;
              }
              // move out the first thread into local variable while holding lock
              thread_to_join = std::move(writer_threads.front());
              // erase first element; vector remains consistent
              writer_threads.erase(writer_threads.begin());
            } // unlock before join
            
            // join outside lock (can block)
            if (thread_to_join.joinable()) {
              thread_to_join.join();
            }
          }
          // finishing_cond.notify_all();
          writer_finishing.store(true);
        } catch(const Rcpp::exception &e){
          throw std::runtime_error(std::string("finisher_thread(): Rcpp Exception : ") + e.what());
        }
        catch(const std::exception &e){
          throw std::runtime_error(std::string("finisher_thread(): C++ Exception : ") + e.what());
        }
        catch(const std::runtime_error &e){
          throw std::runtime_error(std::string("finisher_thread(): C++ Runtime Error : ") + e.what());
        }
        catch(...){
          throw std::runtime_error("finisher_thread(): Unknown Exception");
        }
        // finisher_running.store(false);
      }); //std::move()
      finisher_thread.detach();   
//       finisher_thread = std::thread([this](){ 
//         try{
//           while(writer_running.load() || !writer_threads.empty()){
//             if(!writer_threads.empty()) //.size() crashes if vector is empty
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//               omp_set_lock(&writer_threadsLock);
// #endif
//             auto num_writers = writer_threads.size();
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//             omp_unset_lock(&writer_threadsLock);
// #endif
//             if(num_writers > 0){
//               if(writer_threads.front().joinable())
//                 writer_threads.front().join();
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//               omp_set_lock(&writer_threadsLock);
// #endif
//               writer_threads.erase(writer_threads.begin());
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//               omp_unset_lock(&writer_threadsLock);
// #endif
//             }
//             
//           }
//           // if(writer_threads.front().joinable())
//           //   writer_threads.front().join();
//         }catch(const Rcpp::exception &e){
//           Rcpp::stop(std::string("finisher_thread(): Rcpp Exception : ") + e.what());
//         }
//         catch(const std::exception &e){
//           Rcpp::stop(std::string("finisher_thread(): C++ Exception : ") + e.what());
//         }
//         catch(const std::runtime_error &e){
//           Rcpp::stop(std::string("finisher_thread(): C++ Runtime Error : ") + e.what());
//         }
//         catch(...){
//           Rcpp::stop("finisher_thread(): Unknown Exception");
//         }
//       });
      
    }
  
  }catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("CreateOutputStream(): Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("CreateOutputStream(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("CreateOutputStream(): C++ Runtime Error : ") + e.what());
  }
  catch(...){
    Rcpp::stop("CreateOutputStream(): Unknown Exception");
  }
  return arrow::Status::OK();
}

arrow::Status ArrowWrapper::CreateOutputStream(std::string &outFile, const std::string& outputFormat)
{
  return pImpl->CreateOutputStream(outFile, outputFormat);
}

std::shared_ptr<arrow::DataType> ArrowWrapper::Impl::GetSeqInfoType(void)
{
  return seq_info_type;
}
std::shared_ptr<arrow::DataType> ArrowWrapper::GetSeqInfoType(void)
{
  return pImpl->GetSeqInfoType();
}

std::shared_ptr<arrow::DataType> ArrowWrapper::Impl::GetAlignmentScoresType(void)
{
  return alignment_scores_type;
}
std::shared_ptr<arrow::DataType> ArrowWrapper::GetAlignmentScoresType(void)
{
  return pImpl->GetAlignmentScoresType();
}

std::shared_ptr<arrow::DataType> ArrowWrapper::Impl::GetHSPType(void)
{
  return hsp_type;
}
std::shared_ptr<arrow::DataType> ArrowWrapper::GetHSPType(void)
{
  return pImpl->GetHSPType();
}

std::shared_ptr<arrow::Schema> ArrowWrapper::Impl::GetBLASTSchema(void)
{
  return blast_schema;
}
std::shared_ptr<arrow::Schema> ArrowWrapper::GetBLASTSchema(void)
{
  return pImpl->GetBLASTSchema();
}

std::shared_ptr<arrow::Schema> ArrowWrapper::Impl::GetFASTASchema(void)
{
  return fasta_schema;
}
std::shared_ptr<arrow::Schema> ArrowWrapper::GetFASTASchema(void)
{
  return pImpl->GetFASTASchema();
}

std::shared_ptr<arrow::Schema> ArrowWrapper::Impl::GetSchema(void)
{
  if (!blast_schema) {
    // initialize default options or throw helpful error
    Rcpp::stop("ArrowWrapper::GetSchema(): blast_schema is not initialised");
  }
  return blast_schema;
}

std::shared_ptr<arrow::Schema> ArrowWrapper::GetSchema(void)
{
  return pImpl->GetSchema();
}

/*std::shared_ptr<parquet::WriterProperties> ArrowWrapper::GetParquetWriterProps(void)
{
    return parquet_writer_props;
}*/
/*std::shared_ptr<parquet::ArrowWriterProperties> ArrowWrapper::GetArrowWriterProps(void)
{
    return arrow_writer_props;
}*/
std::shared_ptr<arrow::KeyValueMetadata> ArrowWrapper::Impl::GetBLASTMetadata(void)
{
  return blast_metadata;
}
std::shared_ptr<arrow::KeyValueMetadata> ArrowWrapper::GetBLASTMetadata(void)
{
  return pImpl->GetBLASTMetadata();
}

void ArrowWrapper::Impl::AddFASTAMetadata(const std::string &key, const std::string &value)
{
  blast_metadata->Append(key, value);
}
void ArrowWrapper::AddFASTAMetadata(const std::string &key, const std::string &value)
{
  pImpl->AddFASTAMetadata(key, value);
}

arrow::ipc::IpcWriteOptions ArrowWrapper::Impl::GetArrowIPCOptions(void)
{
  return ipc_options;
}
arrow::ipc::IpcWriteOptions ArrowWrapper::GetArrowIPCOptions(void)
{
  return pImpl->GetArrowIPCOptions();
}

arrow::csv::WriteOptions ArrowWrapper::Impl::GetArrowCSVOptions(void)
{
  return csv_options;
}
arrow::csv::WriteOptions ArrowWrapper::GetArrowCSVOptions(void)
{
  return pImpl->GetArrowCSVOptions();
}

std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> ArrowWrapper::Impl::MMapFile(const std::string_view &filename, const char *delim)
{
  FILE *file_ptr;
  long fileSize;
  char *end_of_file_ptr;
  std::string delim_str(delim);

  file_ptr = fopen(filename.data(), "r");

  if (!file_ptr)
  {
    // std::cerr << "Error: Failed to open file: " << filename.data() << std::endl;
    REprintf("Error: Failed to open file: %s \n", filename.data());
    return nullptr;
  }

  // Get the file size
  fileSize = GetFileSize(file_ptr);
#if defined(linux) || defined(MINGW32)
  char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file_ptr), 0)); // MAP_SHARED
#else
  char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, _fileno(file_ptr), 0)); // MAP_SHARED
#endif // linux

  // char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file_ptr), 0)); // MAP_SHARED

  if (fileData_ptr == MAP_FAILED)
  {
    // std::cerr << "Error: Failed to map file " << filename.data() << std::endl;
    REprintf("Error: Failed to map file : %s \n", filename.data());
    fclose(file_ptr);
    return nullptr;
  }

  end_of_file_ptr = fileData_ptr + fileSize;

  return std::make_shared<std::tuple<FILE *, std::shared_ptr<char>, long, char *>>(std::make_tuple(file_ptr, std::shared_ptr<char>(fileData_ptr, [fileSize, file_ptr](char *ptr)
                                                                                                                                   {
    munmap(ptr, fileSize);
    fclose(file_ptr); }),
                                                                                                   fileSize, end_of_file_ptr));
}
// std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> ArrowWrapper::MMapFile(const std::string_view &filename, const char *delim)
// {
//   return pImpl->MMapFile(filename, delim);
// }

long ArrowWrapper::Impl::GetFileSize(FILE *file_ptr)
{
  fseek(file_ptr, 0, SEEK_END);
  long fileSize = ftell(file_ptr);
  rewind(file_ptr);
  return fileSize;
}

arrow::Status ArrowWrapper::Impl::FinishOutputStream()
{
  try{
    
    // // writer_loop() CODE
    // if (writer_thread.joinable()) {
    //   writer_thread.join();
    // }
    
    // for (auto &t : writer_threads) {
    //   if (t.joinable())
    //     t.join();
    // }
    
    if(!save2file){
      return arrow::Status::OK(); 
    }
    
    writer_running.store(false);

    // finisher_thread.join(); //cannot join detached threads

    {
      std::unique_lock<std::mutex> lk(finishing_mutex);
      finishing_cond.wait(lk, [this]() {
        Rcpp::checkUserInterrupt();
        return writer_finishing.load();
      });
      lk.unlock();
    }
        
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_set_lock(&rbv_batchLock);
  #endif
    auto rbv_size = rbv_batch->size();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_unset_lock(&rbv_batchLock);
#endif

    if (rbv_size > 0)
    {
      ARROW_RETURN_NOT_OK(WriteBatch2File());
      // signal writer to stop
      // writer_running.store(false);
    }
    
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
   omp_set_lock(&rbv_batchLock);
  #endif
    rbv_batch->clear();
    rbv_batch->shrink_to_fit();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_unset_lock(&rbv_batchLock);
  #endif
    Rcpp::Rcout << "Done writing to file." << std::endl <<std::flush; //DEBUG

    if(outFileStream)
      ARROW_RETURN_NOT_OK(outFileStream->Flush());
    
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_set_lock(&rec_writerLock);
#endif
    if(rec_writer){
      ARROW_RETURN_NOT_OK(rec_writer->Close());
      rec_writer.reset();
    }

    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_unset_lock(&rec_writerLock);
#endif
    
    if(parquetFileStream){
      if(!parquetFileStream->closed())
      {
        arrow::Status st2 = parquetFileStream->Close();
        if (!st2.ok()) {
          throw std::runtime_error(std::string("FinishOutputStream(): Error closing parquetFileStream: ") + st2.ToString());
        }
      }
      parquetFileStream.reset();
    }
    
    if (outFileStream) {
      if(!outFileStream->closed())
      {
        arrow::Status st2 = outFileStream->Close();
        if (!st2.ok()) {
          throw std::runtime_error(std::string("FinishOutputStream(): Error closing outFileStream: ") + st2.ToString());
        }
      }
      outFileStream.reset();
    }
    // writer_writing.store(false);
    return arrow::Status::OK();
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("FinishOutputStream(): Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("FinishOutputStream(): C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("FinishOutputStream(): C++ Runtime Error : ") + e.what());
  }
  catch(...){
    Rcpp::stop("FinishOutputStream(): Unknown Exception");
  }
  // std::thread fin_thread([this]()
  //                        {
  //                          // if (this->writer_threads.size() > 0)
  //                          // {
  //                          //   static_cast<void>(this->writer_threads[this->writer_threads.size() - 1].join());
  //                          //   //  for (auto &wrt_thread : this->writer_threads)
  //                          //   //  {
  //                          // 
  //                          //   //    static_cast<void>(wrt_thread.join());
  //                          //   //  }
  //                          // }
  //                          
  //                          for (auto &t : writer_threads) {
  //                            if (t.joinable())
  //                              t.join();
  //                          }
  //  
  //                          // static_cast<void>(this->rec_writer->Close());
  //                          this->writer_threads.clear();
  //                          this->rbv_batch->clear();
  //                          Rcpp::Rcout << "Done writing to file." << std::endl <<std::flush; //DEBUG
  //                          
  //                          Rcpp::Rcout << "here 1.7.2" << std::endl << std::flush; //DEBUG
  //                          // // Close compressed stream
  //                          // arrow::Status st1 = compressed_outstream->Close();
  //                          // Rcpp::Rcout << "here 1.7.3" << std::endl << std::flush; //DEBUG
  //                          // if (!st1.ok()) {
  //                          //   Rcpp::Rcerr << "Error closing compressed_outstream: " << st1.ToString() << std::endl;
  //                          // }
  //                          
  //                          Rcpp::Rcout << "here 1.7.4" << std::endl << std::flush; //DEBUG
  //                          // outFileStream->Flush();
  //                          // Rcpp::Rcout << "here 1.7.5" << std::endl << std::flush; //DEBUG
  //                          // // Optionally check outFileStream closed; but compressed close usually handled it
  //                          // if (outFileStream && !outFileStream->closed()) {
  //                          //   Rcpp::Rcout << "here 1.7.6" << std::endl << std::flush; //DEBUG
  //                          //   arrow::Status st2 = outFileStream->Close();
  //                          //   Rcpp::Rcout << "here 1.7.7" << std::endl << std::flush; //DEBUG
  //                          //   if (!st2.ok()) {
  //                          //     Rcpp::Rcerr << "Error closing outFileStream: " << st2.ToString() << std::endl;
  //                          //   }
  //                          // }
  //                          
  //                          // outFileStream.reset();
  //                        });
  // 
  // // fin_thread.detach();
  // if (fin_thread.joinable())
  //   fin_thread.join();
}

arrow::Status ArrowWrapper::FinishOutputStream()
{
  return pImpl->FinishOutputStream();
}

void ArrowWrapper::Impl::writer_loop() {
  try {
    
    // Create the Arrow RecordBatchWriter once
    auto writer_res = arrow::ipc::MakeFileWriter(outFileStream.get(), GetBLASTSchema(),
                                                 GetArrowIPCOptions(), GetBLASTMetadata());
    if (!writer_res.ok()) {
      {
        std::lock_guard<std::mutex> lk(writer_error_mtx);
        writer_error_msg = "MakeFileWriter failed: " + writer_res.status().ToString();
      }
      writer_failed.store(true);
      return;
    }
    rec_writer = writer_res.ValueOrDie();
    
    while (writer_running.load() || !rbv_batch->empty()) {
      // assert(!Progress::check_abort()); // R API calls inside threads will crash c++ runtime
      // Rcpp::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
      if (rbv_batch->size() > rb_batch_size.load() || rb_batch_size.load() == 0) {
        
        arrow::RecordBatchVector rbv_buffer;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_set_lock(&rbv_batchLock);
#endif
        rbv_buffer.swap(*rbv_batch);
        rbv_batch->clear();
        rbv_batch->shrink_to_fit();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&rbv_batchLock);
#endif
        
        switch(OutputFormat2Enum(output_format)){
        case ArrowWrapper::EOutputFormat::eIPC:
        case ArrowWrapper::EOutputFormat::eCSV: {
          for (std::shared_ptr<arrow::RecordBatch> &local_rb : rbv_buffer)
          {
            // assert(!Progress::check_abort()); // R API calls inside threads will crash c++ runtime
            // Rcpp::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
            if (local_rb->num_rows() <= 0){
              local_rb.reset();
              continue;
            }
            arrow::Status v = local_rb->ValidateFull();
            if (!v.ok()) {
              // log only — do not call Rcpp::stop from this thread
              std::lock_guard<std::mutex> lk(writer_error_mtx);
              writer_error_msg = "Invalid RecordBatch: " + v.ToString();
              writer_failed.store(true);
              break;
            }
            arrow::Status sts = rec_writer->WriteRecordBatch(*local_rb);
            if (!sts.ok()) {
              std::lock_guard<std::mutex> lk(writer_error_mtx);
              writer_error_msg = "WriteRecordBatch failed: " + sts.ToString();
              writer_failed.store(true);
              break;
            }
          }
          break;
        }
        case ArrowWrapper::EOutputFormat::eParquet :{
         auto table = arrow::Table::FromRecordBatches(GetBLASTSchema(), rbv_buffer).ValueOrDie();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
          omp_set_lock(&rec_writerLock);
#endif
          parquet::arrow::WriteTable(*table.get(),
                                     arrow::default_memory_pool(), parquetFileStream,
                                     /*chunk_size=*/ parquet_batch_size, parquet_props, parquet_arrow_props);
          parquetFileStream->Flush();
          arrow::Status st2 = parquetFileStream->Close();
          if (!st2.ok()) {
            std::lock_guard<std::mutex> lk(writer_error_mtx);
            writer_error_msg = "WriteRecordBatch failed: " + st2.ToString();
            writer_failed.store(true);
            break;
          }
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
          omp_unset_lock(&rec_writerLock);
#endif
          
          break;
        }
        case ArrowWrapper::EOutputFormat::unknown:
        default: {
          throw std::runtime_error(std::string("WriteBatch2File() - Unsupported output format. Supported values ipc/csv/parquet"));
        }
        }
        
        if(outFileStream)
          outFileStream->Flush();
        if(parquetFileStream)
          parquetFileStream->Flush();
        
        rbv_buffer.clear();
        rbv_buffer.shrink_to_fit();
        // rbv_not_full.notify_all();
        // rbv_not_full.notify_one();
      }
      else{
        continue;
      }

#ifdef ARROW_HAVE_MEMORY_POOL_RELEASE
      arrow::default_memory_pool()->ReleaseUnused();
#endif
      if (writer_failed.load()) break;
    } // while
    
    // Close writer gracefully
    if (rec_writer) {
      rec_writer->Close();
      rec_writer.reset();
    }
    if (outFileStream) {
      outFileStream->Flush();
      outFileStream->Close(); // if desired; check Arrow version semantics
      outFileStream.reset();
    }
    if(parquetFileStream)
      parquetFileStream->Flush();
  }
  catch (const std::exception &ex) {
    {
      std::lock_guard<std::mutex> lk(writer_error_mtx);
      writer_error_msg = std::string("writer_loop(): C++ exception: ") + ex.what();
    }
    writer_failed.store(true);
  }
  catch (...) {
    std::lock_guard<std::mutex> lk(writer_error_mtx);
    writer_error_msg = "writer_loop(): unknown exception";
    writer_failed.store(true);
  }
}


arrow::Status ArrowWrapper::Impl::WriteBatch2File()
{
  // DO NOT USE R API CALLS FROM WITHIN STD::THREAD AS IT LEADS TO STACK CORRUPTION
  try{
    if(writer_writing.load())
    {
      return arrow::Status::OK();
    }

    if(!save2file)
    {
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_set_lock(&rbv_batchLock);
  #endif
      rbv_batch->clear();
      rbv_batch->shrink_to_fit();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_unset_lock(&rbv_batchLock);
  #endif
      return arrow::Status::OK();
    }

    writer_writing.store(true);
    
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_set_lock(&rbv_batchLock);
  #endif
    auto is_rbv_empty = rbv_batch->empty();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_unset_lock(&rbv_batchLock);
  #endif
    if (is_rbv_empty) {
      // nothing to do
      return arrow::Status::OK();
    }

    arrow::RecordBatchVector rbv_buffer;
    
      // move all pending batches (or some bounded number) to local_buffer
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_set_lock(&rbv_batchLock);
  #endif
      rbv_buffer.reserve(rbv_batch->size());
      rbv_buffer.swap(*rbv_batch);
      rbv_batch->clear();
      rbv_batch->shrink_to_fit();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_unset_lock(&rbv_batchLock);
  #endif
    
    
    switch(OutputFormat2Enum(output_format)){
    case ArrowWrapper::EOutputFormat::eIPC:
    case ArrowWrapper::EOutputFormat::eCSV:
      {
        for (std::shared_ptr<arrow::RecordBatch> &rb : rbv_buffer) //const auto
        {
          // assert(!Progress::check_abort()); // R API calls inside threads will crash c++ runtime
          // Rcpp::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
          if (rb)
          {
            if (rb->num_rows() > 0)
            {
              arrow::Status rb_sts = rb->ValidateFull();
              if (rb_sts.ok())
              {
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                  omp_set_lock(&rec_writerLock);
    #endif
                  arrow::Status sts = rec_writer->WriteRecordBatch(*rb);
                  outFileStream->Flush();
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                  omp_unset_lock(&rec_writerLock);
    #endif
                  if (!sts.ok())
                  {
                    throw std::runtime_error(std::string("WriteBatch2File(): Error writing RB (CSV/IPC): ") + sts.message());
                  }  
                
              }
            }
            rb.reset();
          }
        }
        break;
      }
    case ArrowWrapper::EOutputFormat::eParquet :{
        auto table_ = arrow::Table::FromRecordBatches(GetBLASTSchema(), rbv_buffer);
        if(!table_.ok()){
          throw std::runtime_error(std::string("WriteBatch2File(): Error converting RBV to Table (Parquet): ") + table_.status().message() + std::string(" : ") + table_.status().detail()->ToString());
        }
        auto table = table_.ValueOrDie();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_set_lock(&rec_writerLock);
#endif
        parquet::arrow::WriteTable(*table.get(),
                                   arrow::default_memory_pool(), parquetFileStream,
                                   /*chunk_size=*/ parquet_batch_size, parquet_props, parquet_arrow_props);
        parquetFileStream->Flush();
        arrow::Status st2 = parquetFileStream->Close();
        if (!st2.ok()) {
          throw std::runtime_error(std::string("WriteBatch2File(): Error closing parquetFileStream: ") + st2.ToString());
        }
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&rec_writerLock);
#endif
        break;
      }
    case ArrowWrapper::EOutputFormat::unknown:
    default: {
      throw std::runtime_error(std::string("WriteBatch2File() - Unsupported output format. Supported values ipc/csv/parquet"));
    }
    }
    
    rbv_buffer.clear();
    rbv_buffer.shrink_to_fit();
    
    // release local big containers
    #ifdef ARROW_HAVE_MEMORY_POOL_RELEASE
      arrow::default_memory_pool()->ReleaseUnused();
    #endif

  return arrow::Status::OK();
  }
  catch(const Rcpp::exception &e){
    Rcpp::stop(std::string("WriteBatch2File() - Rcpp Exception : ") + e.what());
  }
  catch(const std::exception &e){
    Rcpp::stop(std::string("WriteBatch2File() - C++ Exception : ") + e.what());
  }
  catch(const std::runtime_error &e){
    Rcpp::stop(std::string("WriteBatch2File(): C++ Runtime Error : ") + e.what());
  }
  catch(...){
    Rcpp::stop( "WriteBatch2File() - Unknown Exception" );
  }
}

void CountCharacter_thread(const std::string &filename, char character, std::atomic<unsigned int> &count, size_t start, size_t end)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    // std::cerr << "Failed to open the file." << std::endl;
    REprintf("Failed to open the file.\n");
    return;
  }

  file.seekg(start, std::ios::beg);
  size_t chunkSize = end - start;
  std::vector<char> buffer(chunkSize);

  file.read(buffer.data(), chunkSize);

  for (size_t i = 0; i < chunkSize; ++i)
  {
    if (buffer[i] == character)
    {
      ++count;
    }
  }
}

int ArrowWrapper::Impl::CountCharacter(std::string filename, char character, int num_threads)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    // std::cerr << "Failed to open the file." << std::endl;
    REprintf("Failed to open the file.\n");
    return 1;
  }

  file.seekg(0, std::ios::end);
  size_t fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  std::atomic<unsigned int> totalCount(0);
  std::vector<std::thread> threads;

  size_t chunkSize = fileSize / num_threads;
  size_t start = 0;
  size_t end = chunkSize;

  for (int i = 0; i < num_threads - 1; ++i)
  {
    threads.emplace_back(CountCharacter_thread, filename, character, std::ref(totalCount), start, end);
    start = end;
    end += chunkSize;
  }

  // The last thread might have a slightly larger chunk
  threads.emplace_back(CountCharacter_thread, filename, character, std::ref(totalCount), start, fileSize);

  // Wait for all threads to finish
  for (auto &thread : threads)
  {
    thread.join();
  }
  
  // max_records = std::max(totalCount, max_records);
  unsigned int totalCount_ = totalCount;
  max_records = (max_records >= totalCount_) ? max_records : totalCount_;
  // unsigned int max_writer_threads_ = static_cast<int>(std::ceil((1 / (max_records / n_threads)) + 1));
  // max_writer_threads =  max_writer_threads_ < n_threads ? max_writer_threads_ : n_threads;
  
  return totalCount;
}
int ArrowWrapper::CountCharacter(std::string filename, char character, int num_threads)
{
  return pImpl->CountCharacter(filename, character, num_threads);
}

int ArrowWrapper::Impl::GetColumnCount(const std::string_view &filename, char delim)
{
  std::ifstream file(filename.data());
  if (!file.is_open())
  {
    // std::cerr << "Failed to open file: " << filename << std::endl;
    REprintf("Failed to open the file: %s \n", filename.data());
    return -1;
  }

  std::string line;
  if (std::getline(file, line))
  {
    std::stringstream ss(line);
    std::string column;
    int count = 0;
    while (std::getline(ss, column, delim))
    {
      count++;
    }
    return count;
  }
  else
  {
    // Rcpp::Rcerr << "File is empty: " << filename << std::endl;
    // std::cerr << "File is empty: " << filename << std::endl;
    REprintf("File is empty: %s \n", filename.data());
    return -1;
  }
}

// Quick guard: return true if entry appears to be a reasonable FASTA entry we can parse.
// - sv must include the leading '>' if it's a headered FASTA chunk
// - perform cheap checks only (no heavy conversions)
// static bool IsValidFastaEntry(std::string_view sv)
// {
//   // // Minimum size: header line '>' + at least 1 char and newline + some sequence
//   // if (sv.size() < 4) return false;
//   // 
//   // // Must start with '>'
//   // if (sv.front() != '>') return false;
//   // 
//   // // Find end of header line
//   size_t nlpos = sv.find_first_of("\r\n");
//   // if (nlpos == std::string_view::npos) return false;        // no newline -> bad
//   // if (nlpos <= 1) return false;                            // header too short (">" or ">\n")
//   
//   // Sequence region (after first EOL)
//   size_t seq_start = nlpos == std::string_view::npos ? 0 : nlpos + 1;
//   // skip additional \r\n sequences
//   while (seq_start < sv.size() && (sv[seq_start] == '\r' || sv[seq_start] == '\n')) ++seq_start;
//   // if (seq_start >= sv.size()) return false;                // header-only, no sequence
//   
//   std::string_view seq_sv = sv.substr(seq_start);
//   
//   // Cheap: detect embedded NULs or VERY high control character density
//   size_t n_null = 0;
//   size_t n_control = 0;
//   size_t n_letters = 0;
//   size_t n_bytes = seq_sv.size();
//   const unsigned char *bytes = reinterpret_cast<const unsigned char*>(seq_sv.data());
//   
//   for (size_t i = 0; i < n_bytes; ++i) {
//     unsigned char b = bytes[i];
//     if (b == 0) {
//       ++n_null;
//     } else if (b < 32 && b != '\t' && b != '\n' && b != '\r') {
//       ++n_control;
//     }
//     if (std::isalpha(b)) ++n_letters;
//   }
//   
//   // Reject if many NULs (likely binary or UTF-16) or many control characters
//   if (n_null > 0) return false;
//   if (n_control > n_bytes / 10) return false; // >10% weird controls -> skip
//   
//   // Reject if sequence has too few letters (not a biological sequence)
//   if (n_letters < 1) return false;
//   
//   // Heuristic UTF-16-like detection: many zero bytes in even or odd offsets
//   size_t zero_even = 0, zero_odd = 0;
//   size_t check_len = std::min<size_t>(n_bytes, 200);
//   for (size_t i = 0; i + 1 < check_len; ++i) {
//     if (bytes[i] == 0 && (i % 2 == 0)) ++zero_even;
//     if (bytes[i] == 0 && (i % 2 == 1)) ++zero_odd;
//   }
//   if ( (zero_even > (size_t) (check_len / 4)) || (zero_odd > (size_t) (check_len / 4)) ) {
//     // Looks like UTF-16/16-bit binary text -> skip (you could convert instead)
//     return false;
//   }
//   
//   // Optionally: max header length sanity check (avoid huge headers)
//   // if (nlpos > 1000) return false;
//   
//   // Looks OK
//   return true;
// }

// Function to trim leading whitespace
static void trim_left(std::string_view& sv) {
  size_t first_char = 0;
  while (first_char < sv.size() && std::isspace(static_cast<unsigned char>(sv[first_char]))) {
    first_char++;
  }
  sv.remove_prefix(first_char);
}

// Function to trim trailing whitespace
static void trim_right(std::string_view& sv) {
  size_t last_char = sv.size();
  while (last_char > 0 && std::isspace(static_cast<unsigned char>(sv[last_char - 1]))) {
    last_char--;
  }
  sv.remove_suffix(sv.size() - last_char);
}

// Function to trim both leading and trailing whitespace
static void trim(std::string_view& sv) {
  trim_left(sv);
  trim_right(sv);
}

std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::Impl::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{

  // Remove embedded NULs which can confuse C++ string handling / toolkit
  auto strip_nuls_and_controls = [](std::string &s) {
    s.erase(std::remove(s.begin(), s.end(), '\0'), s.end()); // remove NULs
    // optionally remove other control chars except newline if present in your flow
    s.erase(std::remove_if(s.begin(), s.end(),
                           [](unsigned char c){ return (c < 32 && c != '\n' && c != '\r' && c != '\t'); }), s.end());
  };
  
  std::string delim_str(delim);

  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> file_ptrs = MMapFile(filename, delim);

  char *start_of_file = std::get<1>(*file_ptrs).get();
  char *end_of_file = std::get<3>(*file_ptrs);

  auto is_trimmable_tail = [](unsigned char b) -> bool {
    // trim final: NUL, CR, LF, TAB. Also optionally trim other controls (be conservative).
    if (b == 0x00) return true;
    if (b == '\r' || b == '\n' || b == '\t') return true;
    // optionally: treat other low controls as trimmable but log them:
    // if (b < 32) return true;
    return false;
  };
  
  // Heuristic UTF-16 detection: many zero bytes on even or odd offsets near the beginning
  auto looks_like_utf16 = [&](const char* s, const char* e)->bool {
    const ptrdiff_t check_len = std::min<ptrdiff_t>(200, e - s);
    if (check_len < 8) return false;
    size_t zero_even = 0, zero_odd = 0;
    const unsigned char* u = reinterpret_cast<const unsigned char*>(s);
    for (ptrdiff_t i = 0; i + 1 < check_len; ++i) {
      if (u[i] == 0) { if ((i & 1) == 0) ++zero_even; else ++zero_odd; }
    }
    // If many zeros on one parity -> probably UTF-16LE/BE
    return (zero_even > (size_t)(check_len / 4)) || (zero_odd > (size_t)(check_len / 4));
  };
  
  // detect BOMs
  auto has_utf8_bom = [&](const char* s, const char* e)->bool {
    return (e - s) >= 3 && 
      (static_cast<unsigned char>(s[0]) == 0xEF &&
      static_cast<unsigned char>(s[1]) == 0xBB &&
      static_cast<unsigned char>(s[2]) == 0xBF);
  };
  auto has_utf16le_bom = [&](const char* s, const char* e)->bool {
    return (e - s) >= 2 &&
      static_cast<unsigned char>(s[0]) == 0xFF &&
      static_cast<unsigned char>(s[1]) == 0xFE;
  };
  auto has_utf16be_bom = [&](const char* s, const char* e)->bool {
    return (e - s) >= 2 &&
      static_cast<unsigned char>(s[0]) == 0xFE &&
      static_cast<unsigned char>(s[1]) == 0xFF;
  };
  
  // lambda: convert UTF-16 (LE or BE) in string_view -> UTF-8 std::string
  auto utf16_to_utf8_iconv = [](std::string_view sv, bool &converted) -> std::string {
    converted = false;
    if (sv.empty()) return std::string();
    
    // helpers to detect BOM / UTF-16 pattern
    auto has_bom_le = [&](const char* s, size_t n)->bool {
      return n >= 2 && (static_cast<unsigned char>(s[0]) == 0xFF)
      && (static_cast<unsigned char>(s[1]) == 0xFE);
    };
    auto has_bom_be = [&](const char* s, size_t n)->bool {
      return n >= 2 && (static_cast<unsigned char>(s[0]) == 0xFE)
      && (static_cast<unsigned char>(s[1]) == 0xFF);
    };
    auto looks_like_utf16 = [&](const char* s, size_t n)->int {
      // return 1 for LE, 2 for BE, 0 for not-utf16
      const size_t check_len = std::min<size_t>(n, 200);
      if (check_len < 8) return 0;
      size_t zero_even = 0, zero_odd = 0;
      const unsigned char* u = reinterpret_cast<const unsigned char*>(s);
      for (size_t i = 0; i + 1 < check_len; ++i) {
        if (u[i] == 0) {
          if ((i & 1) == 0) ++zero_even;
          else ++zero_odd;
        }
      }
      if (zero_odd > check_len/8) return 1; // many zeros on odd positions -> LE
      if (zero_even > check_len/8) return 2; // many zeros on even positions -> BE
      return 0;
    };
    
    const char* data = sv.data();
    size_t len = sv.size();
    size_t skip = 0;
    std::string src_encoding;
    
    if (has_bom_le(data, len))       { src_encoding = "UTF-16LE"; skip = 2; }
    else if (has_bom_be(data, len))  { src_encoding = "UTF-16BE"; skip = 2; }
    else {
      int heuristic = looks_like_utf16(data, len);
      if (heuristic == 1) { src_encoding = "UTF-16LE"; skip = 0; }
      else if (heuristic == 2) { src_encoding = "UTF-16BE"; skip = 0; }
      else {
        // Not UTF-16 according to our heuristics -> return original bytes unchanged
        converted = false;
        return std::string(sv); // copy original data
      }
    }
    
    // copy input into a mutable std::string (iconv wants mutable pointers)
    std::string inbuf;
    if (skip > 0 && len > skip) {
      inbuf.assign(data + skip, len - skip);
    } else if (skip > 0) {
      // only BOM present and no data -> nothing to convert
      converted = true;
      return std::string();
    } else {
      inbuf.assign(data, len);
    }
    
    // Prepare iconv
    iconv_t cd = iconv_open("UTF-8", src_encoding.c_str());
    if (cd == (iconv_t)-1) {
      // iconv not available for that encoding or other error.
      converted = false;
      return std::string(sv);
    }
    
    // Estimate output size: roughly 4 bytes per UTF-16 code unit (over-estimate safe)
    size_t inbytesleft = inbuf.size();
    size_t outbuf_capacity = (inbytesleft + 1) * 3 + 32;
    std::string out;
    out.resize(outbuf_capacity);
    char* inptr = inbuf.empty() ? nullptr : &inbuf[0];
    char* outptr = out.empty() ? nullptr : &out[0];
    size_t outbytesleft = outbuf_capacity;
    
    // iconv signature uses char** on many platforms
    size_t res = iconv(cd, &inptr, &inbytesleft, &outptr, &outbytesleft);
    if (res == (size_t)-1) {
      // conversion error: we could try incremental loop, but for simplicity return original
      // Optionally you can inspect errno (EILSEQ, EINVAL, E2BIG)
      iconv_close(cd);
      converted = false;
      return std::string(sv);
    }
    
    // construct final string from used bytes
    size_t out_used = outbuf_capacity - outbytesleft;
    out.resize(out_used);
    iconv_close(cd);
    
    converted = true;
    return out;
  };
  
  char* file_start = start_of_file;
  char* file_end = end_of_file;   // original one-past-last or your current end pointer
  
  // Quick encoding checks:
  if (has_utf8_bom(file_start, file_end)) {
    Rcpp::Rcout << "Detected UTF-8 BOM; skipping BOM bytes." << std::endl << std::flush;
    file_start += 3; // skip BOM for parsing
  } else if (has_utf16le_bom(file_start, file_end) || has_utf16be_bom(file_start, file_end) ||
    looks_like_utf16(file_start, file_end)) {
    // UTF-16: do NOT try to parse as ASCII. Convert the file to UTF-8 (iconv/boost/ICU) or fail.
    Rcpp::Rcerr << "Detected UTF-16 (or many NULs) - convert file to UTF-8 before processing." << std::endl << std::flush;
    // You can either call a converter here (iconv) or bail out.
    // return error / throw / log for the caller.
  }
  
  // Trim trailing junk from the end of file (do not write to mmap):
  char* adjusted_end = file_end; //const
  while (adjusted_end > file_start) {
    unsigned char last = static_cast<unsigned char>(*(adjusted_end - 1));
    if (is_trimmable_tail(last)) {
      --adjusted_end;
    } else {
      break;
    }
  }
 
  // if (adjusted_end != file_end) { //DEBUG
  //   Rcpp::Rcout << "Trimmed " << (file_end - adjusted_end) << " trailing bytes from file end." << std::endl << std::flush; //DEBUG
  // } //DEBUG
  
  char *p = start_of_file;

  // unsigned char b1 = static_cast<unsigned char>(*(adjusted_end - 3));
  // unsigned char b2 = static_cast<unsigned char>(*(adjusted_end - 2));
  // unsigned char b3 = static_cast<unsigned char>(*(adjusted_end - 1));
  // unsigned char b4 = static_cast<unsigned char>(*(adjusted_end));
  // 
  // Rcpp::Rcout << "HEX tail: "
  //             << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b1) << " "
  //             << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b2) << " "
  //             << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b3)
  //             << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b4)
  //             << std::dec << std::endl << std::flush; // restore decimal
  // Rcpp::Rcout << adjusted_end - 3 << adjusted_end - 2 << adjusted_end - 1 << adjusted_end[0] << std::endl << std::flush; //DEBUG
  
  arrow::RecordBatchVector ret_results;

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_lock_t pLock;
  omp_lock_t ret_resultsLock;
  omp_init_lock(&pLock);
  omp_init_lock(&ret_resultsLock);
#endif

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp parallel num_threads(num_threads) shared(adjusted_end, start_of_file, delim) // entry_ptr_vec
#endif
  {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp for schedule(dynamic) nowait // schedule(dynamic)
#endif
    for (int i = 0; i < num_threads; ++i)
    {

      // assert(!Progress::check_abort());
      Rcpp::checkUserInterrupt();
      //  Get the thread-specific range to process
      size_t chunk_size = (adjusted_end - start_of_file) / num_threads;
      char *thread_start = start_of_file + i * chunk_size;
      char *thread_end = (i == num_threads - 1) ? adjusted_end : (thread_start + chunk_size);

      if (thread_start != start_of_file)
      {
        while (strncmp(thread_start, delim, strlen(delim)) != 0)
        {
          --thread_start;
        }
      }
      if (thread_end != adjusted_end)
      {
        while (strncmp(thread_end, delim, strlen(delim)) != 0)
        {
          ++thread_end;
        }
        thread_end = thread_end - 1;
      }
      // Process the entries within the thread's range
      char *entryStart = nullptr;
      char *entryEnd = nullptr;

      const size_t delim_len = std::strlen(delim);
      
      // p is not used now; we create piter
      const char* piter = thread_start;
      // bool endLoop = false;
      while (piter < thread_end && piter < adjusted_end)
        {
        // find start of next delimiter inside [piter, thread_end)
        const char* entryStart = std::search(piter, static_cast<const char*>(thread_end),
                                             delim, delim + delim_len);
        if (entryStart == thread_end) {
          // no more delimiters in this chunk
          break;
        }
        
        // find next delimiter after this one to mark end-of-entry
        const char* nextDelim = std::search(entryStart + delim_len, static_cast<const char*>(thread_end),
                                            delim, delim + delim_len);
        // entryEnd is one-past-last-byte of the entry region; may be thread_end (>= adjusted_end)
        const char* entryEnd = (nextDelim == thread_end) ? thread_end : nextDelim;
        
        // sanity: entryEnd must be after entryStart
        if (entryEnd <= entryStart) {
          // malformed or empty; advance to avoid infinite loop
          piter = (entryEnd < thread_end) ? entryEnd + 1 : thread_end;
          continue;
        }
        
        // compute length; safe because entryEnd >= entryStart
        size_t entry_len = static_cast<size_t>(entryEnd - entryStart);
        if (entry_len == 0) {
          piter = entryEnd;
          if (entryEnd >= adjusted_end) break;
          continue;
        }
        
        // build string_view and trim trailing CR/LF without ever dereferencing entryEnd
        std::string_view sv_entry(entryStart, entry_len);
        trim(sv_entry);
        while (!sv_entry.empty() && (sv_entry.back() == '\r' || sv_entry.back() == '\n')) {
          sv_entry.remove_suffix(1);
        }
        
        // // Debug printing (do NOT dereference entryEnd unless entryEnd < thread_end)
        // Rcpp::Rcout << "entryStart: " << static_cast<const void*>(entryStart)
        //             << " entryEnd: "  << static_cast<const void*>(entryEnd)
        //             << " adjusted_end: " << static_cast<const void*>(adjusted_end)
        //             << " entryEnd >= adjusted_end? " << (entryEnd >= adjusted_end ? "yes" : "no")
        //             << std::endl << std::flush; // DEBUG
        // Rcpp::Rcout << "HEX: "
        //             << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(static_cast<unsigned char>(*(entryEnd)))
        //             << std::dec << std::endl << std::flush; // restore decimal
        // Rcpp::Rcout << "DEC: " << static_cast<int>(static_cast<unsigned char>(*(entryEnd))) << std::endl << std::flush; // restore decimal
        // 
        // if (!sv_entry.empty()) {
        //   Rcpp::Rcout << "first: " << sv_entry.front()
        //               << " last: "  << sv_entry.back()
        //               << " len: " << sv_entry.size() << std::endl << std::flush;
        // }
        
        bool converted = false;
        std::string utf8_str = utf16_to_utf8_iconv(sv_entry, converted);
        std::string_view to_parse = converted ? std::string_view(utf8_str) : sv_entry;
        
        // Parse and skip empty/malformed sequences
        FastaSequenceData conv_entry = CastToType(to_parse); //sv_entry
        
        // strip_nuls_and_controls(conv_entry.seq);
        
        conv_entry.seq.erase(std::remove_if(conv_entry.seq.begin(), conv_entry.seq.end(),
                                      [](char c){ 
                                        // allow A-Z, a-z, '*' (you may customize for DNA vs protein)
                                        return !((std::isalpha((unsigned char)c) || c=='*') && static_cast<int>(static_cast<unsigned char>(c)) > 10); 
                                      }), conv_entry.seq.end());
        
        if (conv_entry.seq.empty()) {
          piter = entryEnd;
          if (entryEnd >= adjusted_end) break;
          continue;
        }
        
        AddRecordCount();
        if (Entry_callback) {
          std::shared_ptr<arrow::RecordBatchVector> tmp_result =
            Entry_callback(std::make_shared<FastaSequenceData>(conv_entry));
          if (return_values) {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
            omp_set_lock(&ret_resultsLock);
#endif
            ret_results.insert(ret_results.end(), tmp_result->begin(), tmp_result->end());
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
            omp_unset_lock(&ret_resultsLock);
#endif
          }
          // else {
          tmp_result->clear();
          tmp_result->shrink_to_fit();
          // }
        }
        
        // advance to next search position
        piter = entryEnd;
        // if we hit the real file-end, we're done
        if (entryEnd >= adjusted_end) {
          break;
        }
      }
      
//       for (char *p = thread_start; p < thread_end; ++p)
//       {
//         assert(!Progress::check_abort());
// 
//         if (strncmp(p, delim, strlen(delim)) == 0)
//         {
//           entryStart = strstr(p, delim);
//           entryEnd = strstr(p + 1, delim);
// 
//           if (entryStart == nullptr) {
//             // nothing to do; defensive: avoid using null pointer
//             break;
//           }
// 
//           // If no next delimiter was found, treat end_of_file as the end (exclusive)
//           if (entryEnd == nullptr) {
//             entryEnd = end_of_file - 1;
//           }
// 
//           if (entryEnd == "\n" || entryEnd == "\r") {
//             entryEnd--;
//           }
// 
//           if (entryStart == end_of_file - 1 && entryEnd == end_of_file - 1) {
//             Rcpp::Rcout << "EOF" << std::endl << std::flush; //DEBUG
//             break;
//           }
//           
//           // Sanity checks: entryEnd must be strictly after entryStart
//           if (entryEnd <= entryStart) {
//             // malformed or empty entry, skip / stop to avoid infinite loop
//             break;
//           }
// 
//           // if entryStart is at or beyond EOF, stop
//           if (entryStart >= end_of_file && entryEnd >= end_of_file) {
//             break;
//           }
// 
//           // if (entryEnd == nullptr)
//           // {
//           //   // Handle the case where the delimiter is not found within the (final) thread's range
//           //   entryEnd = end_of_file - 1;
//           // }
// 
// 
//           // Process the entry from entryStart to entryEnd
//           // // std::string full_entry(entryStart, entryEnd - entryStart - 1);
//           // std::string_view sv_entry(entryStart, (entryEnd - entryStart - 1));
//           // Now safe to form the string_view: length is entryEnd - entryStart
//           size_t entry_len = static_cast<size_t>(entryEnd - entryStart);
//           std::string_view sv_entry(entryStart, entry_len);
//            trim(sv_entry);
//           Rcpp::Rcout << "+++++++++++++" << std::endl << std::flush; //DEBUG
//           Rcpp::Rcout << entryStart[0] << std::endl << std::flush; //DEBUG
//           Rcpp::Rcout << entryEnd[0] << std::endl << std::flush; //DEBUG
//           Rcpp::Rcout << "HEX: "
//                       << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(static_cast<unsigned char>(*(entryEnd)))
//                       << std::dec << std::endl << std::flush; // restore decimal
//           Rcpp::Rcout << "DEC: " << static_cast<int>(static_cast<unsigned char>(*(entryEnd))) << std::endl << std::flush; // restore decimal
//           Rcpp::Rcout << entry_len << std::endl << std::flush; //DEBUG
//           Rcpp::Rcout << "=============" << std::endl << std::flush; //DEBUG
//           // if (!IsValidFastaEntry(sv_entry)) {
//           //   // Rcpp::Rcout << "Skipping malformed entry (pos " << (entryStart - start_of_file) << ")\n";
//           //   continue;
//           // }
// 
//           if (!sv_entry.empty()) //full_entry
//           {
//             AddRecordCount();
//             const FastaSequenceData conv_entry = CastToType(sv_entry); //full_entry
//             if (Entry_callback != nullptr && !conv_entry.seq.empty())
//             {
//               std::shared_ptr<arrow::RecordBatchVector> tmp_result = Entry_callback(std::make_shared<FastaSequenceData>(conv_entry));
// 
//               if (return_values)
//               {
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                 omp_set_lock(&ret_resultsLock);
// #endif
// 
//                 ret_results.insert(ret_results.end(), tmp_result->begin(), tmp_result->end());
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
//                 omp_unset_lock(&ret_resultsLock);
// #endif
//               }
//               // else
//               // {
//                 tmp_result->clear();
//               // }
//             }
//           }
//           p = entryEnd - 1; // Move to the next position after the delimiter
//         }
// 
//         // // move the outer loop pointer to just before the end delimiter so the for-loop
//         // // increment will continue scanning past this entry
//         // if (entryEnd < end_of_file)
//         //   p = entryEnd - 1;
//         // else
//         //   p = end_of_file - 1;
// 
//       }
    }
  }
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp barrier
#endif

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_destroy_lock(&pLock);
  omp_destroy_lock(&ret_resultsLock);
#endif

  
  return std::make_shared<arrow::RecordBatchVector>(ret_results);
}
std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return pImpl->SplitFilesIntoEntries(filename, delim, num_threads, Entry_callback, return_values);
}