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
    arrow_LFS = arrow::fs::LocalFileSystem();
    std::string username = "";
  #if defined(linux) //|| defined(MINGW32)
    username = get_username_safe(); //getlogin();
  #endif
  
    if (username.empty())
    {
      username = "Unknown";
    }
  
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
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    ipc_options.use_threads = true;
#else
    ipc_options.use_threads = false;
#endif
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
  
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    parquet_arrow_props = ArrowWriterProperties::Builder().store_schema()
                          ->set_use_threads(true)
                          ->set_engine_version(ArrowWriterProperties::V2)
                          ->build();
#else
    parquet_arrow_props = ArrowWriterProperties::Builder().store_schema()
      ->set_use_threads(false)
      ->set_engine_version(ArrowWriterProperties::V2)
      ->build();
#endif
  
    rbv_batch = std::make_shared<arrow::RecordBatchVector>();
    blast_metadata = std::make_shared<arrow::KeyValueMetadata>();
    AddFASTAMetadata("format", "Arrow IPC/Parquet");
    AddFASTAMetadata("Created By", username);
    AddFASTAMetadata("R package", "QuickBLAST");
    fasta_schema = arrow::schema({arrow::field("index", arrow::int64()), arrow::field("header", arrow::utf8()), arrow::field("sequence", arrow::utf8())});
  
    seq_info_type = arrow::struct_({
        arrow::field("num_alignments", arrow::int64()),
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

  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_init_lock(&rec_countLock);
    omp_init_lock(&proc_rec_countLock);
    omp_init_lock(&writer_threadsLock);
    omp_init_lock(&rbv_batchLock);
    omp_init_lock(&rec_writerLock);
  #endif

    Rcpp::Rcout << std::flush;
  }
  catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[ArrowWrapper::Impl()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush;
  }
  catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[ArrowWrapper::Impl()]: Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }
  catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[ArrowWrapper::Impl()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }
  catch(...){
    Rcpp::Rcerr << "[ArrowWrapper::Impl()]: Unknown Exception"  << std::endl << std::flush;
  }
}

ArrowWrapper::Impl::~Impl()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_destroy_lock(&rec_countLock);
  omp_destroy_lock(&proc_rec_countLock);
  omp_destroy_lock(&writer_threadsLock);
  omp_destroy_lock(&rbv_batchLock);
  omp_destroy_lock(&rec_writerLock);
#endif
}

ArrowWrapper::ArrowWrapper()
    : pImpl(std::make_unique<Impl>()) {}

// Destructor: default implementation (pImpl will automatically clean up)
ArrowWrapper::~ArrowWrapper() = default;

void ArrowWrapper::Impl::SetBatchSize(unsigned int batch_size)
{
  this->rb_batch_size.store(batch_size);
  if(batch_size > 0)
  {
    csv_options.batch_size = batch_size;
    parquet_batch_size = batch_size;
  }
  if(verbose)
    Rcpp::Rcout << "Batch Size: " << batch_size << std::endl << std::flush; //DEBUG
}
void ArrowWrapper::SetBatchSize(unsigned int batch_size)
{
  pImpl->SetBatchSize(batch_size);
}

void ArrowWrapper::Impl::SetVerbosity(bool verbose)
{
  this->verbose = verbose;
}
void ArrowWrapper::SetVerbosity(bool verbose)
{
  pImpl->SetVerbosity(verbose);
}

unsigned int ArrowWrapper::Impl::GetBatchSize()
{
  return this->rb_batch_size.load();
}
unsigned int ArrowWrapper::GetBatchSize()
{
  return pImpl->GetBatchSize();
}

unsigned int ArrowWrapper::Impl::GetRecordCount() { return rec_count; }
unsigned int ArrowWrapper::GetRecordCount() { return pImpl->GetRecordCount(); }

unsigned int ArrowWrapper::Impl::GetPendingRecordCount(){
  return max_records - proc_rec_count;
}
unsigned int ArrowWrapper::GetPendingRecordCount(){
  return pImpl->GetPendingRecordCount();
}

void ArrowWrapper::Impl::ResetProcRecordCount()
{
  proc_rec_count = 1;
}
void ArrowWrapper::ResetProcRecordCount()
{
  pImpl->ResetProcRecordCount();
}

void ArrowWrapper::Impl::AddProcRecordCount()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp atomic update
  proc_rec_count++;
#else
  proc_rec_count++;
#endif
}
void ArrowWrapper::AddProcRecordCount()
{
  pImpl->AddProcRecordCount();
}

unsigned int ArrowWrapper::Impl::GetProcRecordCount() { return proc_rec_count; }
unsigned int ArrowWrapper::GetProcRecordCount() { return pImpl->GetProcRecordCount(); }

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
#pragma omp atomic update
  rec_count++;
#else
  rec_count++;
#endif
}
void ArrowWrapper::AddRecordCount()
{
  pImpl->AddRecordCount();
}

void ArrowWrapper::SetBLASTSeqLimit(unsigned int seq_limit){
  pImpl->SetBLASTSeqLimit(seq_limit);
}
void ArrowWrapper::Impl::SetBLASTSeqLimit(unsigned int seq_limit){
  blast_sequence_limit = seq_limit;
}
unsigned int ArrowWrapper::GetBLASTSeqLimit(){
  return pImpl->GetBLASTSeqLimit();
}
unsigned int ArrowWrapper::Impl::GetBLASTSeqLimit(){
  return blast_sequence_limit;
}

void ArrowWrapper::Impl::SetThreadCount(int num_threads)
{
  if (num_threads > 1)
  {
    ipc_options.use_threads = true;
  }
  else
  {
    ipc_options.use_threads = false;
  }

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp atomic read
  n_threads = num_threads;
#else
  n_threads = num_threads;
#endif
}
void ArrowWrapper::SetThreadCount(int num_threads)
{
  pImpl->SetThreadCount(num_threads);
}

FastaSequenceData ArrowWrapper::Impl::CastToType(const std::string_view &full_entry_sv)
{
  // full_entry_sv points inside the mmap buffer; DO NOT store
  // string_view beyond the lifetime of that buffer.
  FastaSequenceData fasta_data;
  fasta_data.rec_no = GetRecordCount(); 
  // Expect format:
  // >header-line\r?\n
  // sequence lines (maybe many, with whitespace)
  // We avoid regex; do simple manual parsing.

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
  seq.shrink_to_fit();
  
  fasta_data.seq.swap(seq); // move efficiently
  
  // std::cout << fasta_data.rec_no << ">" << fasta_data.header << std::endl << std::flush; //DEBUG
  // std::cout << fasta_data.seq << std::endl << std::flush; //DEBUG
  
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

arrow::Status ArrowWrapper::Impl::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
{
  try{
    RcppThread::checkUserInterrupt();
    if(!save2file)
    {
      return arrow::Status::OK();
    }
    
    arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
    if (rb_)
    {
 
 //      std::cout << "\rPRE: AddRB2Batch(): Writer threads : " << writer_threads.size() << " : " << rbv_batch->size() << " : " << rb_batch_size.load(std::memory_order_acquire) << " : " << proc_rec_count << " : " << max_records << " : " << is_writing_pre << "...." << std::flush; //DEBUG //target //std::endl
 
     if(rb_->num_rows() == 0)
       return arrow::Status::OK();

  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_set_lock(&rbv_batchLock);
  #endif
      rbv_batch->emplace_back(std::move(rb_));
      unsigned int ret_size = rbv_batch->size();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
      omp_unset_lock(&rbv_batchLock);
  #endif
     
     auto tmp_batch_size = rb_batch_size.load(std::memory_order_acquire);
      if((ret_size < tmp_batch_size && tmp_batch_size != 0) && blast_sequence_limit != 0){
        return arrow::Status::OK(); 
      }
      
     if(writer_writing.load(std::memory_order_acquire)){
       if(tmp_batch_size > 1)
         tmp_batch_size--;
       rb_batch_size.store(tmp_batch_size, std::memory_order_release);
       return arrow::Status::OK();
     }else{
       unsigned int max_proc_diff = max_records - proc_rec_count;
       max_proc_diff = max_proc_diff == 0 ? 1 : max_proc_diff; //Checking underflow;
       if(tmp_batch_size > max_proc_diff)
         tmp_batch_size = max_proc_diff;
       else
         tmp_batch_size++;
     }
     rb_batch_size.store(tmp_batch_size, std::memory_order_release);
     
     {
        safe_jthread write_thread(std::thread([this]()
          {
              try{
                static_cast<void>(this->WriteBatch2File());
                writer_writing.store(false, std::memory_order_release);
              }catch(const std::runtime_error &e){
                {
                  std::lock_guard<std::mutex> lk(writer_error_mtx);
                  writer_error_msg = std::string("[thread::WriteBatch2File()]: C++ Runtime Exception : ") + e.what();
                }
                writer_failed.store(true, std::memory_order_release);
              }catch(const Rcpp::exception &e){
                {
                  std::lock_guard<std::mutex> lk(writer_error_mtx);
                  writer_error_msg = std::string("[thread::WriteBatch2File()]: Rcpp Exception : ") + e.what();
                }
                writer_failed.store(true, std::memory_order_release);
              }catch(const std::exception &e){
                {
                  std::lock_guard<std::mutex> lk(writer_error_mtx);
                  writer_error_msg = std::string("[thread::WriteBatch2File()]: C++ Exception : ") + e.what();
                }
                writer_failed.store(true, std::memory_order_release);
              }catch(...){
                {
                  std::lock_guard<std::mutex> lk(writer_error_mtx);
                  writer_error_msg = "[thread::WriteBatch2File()]: Unknown Exception";
                }
                writer_failed.store(true, std::memory_order_release);
              }
              writer_writing.store(false, std::memory_order_release);
          }));
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
          omp_set_lock(&writer_threadsLock);
  #endif
          static_cast<void>(writer_threads.emplace_back(std::move(write_thread)));
          // static_cast<void>(writer_threads.emplace_back(write_thread));
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
          omp_unset_lock(&writer_threadsLock);
  #endif
          waiting4writer_cond.notify_all();
      }
      
    }
    
    return arrow::Status::OK();
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[AddRB2Batch()]: C++ Runtime Exception : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[AddRB2Batch()]: Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[AddRB2Batch()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }catch(...){
    Rcpp::Rcerr << "[AddRB2Batch()]: Unknown Exception" << std::endl << std::flush;
  }
  return arrow::Status::Invalid("[AddRB2Batch()]: Caught an Exception.");
}

arrow::Status ArrowWrapper::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
{
  return pImpl->AddRB2Batch(rb_);
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
      save2file = false;
    }else{
      save2file = true;
      output_filename = outFile;
      if(verbose){
        Rcpp::Rcout << "Writing to : " << output_filename << std::endl << std::flush; //DEBUG
        Rcpp::Rcout << "Output Format : " << output_format << std::endl << std::flush; //DEBUG
      }
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
            Rcpp::Rcerr << std::string("[CreateOutputStream()] - Could not open append stream to ") + output_filename<< std::endl << std::flush;
            return outFileStream_res.status();
          }
          outFileStream = outFileStream_res.ValueOrDie();
          auto writer_ = arrow::ipc::MakeFileWriter(outFileStream.get(), GetBLASTSchema(), GetArrowIPCOptions(), GetBLASTMetadata());
          if (!writer_.ok())
          {
            Rcpp::Rcerr << std::string("WriteBatch2File() - Error initiating IPC file writer: ") + writer_.status().message() << std::endl << std::flush;
            return writer_.status();
          }
          rec_writer = writer_.ValueOrDie();
          break;
        }
      case ArrowWrapper::EOutputFormat::unknown:
      default: {
        return arrow::Status::Invalid((std::string("CreateOutputStream() - Unsupported output format. Supported values ipc/csv/parquet")));
      }
      case ArrowWrapper::EOutputFormat::eCSV: {
          auto outFileStream_res = arrow_LFS.OpenAppendStream(output_filename, blast_metadata);
          if(!outFileStream_res.ok()){
            Rcpp::Rcerr << std::string("CreateOutputStream() - Could not open append stream to ") + output_filename << std::endl << std::flush;
            return outFileStream_res.status();
          }
          outFileStream = outFileStream_res.ValueOrDie();
          auto writer_ = arrow::csv::MakeCSVWriter(outFileStream.get(), GetBLASTSchema(), GetArrowCSVOptions());
          if (!writer_.ok())
          {
            Rcpp::Rcerr << std::string("CreateOutputStream() - Error initiating CSV file writer: ") + writer_.status().message() << std::endl << std::flush;
            return writer_.status();
          }
          rec_writer = writer_.ValueOrDie();
          break;
        }
      case ArrowWrapper::EOutputFormat::eParquet: {
          ARROW_ASSIGN_OR_RAISE(parquetFileStream,  arrow::io::FileOutputStream::Open(output_filename));
          parquet_writer = std::move(parquet::arrow::FileWriter::Open(*GetBLASTSchema(),
                                                        arrow::default_memory_pool(), parquetFileStream,
                                                        parquet_props, parquet_arrow_props).ValueOrDie()); 
          break;
        }
      }

      finisher_thread = std::thread([this]() {
        try {
          while ( (writer_running.load() || !writer_threads.empty()) && !writer_failed.load(std::memory_order_acquire) && !writer_finishing.load(std::memory_order_acquire)) {
            // assert(!Progress::check_abort()); // R API calls inside threads will crash c++ runtime
            RcppThread::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
            // 1. Declare an empty, non-joinable thread object (NOT a reference)
            std::thread thread_to_join; 
            
            {
              // scoped lock for safe access to writer_threads
              std::unique_lock<std::mutex> lk(writer_threads_mutex); 
              
              if (writer_threads.empty()) {
                // wait until there is a thread to join, or the writer stops/fails
                waiting4writer_cond.wait(lk, [this]() {
                  return !writer_running.load(std::memory_order_acquire) || 
                    !writer_threads.empty() || 
                    writer_failed.load(std::memory_order_acquire);
                });
                
                // The unique_lock will automatically unlock when 'continue' jumps to the next loop iteration
                continue; 
              }
              
              // 2. Move the thread out of the vector and into our local variable
              // (Assuming writer_threads holds your safe_jthread, so .get() accesses the std::thread)
              thread_to_join = std::move(writer_threads.front().get());
              // thread_to_join = std::move(writer_threads.front());
              
              // 3. Erase the first element. Because we moved the data out, this just safely deletes an empty shell.
              writer_threads.erase(writer_threads.begin());
              
            } // 4. unique_lock goes out of scope and automatically unlocks the mutex here
            
            // 5. Join outside the lock (so we don't block other threads while waiting)
            if (thread_to_join.joinable()) {
              thread_to_join.join();
            }
          }
          writer_finishing.store(true, std::memory_order_release);
          finishing_cond.notify_all();
        }catch(const std::runtime_error &e){
          {
            std::lock_guard<std::mutex> lk(writer_error_mtx);
            writer_error_msg = std::string("finisher_thread(): C++ Runtime Error : ") + e.what();
          }
          writer_failed.store(true, std::memory_order_release);
        }catch(const Rcpp::exception &e){
          {
            std::lock_guard<std::mutex> lk(writer_error_mtx);
            writer_error_msg = std::string("finisher_thread(): Rcpp Exception : ") + e.what();
          }
          writer_failed.store(true, std::memory_order_release);
        }catch(const std::exception &e){
          {
            std::lock_guard<std::mutex> lk(writer_error_mtx);
            writer_error_msg = std::string("finisher_thread(): C++ Exception : ") + e.what();
          }
          writer_failed.store(true, std::memory_order_release);
        }catch(...){
          {
            std::lock_guard<std::mutex> lk(writer_error_mtx);
            writer_error_msg = "finisher_thread(): Unknown Exception";
          }
          writer_failed.store(true, std::memory_order_release);
        }
        writer_finishing.store(true, std::memory_order_release);
        finishing_cond.notify_all();
      }); //std::move()
      finisher_thread.detach();   
    }
    return arrow::Status::OK();
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[CreateOutputStream()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[CreateOutputStream()]: Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[CreateOutputStream()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }catch(...){
    Rcpp::Rcerr << "[CreateOutputStream()]: Unknown Exception" << std::endl << std::flush;
  }
  return arrow::Status::Invalid("[AddRB2Batch()]: Caught an Exception.");
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
    Rcpp::Rcerr << "[MMapFile()] Error: Failed to open file:" << filename.data() << std::endl << std::flush;
    return nullptr;
  }

  // Get the file size
  fileSize = GetFileSize(file_ptr);
#if defined(linux) || defined(MINGW32)
  char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file_ptr), 0)); // MAP_SHARED
#else
  char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, _fileno(file_ptr), 0)); // MAP_SHARED
#endif // linux

  if (fileData_ptr == MAP_FAILED)
  {
    REprintf("Error: Failed to map file : %s \n", filename.data());
    fclose(file_ptr);
    return nullptr;
  }

  end_of_file_ptr = fileData_ptr + fileSize;
  
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
  
  char* file_start = fileData_ptr; //file_ptr;
  char* file_end = end_of_file_ptr;   // original one-past-last or your current end pointer
  
  // Quick encoding checks:
  if (has_utf8_bom(file_start, file_end)) {
    if(verbose){
      Rcpp::Rcout << "Detected UTF-8 BOM; skipping BOM bytes." << std::endl << std::flush;
    }
    file_start += 3; // skip BOM for parsing
  } else if (has_utf16le_bom(file_start, file_end) || has_utf16be_bom(file_start, file_end) ||
    looks_like_utf16(file_start, file_end)) {
    // UTF-16: do NOT try to parse as ASCII. Convert the file to UTF-8 (iconv/boost/ICU) or fail.
    if(verbose){
      Rcpp::Rcerr << "Detected UTF-16 (or many NULs) - convert file to UTF-8 before processing." << std::endl << std::flush;
    }
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
  
  return std::make_shared<std::tuple<FILE *, std::shared_ptr<char>, long, char *>>(std::make_tuple(file_ptr, std::shared_ptr<char>(fileData_ptr, [fileSize, file_ptr](char *ptr)
                                                                                                                                   {
    munmap(ptr, fileSize);
    fclose(file_ptr); }),
    fileSize, 
    // end_of_file_ptr));
    adjusted_end));
}

std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> ArrowWrapper::MMapFile(const std::string_view &filename, const char *delim){
  return pImpl->MMapFile(filename, delim);
}

long ArrowWrapper::GetFileSize(FILE *file_ptr){
  return pImpl->GetFileSize(file_ptr);
}

long ArrowWrapper::Impl::GetFileSize(FILE *file_ptr)
{
  fseek(file_ptr, 0, SEEK_END);
  long fileSize = ftell(file_ptr);
  rewind(file_ptr);
  return fileSize;
}

std::shared_ptr<std::list<FastaSequenceData>> ArrowWrapper::FetchRecordByBatch(const std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> &file_ptr, unsigned int batch_size, unsigned int from_rec, const char *delim){
  return pImpl->FetchRecordByBatch(file_ptr, batch_size, from_rec, delim);
}

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

std::shared_ptr<std::list<FastaSequenceData>> ArrowWrapper::Impl::FetchRecordByBatch(const std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> &file_ptr, unsigned int batch_size, unsigned int from_rec, const char *delim){
  
  std::string delim_str(delim);
  const size_t delim_len = std::strlen(delim);
  using BM = std::boyer_moore_horspool_searcher<const char*>; //std::boyer_moore_searcher<const char*>;
  auto bm = BM(delim, delim + delim_len);
  
  const char *start_fptr = std::get<1>(*file_ptr).get();
  const char *end_fptr = std::get<3>(*file_ptr);
  char *entryStart = nullptr;
  char *entryEnd = nullptr;

  if(batch_size == 0) batch_size++;
  const unsigned int to_rec = from_rec + batch_size;
  unsigned int rec_no = 0;
  
  std::shared_ptr<std::list<FastaSequenceData>> ret_list = std::make_shared<std::list<FastaSequenceData>>();
  
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
  
  // bool endLoop = false;
  // while (start_fptr < thread_end && start_fptr < adjusted_end)
  while (start_fptr < end_fptr)
  {
    // find start of next delimiter inside [start_fptr, end_fptr)
    const char* entryStart = std::search(start_fptr, end_fptr, bm); //delim, delim + delim_len
    if (entryStart == end_fptr) {
      // no more delimiters in this chunk
      break;
    }
    
    // find next delimiter after this one to mark end-of-entry
    const char* nextDelim = std::search(entryStart + delim_len, end_fptr, bm); //delim, delim + delim_len
    // entryEnd is one-past-last-byte of the entry region; may be end_fptr (>= start_fptr)
    const char* entryEnd = (nextDelim == end_fptr) ? end_fptr : nextDelim;
    
    // sanity: entryEnd must be after entryStart
    if (entryEnd <= entryStart) {
      // malformed or empty; advance to avoid infinite loop
      start_fptr = (entryEnd < end_fptr) ? entryEnd + 1 : end_fptr;
      continue;
    }
    
    // compute length; safe because entryEnd >= entryStart
    size_t entry_len = static_cast<size_t>(entryEnd - entryStart);
    if (entry_len == 0) {
      start_fptr = entryEnd;
      // if (entryEnd >= end_fptr) break;
      if (entryEnd >= end_fptr) break;
      continue;
    }
  
    rec_no++;
    if(rec_no < from_rec){
      // advance to next search position
      start_fptr = entryEnd;
      continue;
    }
    if(rec_no > to_rec){
      break;
    }
    // build string_view and trim trailing CR/LF without ever dereferencing entryEnd
    std::string_view sv_entry(entryStart, entry_len);
    trim(sv_entry);
    while (!sv_entry.empty() && (sv_entry.back() == '\r' || sv_entry.back() == '\n')) {
      sv_entry.remove_suffix(1);
    }
        
    bool converted = false;
    std::string utf8_str = utf16_to_utf8_iconv(sv_entry, converted);
    std::string_view to_parse = converted ? std::string_view(utf8_str) : sv_entry;
    
    // Parse and skip empty/malformed sequences
    FastaSequenceData conv_entry = CastToType(to_parse); //sv_entry
    
    conv_entry.seq.erase(std::remove_if(conv_entry.seq.begin(), conv_entry.seq.end(),
                                        [](char c){ 
                                          // allow A-Z, a-z, '*' (you may customize for DNA vs protein)
                                          return !((std::isalpha((unsigned char)c) || c=='*') && static_cast<int>(static_cast<unsigned char>(c)) > 10); 
                                        }), conv_entry.seq.end());
    
    // advance to next search position
    start_fptr = entryEnd;
    if (conv_entry.seq.empty()) {
      // if (entryEnd >= adjusted_end) break;
      if (entryEnd >= end_fptr) break;
      continue;
    }
  
    ret_list->emplace_back(conv_entry);  
    AddRecordCount();
    } // while - End
 
  return ret_list; 
}

arrow::Status ArrowWrapper::Impl::FinishOutputStream()
{
  try{
  
    if(!save2file){
      if(writer_threads.size() > 0)
      for(auto &t : writer_threads){
        if(t.get().joinable())
          t.get().join();
      }
      return arrow::Status::OK(); 
    }
    
    writer_running.store(false, std::memory_order_release);
    waiting4writer_cond.notify_all();
    Rcpp::Rcout << output_filename << std::endl << std::flush; //DEBUG
    if(writer_failed.load(std::memory_order_acquire))
    {
      return arrow::Status::Invalid(std::string("FinishOutputStream(): Writer thread(s) failed: ") + writer_error_msg);
    } 
    
    {
      std::unique_lock<std::mutex> lk(finishing_mutex);
      finishing_cond.wait(lk, [this]() {
        RcppThread::checkUserInterrupt();
        return (!this->save2file || (!writer_writing.load(std::memory_order_acquire) && writer_finishing.load(std::memory_order_acquire)) || writer_failed.load(std::memory_order_acquire));
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
      writer_running.store(true, std::memory_order_release);
      ARROW_RETURN_NOT_OK(WriteBatch2File());
      writer_running.store(false, std::memory_order_release);
    }
    
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
   omp_set_lock(&rbv_batchLock);
  #endif
    rbv_batch->clear();
    rbv_batch->shrink_to_fit();
  #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_unset_lock(&rbv_batchLock);
  #endif
    if(verbose)
      Rcpp::Rcout << "Done writing to file." << std::endl <<std::flush; 

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
    
    if(parquet_writer){
      arrow::Status st2 = parquet_writer->Close();
      if (!st2.ok()) {
        Rcpp::Rcerr << std::string("FinishOutputStream(): Error closing output writer stream (Parquet): ") + st2.message() + std::string(" : ") + st2.detail()->ToString() << std::endl << std::flush;
        return st2;
      }
      parquet_writer.reset();
    }
    
    if(parquetFileStream){
      if(!parquetFileStream->closed())
      {
        arrow::Status st2 = parquetFileStream->Close();
        if (!st2.ok()) {
          Rcpp::Rcerr << std::string("FinishOutputStream(): Error closing parquetFileStream: ") + st2.ToString() << std::endl << std::flush;
          return st2;
        }
      }
      parquetFileStream.reset();
    }
    
    if (outFileStream) {
      if(!outFileStream->closed())
      {
        arrow::Status st2 = outFileStream->Close();
        if (!st2.ok()) {
          Rcpp::Rcerr << std::string("[FinishOutputStream()]: Error closing outFileStream: ") + st2.ToString() << std::endl << std::flush;
          return st2;
        }
      }
      outFileStream.reset();
    }
    // writer_writing.store(false);
    return arrow::Status::OK();
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[FinishOutputStream()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[FinishOutputStream()]: Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[FinishOutputStream()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }catch(...){
    Rcpp::Rcerr << "[FinishOutputStream()]: Unknown Exception" << std::endl << std::flush;
  }
  return arrow::Status::Invalid("[FinishOutputStream()]: Caught an Exception.");
}

arrow::Status ArrowWrapper::FinishOutputStream()
{
  return pImpl->FinishOutputStream();
}

arrow::Status ArrowWrapper::Impl::WriteBatch2File()
{
  // DO NOT USE R API CALLS FROM WITHIN STD::(J)THREAD AS IT LEADS TO STACK CORRUPTION
  try{
    if(writer_writing.load(std::memory_order_acquire) || writer_failed.load(std::memory_order_acquire))
    {
      auto tmp_batch_size = rb_batch_size.load(std::memory_order_acquire);
      if(tmp_batch_size > 1)
        tmp_batch_size--;
      rb_batch_size.store(tmp_batch_size, std::memory_order_release);
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
  
    // mark the writer as active BEFORE spawning to avoid producers racing to increase rb_batch_size
    writer_writing.store(true, std::memory_order_release);
    
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
          RcppThread::checkUserInterrupt(); // R API calls inside threads will crash c++ runtime
          if (rb)
          {
            if (rb->num_rows() > 0)
            {
              arrow::Status rb_sts = rb->ValidateFull();
              arrow::Status sts;
              if (rb_sts.ok())
              {
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                  omp_set_lock(&rec_writerLock);
    #endif
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp critical (arrow_write_lock)
#endif
{
                  sts = rec_writer->WriteRecordBatch(*rb);
                  static_cast<void>(outFileStream->Flush());
}
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                  omp_unset_lock(&rec_writerLock);
    #endif
                  
                  if (!sts.ok())
                  {
                    {
                      std::lock_guard<std::mutex> lk(writer_error_mtx);
                      writer_error_msg = std::string("WriteBatch2File(): Error writing RB (CSV/IPC): ") + sts.message();
                    }
                    writer_failed.store(true, std::memory_order_release);
                    break;
                  }  
                
              }
            }
            rb.reset();
          }
        }
        break;
      }
    case ArrowWrapper::EOutputFormat::eParquet :{
        for(const auto &rb: rbv_buffer){
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
            omp_set_lock(&rec_writerLock);
#endif
        arrow::Status st1 = parquet_writer->WriteRecordBatch(*rb);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&rec_writerLock);
#endif
          if (!st1.ok()) {
            {
              std::lock_guard<std::mutex> lk(writer_error_mtx);
              writer_error_msg = std::string("WriteBatch2File(): Error writing table to file (Parquet): ") + st1.message() + std::string(" : ") + st1.detail()->ToString();
            }
            writer_failed.store(true, std::memory_order_release);
            Rcpp::Rcerr << writer_error_msg << std::endl << std::flush; //DEBUG
            break;
          }
          
        }
        break;
      }
    case ArrowWrapper::EOutputFormat::unknown:
    default: {
      {
        std::lock_guard<std::mutex> lk(writer_error_mtx);
        writer_error_msg = std::string("WriteBatch2File() - Unsupported output format. Supported values ipc/csv/parquet");
      }
      writer_failed.store(true, std::memory_order_release);
      break;
    }
    }
    
    rbv_buffer.clear();
    rbv_buffer.shrink_to_fit();
    
    // release local big containers
    #ifdef ARROW_HAVE_MEMORY_POOL_RELEASE
      arrow::default_memory_pool()->ReleaseUnused();
    #endif

    finishing_cond.notify_all();
      
  return arrow::Status::OK();
  }catch(const std::runtime_error &e){
    {
      std::lock_guard<std::mutex> lk(writer_error_mtx);
      writer_error_msg = std::string("WriteBatch2File(): C++ Runtime Error : ") + e.what();
    }
    writer_failed.store(true, std::memory_order_release);
  }catch(const Rcpp::exception &e){
    {
      std::lock_guard<std::mutex> lk(writer_error_mtx);
      writer_error_msg = std::string("WriteBatch2File() - Rcpp Exception : ") + e.what();
    }
    writer_failed.store(true, std::memory_order_release);
  }catch(const std::exception &e){
    {
      std::lock_guard<std::mutex> lk(writer_error_mtx);
      writer_error_msg = std::string("WriteBatch2File() - C++ Exception : ") + e.what();
    }
    writer_failed.store(true, std::memory_order_release);
  }catch(...){
    {
      std::lock_guard<std::mutex> lk(writer_error_mtx);
      writer_error_msg = "WriteBatch2File() - Unknown Exception";
    }
    writer_failed.store(true, std::memory_order_release);
  }
  writer_finishing.store(true, std::memory_order_release);
  finishing_cond.notify_all();
  return arrow::Status::OK();
}

void CountCharacter_thread(const std::string &filename, char character, std::atomic<unsigned int> &count, size_t start, size_t end)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
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

int ArrowWrapper::Impl::CountCharacter(std::string filename, char character, unsigned int num_threads)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
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
  
  unsigned int totalCount_ = totalCount.load();
  max_records = (max_records >= totalCount_) ? max_records : totalCount_;
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
    REprintf("File is empty: %s \n", filename.data());
    return -1;
  }
}

std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::Impl::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  try{
  
    std::string delim_str(delim);
    const size_t delim_len = std::strlen(delim);
    using BM = std::boyer_moore_horspool_searcher<const char*>; //std::boyer_moore_searcher<const char*>;
    auto bm = BM(delim, delim + delim_len);

  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> file_ptrs = MMapFile(filename, delim);

  char *start_of_file = std::get<1>(*file_ptrs).get();
  char *end_of_file = std::get<3>(*file_ptrs);
  
  char *p = start_of_file;

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
  
  arrow::RecordBatchVector ret_results;

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_lock_t pLock;
  omp_lock_t ret_resultsLock;
  omp_init_lock(&pLock);
  omp_init_lock(&ret_resultsLock);
#endif

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp parallel num_threads(num_threads) shared(end_of_file, start_of_file, bm, utf16_to_utf8_iconv)
#endif
  {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp for schedule(dynamic) nowait 
#endif
    for (int i = 0; i < num_threads; ++i)
    {

      size_t chunk_size = (end_of_file - start_of_file) / num_threads;
      char *thread_start = start_of_file + i * chunk_size;
      char *thread_end = (i == num_threads - 1) ? end_of_file : (thread_start + chunk_size);
      
      // Cache the length so we don't recalculate it in the loop
      size_t delim_len = strlen(delim); 
      
      // 1. SAFE BACKWARD SEARCH: Prevent underflowing past start_of_file
      if (i != 0) {
        while (thread_start > start_of_file && strncmp(thread_start, delim, delim_len) != 0) {
          --thread_start;
        }
      }
      
      // 2. SAFE FORWARD SEARCH: Prevent overflowing past end_of_file
      if (i != num_threads - 1) {
        while (thread_end < (end_of_file - delim_len) && strncmp(thread_end, delim, delim_len) != 0) {
          ++thread_end;
        }
        
        // Only step back if we actually found the delimiter (and didn't just hit EOF)
        if (thread_end < end_of_file) {
          // Depending on how your parser works, you might not even need this -1. 
          // If thread_end is pointing exactly AT the '>', that is usually the perfect boundary.
          thread_end = thread_end - 1; 
        }
      }
      // //  Get the thread-specific range to process
      // size_t chunk_size = (end_of_file - start_of_file) / num_threads;
      // char *thread_start = start_of_file + i * chunk_size;
      // char *thread_end = (i == num_threads - 1) ? end_of_file : (thread_start + chunk_size);
      // 
      // if (thread_start != start_of_file)
      // {
      //   while (strncmp(thread_start, delim, strlen(delim)) != 0)
      //   {
      //     --thread_start;
      //   }
      // }
      // if (thread_end != end_of_file)
      // {
      //   while (strncmp(thread_end, delim, strlen(delim)) != 0)
      //   {
      //     ++thread_end;
      //   }
      //   thread_end = thread_end - 1;
      // }
      
      // Process the entries within the thread's range
      char *entryStart = nullptr;
      char *entryEnd = nullptr;

      // const size_t delim_len = std::strlen(delim);
      
      // p is not used now; we create piter
      const char* piter = thread_start;
      while (piter < thread_end && piter < end_of_file)
        {
        // find start of next delimiter inside [piter, thread_end)
        const char* entryStart = std::search(piter, static_cast<const char*>(thread_end), bm); //delim, delim + delim_len
        if (entryStart == thread_end) {
          // no more delimiters in this chunk
          break;
        }
        
        // find next delimiter after this one to mark end-of-entry
        const char* nextDelim = std::search(entryStart + delim_len, static_cast<const char*>(thread_end), bm); 
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
          if (entryEnd >= end_of_file) break;
          continue;
        }
        
        // build string_view and trim trailing CR/LF without ever dereferencing entryEnd
        std::string_view sv_entry(entryStart, entry_len);
        trim(sv_entry);
        while (!sv_entry.empty() && (sv_entry.back() == '\r' || sv_entry.back() == '\n')) {
          sv_entry.remove_suffix(1);
        }
        
        bool converted = false;
        std::string utf8_str = utf16_to_utf8_iconv(sv_entry, converted);
        std::string_view to_parse = converted ? std::string_view(utf8_str) : sv_entry;
        
        // Parse and skip empty/malformed sequences
        FastaSequenceData conv_entry = CastToType(to_parse); 
        
        conv_entry.seq.erase(std::remove_if(conv_entry.seq.begin(), conv_entry.seq.end(),
                                      [](char c){ 
                                        // allow A-Z, a-z, '*' (you may customize for DNA vs protein)
                                        return !((std::isalpha((unsigned char)c) || c=='*') && static_cast<int>(static_cast<unsigned char>(c)) > 10); 
                                      }), conv_entry.seq.end());
        
        if (conv_entry.seq.empty()) {
          piter = entryEnd;
          if (entryEnd >= end_of_file) break;
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
          tmp_result->clear();
          tmp_result->shrink_to_fit();
        }
        
        // advance to next search position
        piter = entryEnd;
        // if we hit the real file-end, we're done
        if (entryEnd >= end_of_file) {
          break;
        }
      }
      
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
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[SplitFilesIntoEntries()] - C++ Runtime Exception : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[SplitFilesIntoEntries()] - Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch (const std::exception &e) {
    Rcpp::Rcerr << std::string("[SplitFilesIntoEntries()] - C++ exception : ") + e.what() << std::endl << std::flush;
  }catch (...) {
    Rcpp::Rcerr << "[SplitFilesIntoEntries()]: Unknown exception" << std::endl << std::flush;
  }
  return std::make_shared<arrow::RecordBatchVector>();
}
std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return pImpl->SplitFilesIntoEntries(filename, delim, num_threads, Entry_callback, return_values);
}