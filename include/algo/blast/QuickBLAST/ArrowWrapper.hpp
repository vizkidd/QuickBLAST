#pragma once
#include <algo/blast/QuickBLAST/commons.hpp>

#include <iostream>
#include <regex>
#include <string>
#include <future>
#include <memory>
#include <arrow/api.h>
#include <arrow/ipc/options.h>
// #include <parquet/properties.h>
#include <arrow/util/type_fwd.h>
#include <arrow/builder.h>
#include <arrow/record_batch.h>
// #include <arrow/util/string.h>

#include <arrow/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/ipc/writer.h>
#include <arrow/ipc/options.h>
#include <arrow/io/api.h>
#include <arrow/filesystem/filesystem.h>
// #include <parquet/arrow/reader.h>
// #include <parquet/arrow/writer.h>
// #include <parquet/properties.h>
// #include <boost/lexical_cast.hpp>

#include <arrow/csv/api.h>
#include <parquet/api/writer.h>
#include <parquet/arrow/writer.h>
#include <arrow/util/type_fwd.h>

using parquet::WriterProperties;
using parquet::ArrowWriterProperties;
using parquet::ParquetVersion;
using parquet::ParquetDataPageVersion;
using arrow::Compression;

#ifndef ARROWWRAPPER_HPP
#define ARROWWRAPPER_HPP

class ArrowWrapper
{
private:
  struct Impl;
  std::unique_ptr<Impl> pImpl;
  
public:
  enum EOutputFormat{
    eIPC = 0,
    eCSV = 1,
    eParquet = 2,
    unknown = 3
  };
  
  ~ArrowWrapper();
  ArrowWrapper();
  void SetBatchSize(unsigned int batch_size);
  void SetVerbosity(bool verbose);
  unsigned int GetBatchSize(void);
  arrow::Status FinishOutputStream();
  int CountCharacter(std::string filename, char character, int num_threads);
  // template <typename T1>
  // std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values = false);
  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> MMapFile(const std::string_view &filename, const char *delim);
  std::shared_ptr<std::list<FastaSequenceData>> FetchRecordByBatch(const std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> &file_ptr, unsigned int batch_size, unsigned int from_rec, const char *delim);
  long GetFileSize(FILE *file_ptr);
  std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values = false);
  FastaSequenceData CastToType(const std::string_view &full_entry_sv);
  // template <typename T>
  // T CastToType(const std::string &full_entry);
  unsigned int GetRecordCount();
  void ResetRecordCount();
  void AddRecordCount();
  unsigned int GetProcRecordCount();
  void ResetProcRecordCount();
  void AddProcRecordCount();
  unsigned int GetPendingRecordCount();
  void SetBLASTSeqLimit(unsigned int);
  unsigned int GetBLASTSeqLimit();
  void SetThreadCount(int num_threads);
  arrow::Status AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_);
  // arrow::Status AddRBV2Batch(arrow::RecordBatchVector &rbv_);
  arrow::Status CreateOutputStream(std::string &outFile, const std::string& outputFormat);
  
  std::string GetOutputFormat(void);
  std::shared_ptr<arrow::DataType> GetSeqInfoType(void);
  std::shared_ptr<arrow::DataType> GetAlignmentScoresType(void);
  std::shared_ptr<arrow::DataType> GetHSPType(void);
  std::shared_ptr<arrow::Schema> GetBLASTSchema(void);
  std::shared_ptr<arrow::Schema> GetFASTASchema(void);
  std::shared_ptr<arrow::Schema> GetSchema(void);
  // std::shared_ptr<parquet::WriterProperties> GetParquetWriterProps(void);
  // std::shared_ptr<parquet::ArrowWriterProperties> GetArrowWriterProps(void);
  std::shared_ptr<arrow::KeyValueMetadata> GetBLASTMetadata(void);
  void AddFASTAMetadata(const std::string &key, const std::string &value);
  arrow::ipc::IpcWriteOptions GetArrowIPCOptions(void);
  arrow::csv::WriteOptions GetArrowCSVOptions(void);
};

struct ArrowWrapper::Impl
{
  ~Impl();
  Impl();
  std::shared_ptr<arrow::Schema> fasta_schema, blast_schema;
  std::shared_ptr<arrow::DataType> alignment_scores_type, seq_info_type, hsp_type;
  std::shared_ptr<arrow::KeyValueMetadata> blast_metadata;
  arrow::fs::LocalFileSystem arrow_LFS;
  std::shared_ptr<arrow::io::OutputStream> outFileStream;
  std::shared_ptr<arrow::io::FileOutputStream> parquetFileStream;
  std::shared_ptr<arrow::ipc::RecordBatchWriter> rec_writer;
  std::unique_ptr<parquet::arrow::FileWriter> parquet_writer;
  std::shared_ptr<arrow::RecordBatchVector> rbv_batch;
  std::vector<safe_jthread> writer_threads;
  std::thread finisher_thread;
  std::string output_filename, output_format;
  bool save2file;
  bool verbose;
  
  unsigned int parquet_batch_size = 1024, rec_count = 1, blast_sequence_limit = 0, proc_rec_count = 0, n_threads = 1, max_records = 1024; 
  std::atomic<unsigned int> itr_add{1};
  std::atomic<unsigned int> itr_mul{1};
  std::atomic<unsigned int> rb_batch_size{1024};
  
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_lock_t rec_countLock;
  omp_lock_t proc_rec_countLock;
  omp_lock_t writer_threadsLock;
  omp_lock_t rec_writerLock;
  omp_lock_t rbv_batchLock;
#endif
  
  std::atomic<bool> writer_running{false};     // signal to stop/wakeup waits
  std::atomic<bool> writer_writing{false};     // signal to stop/wakeup waits
  std::atomic<bool> writer_waiting{false};     // signal to stop/wakeup waits
  std::atomic<bool> writer_finishing{false};     // signal to stop/wakeup waits
  std::atomic<bool> writer_failed{false};
  std::string writer_error_msg;
  std::mutex writer_error_mtx;
  std::mutex writer_threads_mutex;
  std::condition_variable finishing_cond; 
  std::mutex finishing_mutex;
  std::condition_variable waiting4writer_cond;
  std::mutex waiting4writer_mutex;
  arrow::ipc::IpcWriteOptions ipc_options;
  arrow::csv::WriteOptions csv_options;
  std::shared_ptr<parquet::WriterProperties> parquet_props;
  std::shared_ptr<parquet::ArrowWriterProperties> parquet_arrow_props;
  // std::shared_ptr<parquet::WriterProperties> parquet_writer_props;
  // std::shared_ptr<parquet::ArrowWriterProperties> arrow_writer_props;
  // parquet::WriterProperties::Builder props_bldr;
  // parquet::ArrowWriterProperties::Builder arrow_props_bldr;
  // std::shared_ptr<std::ostringstream> outputStream;
  // Rcpp::XPtr<std::ostringstream> outputStream;
  
  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> MMapFile(const std::string_view &filename, const char *delim);
  std::shared_ptr<std::list<FastaSequenceData>> FetchRecordByBatch(const std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> &file_ptr, unsigned int batch_size, unsigned int from_rec, const char *delim);
  long GetFileSize(FILE *file_ptr);
  
  void SetBatchSize(unsigned int batch_size);
  void SetVerbosity(bool verbose);
  unsigned int GetBatchSize(void);
  arrow::Status FinishOutputStream();
  arrow::Status WriteBatch2File();
  void writer_loop(void);
  int GetColumnCount(const std::string_view &filename, char delim = '\t');
  int CountCharacter(std::string filename, char character, unsigned int num_threads);
  // template <typename T1>
  // std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values = false);
  std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values = false);
  // template <typename T>
  // std::string CastToType(const std::string &full_entry);
  // FastaSequenceData CastToType(const std::string &full_entry);
  FastaSequenceData CastToType(const std::string_view &full_entry_sv);
  // T CastToType(const std::string &full_entry);
  unsigned int GetRecordCount();
  void ResetRecordCount();
  void AddRecordCount();
  unsigned int GetProcRecordCount();
  void ResetProcRecordCount();
  void AddProcRecordCount();
  unsigned int GetPendingRecordCount();
  void SetBLASTSeqLimit(unsigned int);
  unsigned int GetBLASTSeqLimit();
  void SetThreadCount(int num_threads);
  arrow::Status AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_);
  // arrow::Status AddRBV2Batch(arrow::RecordBatchVector &rbv_);
  arrow::Status CreateOutputStream(std::string &outFile, const std::string& outputFormat);
  ArrowWrapper::EOutputFormat OutputFormat2Enum(const std::string& str);
  std::string GetOutputFormat(void);
  std::shared_ptr<arrow::DataType> GetSeqInfoType(void);
  std::shared_ptr<arrow::DataType> GetAlignmentScoresType(void);
  std::shared_ptr<arrow::DataType> GetHSPType(void);
  std::shared_ptr<arrow::Schema> GetBLASTSchema(void);
  std::shared_ptr<arrow::Schema> GetFASTASchema(void);
  std::shared_ptr<arrow::Schema> GetSchema(void);
  // std::shared_ptr<parquet::WriterProperties> GetParquetWriterProps(void);
  // std::shared_ptr<parquet::ArrowWriterProperties> GetArrowWriterProps(void);
  std::shared_ptr<arrow::KeyValueMetadata> GetBLASTMetadata(void);
  void AddFASTAMetadata(const std::string &key, const std::string &value);
  arrow::ipc::IpcWriteOptions GetArrowIPCOptions(void);
  arrow::csv::WriteOptions GetArrowCSVOptions(void);
};

#endif // ARROWWRAPPER_HPP