#ifdef _OPENMP
#include "omp.h"
#endif
#include <arrow/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/ipc/writer.h>
#include <arrow/ipc/options.h>
#include <arrow/io/api.h>
#include <arrow/filesystem/filesystem.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <boost/lexical_cast.hpp>
#include <tuple>
#include <sys/mman.h>
#include <string_view>
#include <map>
#include <regex>
#include <ncbi-tools++/algo/blast/api/blast_types.hpp>
#include <ncbi-tools++/objects/seqalign/Seq_align_set.hpp>
#include <cassert>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unistd.h>
#include <parquet/properties.h>
#include <arrow/util/type_fwd.h>
#include <future>
#include <ncbi-tools++/algo/blast/api/blast_options_handle.hpp>
#include <ncbi-tools++/algo/blast/format/blastfmtutil.hpp>
#include <ncbi-tools++/objtools/align_format/align_format_util.hpp>
#include <ncbi-tools++/objects/seqset/Seq_entry.hpp>
#include <ncbi-tools++/serial/iterator.hpp>
#include <ncbi-tools++/algo/blast/api/bl2seq.hpp>
#include <ncbi-tools++/objects/seq/Seq_literal.hpp>
#include <ncbi-tools++/objtools/simple/simple_om.hpp>
#include <ncbi-tools++/algo/blast/api/sseqloc.hpp>
#include <ncbi-tools++/objmgr/seq_vector.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

struct FastaSequenceData
{
  int rec_no = 1;
  std::string header;
  std::string seq;
};
struct BLASTHitData
{
  int rec_no = 1;
  std::vector<std::string_view> col_names();
  std::vector<std::string_view> col_values();
};

class ArrowWrapper
{
private:
  std::shared_ptr<arrow::Schema> fasta_schema, blast_schema;
  std::shared_ptr<arrow::DataType> alignment_scores_type, seq_info_type, hsp_type;
  std::shared_ptr<arrow::KeyValueMetadata> blast_metadata;
  std::promise<arrow::Status> ok_promise;
  arrow::fs::LocalFileSystem arrow_LFS;
  std::shared_ptr<arrow::io::OutputStream> outFileStream;
  std::shared_ptr<arrow::ipc::RecordBatchWriter> rec_writer;
  std::shared_ptr<arrow::RecordBatchVector> rbv_batch;

  int rb_batch_size = 1024, rec_count = 0;

#ifdef _OPENMP
  omp_lock_t rec_countLock;
  omp_lock_t rbv_batchLock;
#endif

  arrow::ipc::IpcWriteOptions ipc_options;
  std::shared_ptr<parquet::WriterProperties> parquet_writer_props;
  std::shared_ptr<parquet::ArrowWriterProperties> arrow_writer_props;
  std::shared_ptr<std::ostringstream> outputStream;

  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> MMapFile(const std::string_view &filename, const char *delim);
  void CloseFilePtrs(std::tuple<FILE *, char *, long, char *> &file_ptrs);
  long GetFileSize(FILE *file_ptr);

public:
  ~ArrowWrapper();
  ArrowWrapper();
  void SetBatchSize(int batch_size)
  {
    assert(batch_size > 0);
    this->rb_batch_size = batch_size;
  }
  void FinishOutputStream();
  arrow::Status WriteBatch2File();
  std::shared_ptr<arrow::RecordBatch> ReadRecordBatchVector(const std::string &file);
  int GetColumnCount(const std::string_view &filename, char delim = '\t');
  int CountCharacter(std::string filename, char character, int num_threads);
  template <typename T1>
  std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback);
  template <typename T>
  T CastToType(const std::string &full_entry);
  int GetRecordCount() { return rec_count; }
  void ResetRecordCount()
  {
    rec_count = 0;
  }
  void AddRecordCount()
  {
#ifdef _OPENMP
    omp_set_lock(&rec_countLock);
#endif
    rec_count++;
#ifdef _OPENMP
    omp_unset_lock(&rec_countLock);
#endif
  }
  void SetThreadCount(int num_threads)
  {
    if (num_threads > 1)
    {
      ipc_options.use_threads = true;
#ifdef _OPENMP
      omp_set_num_threads(num_threads);
#endif
    }
    else
    {
      ipc_options.use_threads = false;
#ifdef _OPENMP
      omp_set_num_threads(1);
#endif
    }
  }
  arrow::Result<std::shared_ptr<arrow::RecordBatch>> AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
  {
    arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
    if (rb_)
    {
#ifdef _OPENMP
      omp_set_lock(&rbv_batchLock);
#endif
      rbv_batch->emplace_back(std::move(rb_));
      int ret_size = rbv_batch->size();
#ifdef _OPENMP
      omp_unset_lock(&rbv_batchLock);
#endif
      if (ret_size >= rb_batch_size)
      {
        std::thread write_thread([this]()
                                 { this->WriteBatch2File(); });
        write_thread.detach();
      }
      return arrow::Result<std::shared_ptr<arrow::RecordBatch>>(rb_);
    }
    return arrow::Result<std::shared_ptr<arrow::RecordBatch>>(error_sts);
  }
  arrow::Result<std::shared_ptr<arrow::RecordBatchVector>> AddRBV2Batch(const arrow::RecordBatchVector &rbv_)
  {

    arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
    if (rbv_.size() > 0)
    {
#ifdef _OPENMP
      omp_set_lock(&rbv_batchLock);
#endif
      rbv_batch->insert(rbv_batch->end(), rbv_.begin(), rbv_.end());
      int ret_size = rbv_batch->size();
#ifdef _OPENMP
      omp_unset_lock(&rbv_batchLock);
#endif
      if (ret_size >= rb_batch_size)
      {
        std::thread write_thread([this]()
                                 { this->WriteBatch2File(); });
        write_thread.detach();
      }
      return arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(std::make_shared<arrow::RecordBatchVector>(rbv_));
    }
    //}
    return arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(error_sts);
  }
  arrow::Status CreateOutputStream(const std::string &outFile)
  {
    outFileStream = arrow_LFS.OpenAppendStream(outFile, GetBLASTMetadata()).ValueOrDie();
    auto writer_ = arrow::ipc::MakeFileWriter(outFileStream.get(), GetBLASTSchema(), GetArrowIPCOptions(), GetBLASTMetadata());
    if (!writer_.ok())
    {
      return writer_.status();
    }
    rec_writer = writer_.ValueOrDie();
    return arrow::Status::OK();
  }

  std::shared_ptr<arrow::DataType> GetSeqInfoType(void)
  {
    return seq_info_type;
  }
  std::shared_ptr<arrow::DataType> GetAlignmentScoresType(void)
  {
    return alignment_scores_type;
  }
  std::shared_ptr<arrow::DataType> GetHSPType(void)
  {
    return hsp_type;
  }
  std::shared_ptr<arrow::Schema> GetBLASTSchema(void)
  {
    return blast_schema;
  }
  std::shared_ptr<arrow::Schema> GetFASTASchema(void)
  {
    return fasta_schema;
  }
  std::shared_ptr<parquet::WriterProperties> GetParquetWriterProps(void)
  {
    return parquet_writer_props;
  }
  std::shared_ptr<parquet::ArrowWriterProperties> GetArrowWriterProps(void)
  {
    return arrow_writer_props;
  }
  std::shared_ptr<arrow::KeyValueMetadata> GetBLASTMetadata(void)
  {
    return blast_metadata;
  }
  void AddFASTAMetadata(const std::string &key, const std::string &value)
  {
    blast_metadata->Append(key, value);
  }
  arrow::ipc::IpcWriteOptions GetArrowIPCOptions(void)
  {
    return ipc_options;
  }
};
