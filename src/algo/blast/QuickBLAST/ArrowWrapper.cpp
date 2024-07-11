#if defined(WIN32) || defined(MINGW32)
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

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

ArrowWrapper::ArrowWrapper()
{

  ok_promise.set_value(arrow::Status::OK());

  outputStream = std::make_shared<std::ostringstream>();

  arrow_LFS = arrow::fs::LocalFileSystem();
  std::string username;
#if defined(linux) || defined(MINGW32)
  username = getlogin();
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

  ipc_options.allow_64bit = false;
  ipc_options.use_threads = true;
  ipc_options.metadata_version = arrow::ipc::MetadataVersion::V5;
  ipc_options.codec = arrow::util::Codec::Create(arrow::Compression::LZ4_FRAME).ValueOrDie();

  ipc_options.write_legacy_ipc_format = false;

  rbv_batch = std::make_shared<arrow::RecordBatchVector>();
  blast_metadata = std::make_shared<arrow::KeyValueMetadata>();
  AddFASTAMetadata("format", "Arrow Parquet");
  AddFASTAMetadata("Created By", username);
  fasta_schema = arrow::schema({arrow::field("index", arrow::int32()), arrow::field("header", arrow::utf8()), arrow::field("sequence", arrow::utf8())});

  seq_info_type = arrow::struct_({
      arrow::field("num_alignments", arrow::int8()),
      arrow::field("seqids", arrow::struct_({arrow::field("qseqid", arrow::utf8()),
                                             arrow::field("sseqid", arrow::utf8())})),
      arrow::field("seqs", arrow::struct_({arrow::field("qseq", arrow::large_utf8()),
                                           arrow::field("sseq", arrow::large_utf8())})),
      arrow::field("strands", arrow::utf8()),
      arrow::field("lengths", arrow::struct_({arrow::field("qlen", arrow::int8()),
                                              arrow::field("slen", arrow::int8())})),
  });
  this->hsp_type = arrow::struct_({arrow::field("pident", arrow::float64()),
                                   arrow::field("pident_gap", arrow::float64()),
                                   arrow::field("frames", arrow::utf8()),
                                   arrow::field("evalue", arrow::float64()),
                                   arrow::field("length", arrow::int8()),
                                   arrow::field("length01", arrow::float64()),
                                   arrow::field("qstart", arrow::int8()),
                                   arrow::field("qend", arrow::int8()),
                                   arrow::field("sstart", arrow::int8()),
                                   arrow::field("send", arrow::int8()),
                                   arrow::field("bitscore", arrow::float64()),
                                   arrow::field("score", arrow::float64()),
                                   arrow::field("qcovhsp", arrow::float64()),
                                   arrow::field("blast_score", arrow::float64()),
                                   arrow::field("gaps", arrow::int8()),
                                   arrow::field("nident", arrow::int8()),
                                   arrow::field("mismatch", arrow::int8()),
                                   arrow::field("positive", arrow::int8()),
                                   arrow::field("n_splices", arrow::int8()),
                                   arrow::field("hsp_num", arrow::int8()),
                                   arrow::field("sum_evalue", arrow::float64()),
                                   arrow::field("product_coverage", arrow::float64()),
                                   arrow::field("overall_identity", arrow::float64()),
                                   arrow::field("negative_count", arrow::int8()),
                                   arrow::field("matches", arrow::float64()),
                                   arrow::field("high_quality_percent_coverage", arrow::float64()),
                                   arrow::field("exon_identity", arrow::float64()),
                                   arrow::field("consensus_splices", arrow::float64()),
                                   arrow::field("comp_adj_method", arrow::float64())});
  this->alignment_scores_type = arrow::list({hsp_type});

  blast_schema = arrow::schema({arrow::field("seq_info", seq_info_type),
                                arrow::field("hsps", hsp_type)});
#if defined(_OPENMP) || defined(WIN32)
  omp_init_lock(&rec_countLock);
  omp_init_lock(&rbv_batchLock);
  omp_init_lock(&rec_writerLock);
#endif
}

ArrowWrapper::~ArrowWrapper()
{
#if defined(_OPENMP) || defined(WIN32)
  omp_destroy_lock(&rec_countLock);
  omp_destroy_lock(&rbv_batchLock);
  omp_destroy_lock(&rec_writerLock);
#endif

  std::cout << "~ArrowWrapper " << std::endl;
}

void ArrowWrapper::SetBatchSize(int batch_size)
{
  assert(batch_size > 0);
  this->rb_batch_size = batch_size;
}

int ArrowWrapper::GetRecordCount() { return rec_count; }

void ArrowWrapper::ResetRecordCount()
{
  rec_count = 0;
}
void ArrowWrapper::AddRecordCount()
{
#if defined(_OPENMP) || defined(WIN32)
  omp_set_lock(&rec_countLock);
#endif
  rec_count++;
#if defined(_OPENMP) || defined(WIN32)
  omp_unset_lock(&rec_countLock);
#endif
}
void ArrowWrapper::SetThreadCount(int num_threads)
{
  if (num_threads > 1)
  {
    ipc_options.use_threads = true;
#if defined(_OPENMP) || defined(WIN32)
    omp_set_num_threads(num_threads);
#endif
  }
  else
  {
    ipc_options.use_threads = false;
#if defined(_OPENMP) || defined(WIN32)
    omp_set_num_threads(1);
#endif
  }
  n_threads = num_threads;
}
arrow::Result<std::shared_ptr<arrow::RecordBatch>> ArrowWrapper::AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_)
{
  arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
  if (rb_)
  {
#if defined(_OPENMP) || defined(WIN32)
    omp_set_lock(&rbv_batchLock);
#endif
    rbv_batch->emplace_back(std::move(rb_));
    unsigned int ret_size = rbv_batch->size();
#if defined(_OPENMP) || defined(WIN32)
    omp_unset_lock(&rbv_batchLock);
#endif
    if (ret_size >= rb_batch_size)
    {
      if (writer_threads.size() >= n_threads)
      {
        static_cast<void>(writer_threads[writer_threads.size() - 1].join());
        writer_threads.erase(writer_threads.begin() + writer_threads.size() - 1);
      }
      std::thread write_thread([this]()
                               { static_cast<void>(this->WriteBatch2File()); });
      // write_thread.detach();
      static_cast<void>(writer_threads.emplace_back(std::move(write_thread)));
    }
    return arrow::Result<std::shared_ptr<arrow::RecordBatch>>(rb_);
  }
  return arrow::Result<std::shared_ptr<arrow::RecordBatch>>(error_sts);
}
arrow::Result<std::shared_ptr<arrow::RecordBatchVector>> ArrowWrapper::AddRBV2Batch(const arrow::RecordBatchVector &rbv_)
{

  arrow::Status error_sts(arrow::StatusCode::Invalid, "Error Writing to File!");
  if (rbv_.size() > 0)
  {
#if defined(_OPENMP) || defined(WIN32)
    omp_set_lock(&rbv_batchLock);
#endif
    rbv_batch->insert(rbv_batch->end(), rbv_.begin(), rbv_.end());
    unsigned int ret_size = rbv_batch->size();
#if defined(_OPENMP) || defined(WIN32)
    omp_unset_lock(&rbv_batchLock);
#endif
    if (ret_size >= rb_batch_size)
    {
      if (writer_threads.size() >= n_threads)
      {
        static_cast<void>(writer_threads[writer_threads.size() - 1].join());
        writer_threads.erase(writer_threads.begin() + writer_threads.size() - 1);
      }
      std::thread write_thread([this]()
                               { static_cast<void>(this->WriteBatch2File()); });
      // write_thread.detach();
      static_cast<void>(writer_threads.emplace_back(std::move(write_thread)));
    }
    return arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(std::make_shared<arrow::RecordBatchVector>(rbv_));
  }
  //}
  return arrow::Result<std::shared_ptr<arrow::RecordBatchVector>>(error_sts);
}
arrow::Status ArrowWrapper::CreateOutputStream(const std::string &outFile)
{
  outFileStream = arrow_LFS.OpenAppendStream(outFile, GetBLASTMetadata()).ValueOrDie();
  auto codec = arrow::util::Codec::Create(arrow::Compression::GZIP).ValueOrDie();

  compressed_outstream = arrow::io::CompressedOutputStream::Make(codec.get(), outFileStream).ValueOrDie();

  return arrow::Status::OK();
}

std::shared_ptr<arrow::DataType> ArrowWrapper::GetSeqInfoType(void)
{
  return seq_info_type;
}
std::shared_ptr<arrow::DataType> ArrowWrapper::GetAlignmentScoresType(void)
{
  return alignment_scores_type;
}
std::shared_ptr<arrow::DataType> ArrowWrapper::GetHSPType(void)
{
  return hsp_type;
}
std::shared_ptr<arrow::Schema> ArrowWrapper::GetBLASTSchema(void)
{
  return blast_schema;
}
std::shared_ptr<arrow::Schema> ArrowWrapper::GetFASTASchema(void)
{
  return fasta_schema;
}
/*std::shared_ptr<parquet::WriterProperties> ArrowWrapper::GetParquetWriterProps(void)
{
    return parquet_writer_props;
}*/
/*std::shared_ptr<parquet::ArrowWriterProperties> ArrowWrapper::GetArrowWriterProps(void)
{
    return arrow_writer_props;
}*/
std::shared_ptr<arrow::KeyValueMetadata> ArrowWrapper::GetBLASTMetadata(void)
{
  return blast_metadata;
}
void ArrowWrapper::AddFASTAMetadata(const std::string &key, const std::string &value)
{
  blast_metadata->Append(key, value);
}
arrow::ipc::IpcWriteOptions ArrowWrapper::GetArrowIPCOptions(void)
{
  return ipc_options;
}

std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> ArrowWrapper::MMapFile(const std::string_view &filename, const char *delim)
{
  FILE *file_ptr;
  long fileSize;
  char *end_of_file_ptr;
  std::string delim_str(delim);

  file_ptr = fopen(filename.data(), "r");

  if (!file_ptr)
  {
    std::cerr << "Error: Failed to open file: " << filename.data() << std::endl;
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
    std::cerr << "Error: Failed to map file " << filename.data() << std::endl;
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

long ArrowWrapper::GetFileSize(FILE *file_ptr)
{
  fseek(file_ptr, 0, SEEK_END);
  long fileSize = ftell(file_ptr);
  rewind(file_ptr);
  return fileSize;
}

void ArrowWrapper::FinishOutputStream()
{
  if (rbv_batch->size() > 0)
  {
    static_cast<void>(WriteBatch2File());
  }

  std::thread fin_thread([this]()
                         {
                           if (this->writer_threads.size() > 0)
                           {
                             static_cast<void>(this->writer_threads[this->writer_threads.size() - 1].join());
                             //  for (auto &wrt_thread : this->writer_threads)
                             //  {

                             //    static_cast<void>(wrt_thread.join());
                             //  }
                           }
                           // static_cast<void>(this->rec_writer->Close());
                           this->writer_threads.clear();
                           this->rbv_batch->clear();
                           // Rcpp::Rcout << "Done writing to file." << std::endl;
                         });

  fin_thread.detach();
}

arrow::Status ArrowWrapper::WriteBatch2File()
{
  assert(compressed_outstream);

  auto writer_ = arrow::ipc::MakeFileWriter(compressed_outstream.get(), GetBLASTSchema(), GetArrowIPCOptions(), GetBLASTMetadata());
  if (!writer_.ok())
  {
    return writer_.status();
  }
  std::shared_ptr<arrow::ipc::RecordBatchWriter> rec_writer = writer_.ValueOrDie();

  arrow::RecordBatchVector rbv_buffer;
#if defined(_OPENMP) || defined(WIN32)
  omp_set_lock(&rbv_batchLock);
#endif

  rbv_buffer.swap(*rbv_batch);
#if defined(_OPENMP) || defined(WIN32)
  omp_unset_lock(&rbv_batchLock);
#endif

  for (const auto &rb : rbv_buffer)
  {
    if (rb)
    {
      if (rb->num_rows() > 0)
      {
        arrow::Status rb_sts = rb->ValidateFull();
        if (rb_sts.ok())
        {
#if defined(_OPENMP) || defined(WIN32)
          omp_set_lock(&rec_writerLock);
#endif
          arrow::Status sts = rec_writer->WriteRecordBatch(*rb);
#if defined(_OPENMP) || defined(WIN32)
          omp_unset_lock(&rec_writerLock);
#endif

          if (!sts.ok())
          {
            std::cout << "ERROR : Could not write RB" << sts.detail() << std::endl
                      << sts.message() << std::endl
                      << rb->schema()->ToString() << std::endl;
            return sts;
          }
        }
        else
        {
          std::cerr << "Warn : Invalid Alignment RB (Not Writing) : " << rb_sts.detail() << std::endl
                    << rb_sts.message() << std::endl
                    << rb->schema()->ToString() << std::endl;
        }
      }
    }
    else
    {
      std::cerr << "ERR : Invalid Alignment RB Ptr (Not Writing)..." << std::endl;
    }
  }

  static_cast<void>(rec_writer->Close());

  return arrow::Status::OK();
}

void CountCharacter_thread(const std::string &filename, char character, std::atomic<int> &count, size_t start, size_t end)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    std::cerr << "Failed to open the file." << std::endl;
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

int ArrowWrapper::CountCharacter(std::string filename, char character, int num_threads)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    std::cerr << "Failed to open the file." << std::endl;
    return 1;
  }

  file.seekg(0, std::ios::end);
  size_t fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  std::atomic<int> totalCount(0);
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

  return totalCount;
}

int ArrowWrapper::GetColumnCount(const std::string_view &filename, char delim)
{
  std::ifstream file(filename.data());
  if (!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
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
    std::cerr << "File is empty: " << filename << std::endl;
    return -1;
  }
}
