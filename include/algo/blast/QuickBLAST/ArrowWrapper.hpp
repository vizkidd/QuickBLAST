#pragma once

#include <iostream>
#include <regex>
#include <string>
#include <future>
#include <arrow/api.h>
#include <arrow/ipc/options.h>
// #include <parquet/properties.h>
#include <arrow/util/type_fwd.h>
#include <arrow/builder.h>
#include <arrow/record_batch.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include <arrow/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/ipc/writer.h>
#include <arrow/ipc/options.h>
#include <arrow/io/api.h>
#include <arrow/filesystem/filesystem.h>
// #include <parquet/arrow/reader.h>
// #include <parquet/arrow/writer.h>
// #include <parquet/properties.h>
#include <arrow/util/type_fwd.h>
// #include <boost/lexical_cast.hpp>

class QBLIBRARY_API ArrowWrapper
{
private:
    std::shared_ptr<arrow::Schema> fasta_schema, blast_schema;
    std::shared_ptr<arrow::DataType> alignment_scores_type, seq_info_type, hsp_type;
    std::shared_ptr<arrow::KeyValueMetadata> blast_metadata;
    std::promise<arrow::Status> ok_promise;
    arrow::fs::LocalFileSystem arrow_LFS;
    std::shared_ptr<arrow::io::OutputStream> outFileStream;
    std::shared_ptr<arrow::io::CompressedOutputStream> compressed_outstream;
    // std::shared_ptr<arrow::ipc::RecordBatchWriter> rec_writer;
    std::shared_ptr<arrow::RecordBatchVector> rbv_batch;
    std::vector<std::thread> writer_threads;

    unsigned int rb_batch_size = 1024, rec_count = 0, n_threads = 1;

#if defined(_OPENMP) || defined(WIN32)
    omp_lock_t rec_countLock;
    omp_lock_t rec_writerLock;
    omp_lock_t rbv_batchLock;
#endif

    arrow::ipc::IpcWriteOptions ipc_options;
    // std::shared_ptr<parquet::WriterProperties> parquet_writer_props;
    // std::shared_ptr<parquet::ArrowWriterProperties> arrow_writer_props;
    // parquet::WriterProperties::Builder props_bldr;
    // parquet::ArrowWriterProperties::Builder arrow_props_bldr;
    std::shared_ptr<std::ostringstream> outputStream;

    std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> MMapFile(const std::string_view &filename, const char *delim);
    void CloseFilePtrs(std::tuple<FILE *, char *, long, char *> &file_ptrs);
    long GetFileSize(FILE *file_ptr);

public:
    ~ArrowWrapper();
    ArrowWrapper();
    void SetBatchSize(int batch_size);
    void FinishOutputStream();
    arrow::Status WriteBatch2File();
    int GetColumnCount(const std::string_view &filename, char delim = '\t');
    int CountCharacter(std::string filename, char character, int num_threads);
    template <typename T1>
    std::shared_ptr<arrow::RecordBatchVector> SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values = false);
    template <typename T>
    T CastToType(const std::string &full_entry);
    int GetRecordCount();
    void ResetRecordCount();
    void AddRecordCount();
    void SetThreadCount(int num_threads);
    arrow::Result<std::shared_ptr<arrow::RecordBatch>> AddRB2Batch(std::shared_ptr<arrow::RecordBatch> rb_);
    arrow::Result<std::shared_ptr<arrow::RecordBatchVector>> AddRBV2Batch(const arrow::RecordBatchVector &rbv_);
    arrow::Status CreateOutputStream(const std::string &outFile);

    std::shared_ptr<arrow::DataType> GetSeqInfoType(void);
    std::shared_ptr<arrow::DataType> GetAlignmentScoresType(void);
    std::shared_ptr<arrow::DataType> GetHSPType(void);
    std::shared_ptr<arrow::Schema> GetBLASTSchema(void);
    std::shared_ptr<arrow::Schema> GetFASTASchema(void);
    std::shared_ptr<arrow::Schema> GetSchema() { return blast_schema; };
    // std::shared_ptr<parquet::WriterProperties> GetParquetWriterProps(void);
    // std::shared_ptr<parquet::ArrowWriterProperties> GetArrowWriterProps(void);
    std::shared_ptr<arrow::KeyValueMetadata> GetBLASTMetadata(void);
    void AddFASTAMetadata(const std::string &key, const std::string &value);
    arrow::ipc::IpcWriteOptions GetArrowIPCOptions(void);
};

template <typename T>
T ArrowWrapper::CastToType(const std::string &full_entry)
{

    if constexpr (std::is_same_v<T, std::string>)
    {
        return full_entry;
    }
    else if constexpr (std::is_same_v<T, FastaSequenceData>)
    {
        std::regex pattern("^>(\\w+)[\\r\\n]+([\\w\\W]+)$");
        std::regex punct_pattern("[[:punct:][:space:]\\n\\t]");
        std::smatch match;

        std::string full_entry_str(static_cast<std::string>(full_entry));
        FastaSequenceData fasta_data;
        fasta_data.rec_no = GetRecordCount();
        if (std::regex_match(full_entry, match, pattern))
        {
            fasta_data.header = match[1].str();
            fasta_data.seq = std::regex_replace(match[2].str(), punct_pattern, "");
        }
        else
        {

            fasta_data.seq = std::regex_replace(full_entry, punct_pattern, "");
        }
        fasta_data.header = fasta_data.header.empty() ? std::to_string(fasta_data.rec_no) : fasta_data.header;
        AddRecordCount();
        return fasta_data;
    }
    else if constexpr (std::is_same_v<T, int>)
    {
        return std::stoi(full_entry);
    }
    else if constexpr (std::is_same_v<T, long>)
    {
        return std::stol(full_entry);
    }
    else if constexpr (std::is_same_v<T, char>)
    {
        return static_cast<T>(full_entry.front());
    }
    else
    {
        // Handle unsupported types here
        static_assert(std::is_same_v<T, T>, "Unsupported type");
    }
}

template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
{

    if constexpr (!std::is_same_v<T1, std::string> || !std::is_same_v<T1, FastaSequenceData>)
    {
        static_assert(std::is_same_v<T1, T1>, "Unsupported type, only std::string & FastaSequenceData are supported");
    }

    std::string delim_str(delim);

    std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> file_ptrs = MMapFile(filename, delim);

    char *start_of_file = std::get<1>(*file_ptrs).get();
    char *end_of_file = std::get<3>(*file_ptrs);

    char *p = start_of_file;

    arrow::RecordBatchVector ret_results;

#if defined(_OPENMP) || defined(WIN32)
    omp_lock_t pLock;
    omp_lock_t ret_resultsLock;
    omp_init_lock(&pLock);
    omp_init_lock(&ret_resultsLock);
#endif

#if defined(_OPENMP) || defined(WIN32)
#pragma omp parallel num_threads(num_threads) shared(end_of_file, start_of_file, delim) // entry_ptr_vec
#endif
    {
#if defined(_OPENMP) || defined(WIN32)
#pragma omp for schedule(dynamic) nowait // schedule(dynamic)
#endif
        for (int i = 0; i < num_threads; ++i)
        {

            // assert(!Progress::check_abort());
            // Rcpp::checkUserInterrupt();
            //  Get the thread-specific range to process
            size_t chunk_size = (end_of_file - start_of_file) / num_threads;
            char *thread_start = start_of_file + i * chunk_size;
            char *thread_end = (i == num_threads - 1) ? end_of_file : (thread_start + chunk_size);

            if (thread_start != start_of_file)
            {
                while (strncmp(thread_start, delim, strlen(delim)) != 0)
                {
                    --thread_start;
                }
            }
            if (thread_end != end_of_file)
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

            for (char *p = thread_start; p < thread_end; ++p)
            {
                // assert(!Progress::check_abort());

                if (strncmp(p, delim, strlen(delim)) == 0)
                {
                    entryStart = strstr(p, delim);
                    entryEnd = strstr(p + 1, delim);

                    if (entryEnd == nullptr)
                    {
                        // Handle the case where the delimiter is not found within the thread's range
                        entryEnd = end_of_file - 1;
                    }

                    // Process the entry from entryStart to entryEnd
                    std::string full_entry(entryStart, entryEnd - entryStart - 1);
                    if (!full_entry.empty())
                    {
                        AddRecordCount();
                        const T1 conv_entry = CastToType<T1>(full_entry);
                        if (Entry_callback != nullptr)
                        {
                            std::shared_ptr<arrow::RecordBatchVector> tmp_result = Entry_callback(std::make_shared<T1>(conv_entry));

                            if (return_values)
                            {
#if defined(_OPENMP) || defined(WIN32)
                                omp_set_lock(&ret_resultsLock);
#endif

                                ret_results.insert(ret_results.end(), tmp_result->begin(), tmp_result->end());
#if defined(_OPENMP) || defined(WIN32)
                                omp_unset_lock(&ret_resultsLock);
#endif
                            }
                            else
                            {
                                tmp_result->clear();
                            }
                        }
                    }
                    p = entryEnd - 1; // Move to the next position after the delimiter
                }
            }
        }
    }
#if defined(_OPENMP) || defined(WIN32)
#pragma omp barrier
#endif

#if defined(_OPENMP) || defined(WIN32)
    omp_destroy_lock(&pLock);
    omp_destroy_lock(&ret_resultsLock);
#endif

    return std::make_shared<arrow::RecordBatchVector>(ret_results);
}
