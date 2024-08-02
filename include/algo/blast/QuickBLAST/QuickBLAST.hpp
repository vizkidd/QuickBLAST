#pragma once

// // #include "arrow_wrapper.hpp"
// #include <inttypes.h>
// #include <cstdint>
// #include <typeinfo>

// For portability - between mingw(gcc-12.3.0) and msys2(13.3.0) gcc versions
// #if typeid(int64_t).name() == typeid(long int).name()
// typedef int64_t long int
// #elif typeid(int64_t).name() == typeid(long long).name()
// typedef int64_t long long
// #endif
//// #define long int int64_t;

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#include "omp.h"
#endif

// TO make sure gcc doesnt interpret int64_t as long
// #ifdef MINGW32
// #include <inttypes.h>
// typedef int64_t long;
// // typedef int64_t LWT_INT64;
// // typedef uint64_t LWT_UINT64;
// // typedef int32_t LWT_INT32;
// // typedef uint32_t LWT_UINT32;
// #endif

// #if defined(linux) || defined(MINGW32)

// #include <RcppCommon.h>
// #include <Rcpp.h>
// #include "arrow_wrapper.hpp"
// #include "arrow_wrapper-functions.hpp"
// #elif defined(WIN32)
// #include <sdkddkver.h>
// #include <mman.h>
// #include <windows.h>
// #include <lmcons.h>
// // #include <pthread.h>
// #ifdef _OPENMP
// #include "omp.h"
// #endif
// #endif

#include <iostream>
#include <memory>
#include <ncbi_pch.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <objects/seqalign/Seq_align_set.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/format/blastfmtutil.hpp>
#include <objtools/align_format/align_format_util.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <serial/iterator.hpp>
#include <algo/blast/api/bl2seq.hpp>
#include <objects/seq/Seq_literal.hpp>
#include <objtools/simple/simple_om.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <objmgr/seq_vector.hpp>

#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/format/blastfmtutil.hpp>
#include <objtools/align_format/align_format_util.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <serial/iterator.hpp>
#include <algo/blast/api/bl2seq.hpp>
#include <objects/seq/Seq_literal.hpp>
#include <objtools/simple/simple_om.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <objmgr/seq_vector.hpp>

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

/*#if defined(MINGW32)
namespace arrow
{

    // Overload for compatability betwen int64_t in MinGW-Windows (long long) and int64_t in MSYS2-Linux(long)
    Status ArrayBuilder::AppendScalar(const Scalar &scalar, long value) // auto
    {
        // int64_t is long in MSYS2/Linux and long long in MinGW/Windows
        //  Call the original function with a cast
        return AppendScalar(scalar, static_cast<long long>(value));
    }

} // namespace arrow
#endif */

#ifndef QUICKBLAST_HPP
#define QUICKBLAST_HPP

class QuickBLAST
{
public:
    // operator SEXP();
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

    // int num_threads = 4;
    // std::string_view run_name;
    // int obj_id;

private:
    struct Impl;                 // Forward declaration of the implementation struct
    std::unique_ptr<Impl> pImpl; // Pointer to the implementation

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

public:
    // #ifdef linux
    // QBLIBRARY_API QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, Rcpp::List options, bool save_sequences = false);
    // #endif
    QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences = false);
    ~QuickBLAST();
    // void PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream);

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

    std::shared_ptr<arrow::Schema> GetSchema();
    void SetThreadCount(int num_threads);
    int GetHitCount();
    void AddHitCount(int val = 1);
    ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions();
    void ResetHitCount();

    unsigned int GetObjectID();
    void SetObjectID(unsigned int id);

    // template <typename T1>
    std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim = "\n", const int &num_threads = 1, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback = {}, bool return_values = false);
    // template <typename OptionsType>
    // ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options);
    // // ncbi::blast::CBlastOptionsHandle* SetQuickBLASTOptions(const std::string& program_name, const std::string& options);
    ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const std::string &options);

    // Rcpp::List BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress = true);
    auto BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false, int batch_size = 1024);
    std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject);
    SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb);
    SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector);

    template <typename T>
    std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope);
    // private:
    // std::vector<std::pair<std::string, std::string>> BLASTOptionsFromString(const std::string &input);
    // int CountCharacter(std::string filename, char character, int num_threads);
    // void PrintProgressBar(int current, int total, int barWidth = 50);
    // SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name);
    // int GetFrame(int start, int length, ncbi::objects::ENa_strand strand);
    // std::string GetSSeqLocSequence(const SSeqLoc &seq_loc);
    // std::shared_ptr<arrow::RecordBatch> ExtractFASTA(const FastaSequenceData &fasta_data);
    // template <typename T>
    // std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope);
    // template <typename T>
    // SSeqLoc *CreateSSeqLocFromType(T _t, CRef<ncbi::CScope> parent_scope);
    // SSeqLoc *CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope);
};

struct QuickBLAST::Impl
{
    int num_threads = 4;
    std::string_view run_name;
    unsigned int obj_id;

    std::string program;
    ncbi::CRef<ncbi::blast::CBlastOptionsHandle> opts;
    // Rcpp::List blast_options_list;
    // std::string blast_options_str;
    // SEXP blast_options;
    std::string blast_options;
    ESeqType seq_type;
    EStrand strand;
    std::shared_ptr<ArrowWrapper> arrow_wrapper;
    // ArrowWrapper arrow_wrapper;
    // Rcpp::XPtr<ArrowWrapper> arrow_wrapper;
    int hit_count = 0;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_lock_t hit_countLock;
#endif
    bool save_sequences = false;
    int blast_sequence_limit = 1000;
    // bool db_scan_mode = false;
    // std::promise<arrow::Status> ok_promise;

    std::shared_ptr<ArrowWrapper> GetArrowWrapper();

    std::shared_ptr<arrow::Schema> GetSchema();
    void SetThreadCount(int num_threads);
    int GetHitCount();
    void AddHitCount(int val = 1);
    ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions();
    void ResetHitCount();
    unsigned int GetObjectID();
    void SetObjectID(unsigned int id);

    ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, const std::string &options);
    // QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences)
    Impl(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences);

    // QuickBLAST::~QuickBLAST()
    ~Impl();

    // int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
    int GetFrame(int start, int length, ncbi::objects::ENa_strand strand);

    // std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data)
    std::shared_ptr<arrow::RecordBatch> ExtractFASTA(const FastaSequenceData &fasta_data);

    // std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc)
    std::string GetSSeqLocSequence(const SSeqLoc &seq_loc);

    // SEXP QuickBLAST::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
    SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name);

    // SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
    SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb);

    // SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
    SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector);

    // std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input)
    std::vector<std::pair<std::string, std::string>> BLASTOptionsFromString(const std::string &input);

    // template <typename T1>
    // std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
    // template <typename T1>
    std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values);

    // std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
    std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size);

    std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject);

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
    auto BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress);

    std::shared_ptr<arrow::RecordBatchVector> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const TSeqLocVector &sloc, const CScope &scope);
    std::shared_ptr<arrow::RecordBatch> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, const CScope &scope);

    SSeqLoc *CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope);
};

// namespace QuickBLASTUtils{
//         template <typename T>
//     std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope);
//     template <typename T>
//     SSeqLoc *CreateSSeqLocFromType(T _t, CRef<ncbi::CScope> parent_scope);
//     SSeqLoc *CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope);

// }

// template <typename T>
// SSeqLoc *QuickBLAST::CreateSSeqLocFromType(T _t, CRef<ncbi::CScope> parent_scope)
// {
//     return CreateSSeqLocFromType<FastaSequenceData>(arrow_wrapper->CastToType<FastaSequenceData>(_t), parent_scope);
// }

// template <typename T>
// std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> QuickBLAST::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope)
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

#endif // QUICKBLAST_HPP
