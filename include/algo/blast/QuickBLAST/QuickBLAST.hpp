#pragma once


#include <objmgr/bioseq_ci.hpp>

#include <iostream>
#include <memory>
#include <ncbi_pch.hpp>
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

#include <objects/seq/Seqdesc.hpp>    // CSeqdesc, CSeqdesc::SetTitle()
#include <objects/seq/Seq_descr.hpp>  // CSeq_descr
#include <objects/biblio/Title.hpp>   // CTitle (only required if you explicitly create a CTitle)
#include <objmgr/seqdesc_ci.hpp>

#include <algo/blast/api/blast_results.hpp>
#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/data_loaders/blastdb/bdbloader.hpp>
#include <algo/align/util/score_builder.hpp> // CScoreBuilder, AddScore, AddSplignScores
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/format/blastfmtutil.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <objtools/align_format/align_format_util.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <serial/iterator.hpp>
#include <algo/blast/api/bl2seq.hpp>
#include <objects/seq/Seq_literal.hpp>
#include <objtools/simple/simple_om.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <objmgr/seq_vector.hpp>
#include <objmgr/util/sequence.hpp>

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/DebugHelper.hpp>

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

 struct InterruptContext {
   std::atomic<bool> stop{false};
 };
 
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
        eFolder = 2,
        eBLASTDB = 3
    };


private:
    struct Impl;                 // Forward declaration of the implementation struct
    std::unique_ptr<Impl> pImpl; // Pointer to the implementation

public:
    QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences = false, bool save_hsp_sequences = false);
    ~QuickBLAST();
    
    std::shared_ptr<arrow::Schema> GetSchema();
    QuickBLAST::ESeqType GetSeqType();
    std::string GetProgram(void);
    void SetThreadCount(unsigned int num_threads);
    unsigned int GetThreadCount(void);
    int GetHitCount(void);
    void AddHitCount(int val = 1);
    unsigned int GetProcRecordCount();
    ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions();
    void ResetHitCount();

    unsigned int GetObjectID(void);
    void SetObjectID(unsigned int id);
    std::string GetQuickBLASTOptionString(void);
    // template <typename T1>
    std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim = "\n", const int &num_threads = 1, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback = {}, bool return_values = false);
    CRef<ncbi::blast::CBlastOptionsHandle> SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality );

    auto BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, bool verbose = true); //const bool show_progress = true
    std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values = false, int batch_size = 0, bool verbose = true); // const bool show_progress = true
    std::shared_ptr<arrow::RecordBatchVector> BLAST_f2db(const std::string &queryFile, const std::string &subjectDB, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values = false, unsigned int batch_size = 0, bool verbose = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_dbs(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values = false, unsigned int batch_size = 0, bool verbose = true); // const bool show_progress = true
    std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject, bool verbose = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_remote(const std::string &program, const std::string &database, const Rcpp::List &query_input, const QuickBLAST::EInputType input_type, std::string outFile, std::string outFormat,const bool return_values, const unsigned int max_poll_seconds, const unsigned int poll_interval_ms, bool verbose); //const unsigned int max_results
    unsigned int SizeOfDB(const std::string &dbName);
    SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb);
    SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector);
      
    template <typename T>
    std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values); 
    std::shared_ptr<arrow::RecordBatchVector> ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values);
};

struct QuickBLAST::Impl
{
    int num_threads = 4;
    std::string_view run_name;
    unsigned int obj_id;
    std::atomic<bool> quickblast_running{false};
    std::mutex cbl2seq_mutex;
    std::mutex builder_mutex;
    
    std::string program;
    ncbi::CRef<ncbi::blast::CBlastOptionsHandle> opts;
    std::string blast_options;
    ESeqType seq_type;
    EStrand strand;
    std::shared_ptr<ArrowWrapper> arrow_wrapper;
    std::shared_ptr<arrow::RecordBatch> empty_rb;
    int hit_count = 0;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_lock_t hit_countLock;
    omp_lock_t cleaner_threadsLock;
#endif
    bool save_sequences = false, save_hsp_sequences = false;
    std::vector<std::thread> cleaner_threads;
    std::shared_ptr<ArrowWrapper> GetArrowWrapper(void);
    std::shared_ptr<arrow::Schema> GetSchema(void);
    QuickBLAST::ESeqType GetSeqType();
    std::string GetProgram(void);
    void SetThreadCount(unsigned int num_threads);
    unsigned int GetThreadCount(void);
    int GetHitCount(void);
    void AddHitCount(int val = 1);
    unsigned int GetProcRecordCount();
    ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions(void);
    void ResetHitCount(void);
    unsigned int GetObjectID(void);
    void SetObjectID(unsigned int id);
    std::string GetQuickBLASTOptionString(void);
    CRef<ncbi::blast::CBlastOptionsHandle> SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality);
    
    Impl(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences, bool save_hsp_sequences);
    ~Impl();
    
    int GetFrame(int start, int length, ncbi::objects::ENa_strand strand);
    std::shared_ptr<arrow::RecordBatch> ExtractFASTA(const FastaSequenceData &fasta_data);
    std::string GetSSeqLocSequence(const SSeqLoc &seq_loc);
    SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name);
    SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb);
    SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector);
    std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, int batch_size = 0, bool verbose = true); // const bool show_progress
    std::shared_ptr<arrow::RecordBatchVector> BLAST_f2db(const std::string &queryFile, const std::string &subjectDB, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values = false, unsigned int batch_size = 0, bool verbose = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_dbs(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values = false, unsigned int batch_size = 0, bool verbose = true); // const bool show_progress = true
    std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject, bool verbose = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_remote(const std::string &program, const std::string &database, const Rcpp::List &query_input, const QuickBLAST::EInputType input_type, std::string outFile, std::string outFormat, const bool return_values,  const unsigned int max_poll_seconds, const unsigned int poll_interval_ms, bool verbose); 
    unsigned int SizeOfDB(const std::string &dbName);
    auto BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, bool verbose = true); 
    std::shared_ptr<arrow::RecordBatchVector> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const TSeqLocVector &sloc, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values);
    std::shared_ptr<arrow::RecordBatch> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values); 
    std::shared_ptr<arrow::RecordBatchVector> ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values);

private:
    std::pair<std::shared_ptr<ncbi::blast::SSeqLoc>, ncbi::CSeq_entry_Handle>  CreateSSeqLocFromType(const FastaSequenceData& fasta_data, ncbi::CRef<ncbi::CScope> parent_scope);
  bool IsProteinBioseq(const CBioseq_Handle &bh);
  bool GetFullSequenceString(CRef<CSeq_id> id, ncbi::CRef<ncbi::objects::CScope> scope, std::string &out_seq);
  bool GetHSPSequencesFromDenseg(const CDense_seg& dseg, ncbi::CRef<ncbi::objects::CScope> scope, std::string &q_hsp_ungapped, std::string &s_hsp_ungapped, std::string *q_aligned_with_gaps = nullptr, std::string *s_aligned_with_gaps = nullptr);
};

#endif // QUICKBLAST_HPP
