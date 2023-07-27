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
#include "arrow_wrapper.hpp"
#include <thread>
#include <iostream>
#include <iomanip>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rcpp.h>
#include <cassert>
#include <future>
#include <thread>
#ifdef _OPENMP
#include "omp.h"
#endif

USING_NCBI_SCOPE;
USING_SCOPE(blast);

class QuickBLAST
{
public:
    operator SEXP();
    enum ESeqType
    {
        eNucleotide = 0,
        eProtein = 1
    };
    enum EStrand
    {
        ePlus = 0,
        eMinus = 1
    };
    enum EInputType
    {
        eFile = 0,
        eSequenceString = 1
    };

    int num_threads = 4;
    std::string_view run_name;
    std::vector<std::string_view> queryFiles;
    std::vector<std::string_view> subjectFiles;
    std::string_view blastOptions = "-strand plus -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive\"";

private:
    std::string program;
    ncbi::blast::CBlastOptionsHandle *opts;
    ESeqType seq_type;
    EStrand strand;
    std::shared_ptr<ArrowWrapper> arrow_wrapper;
    int hit_count = 0;
    #ifdef _OPENMP
        omp_lock_t hit_countLock;
    #endif
    bool save_sequences = false;
    int blast_sequence_limit = 1000;
    bool db_scan_mode = false;
    std::promise<arrow::Status> ok_promise;

public:
    QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, Rcpp::List options, bool save_sequences = false);
    QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences = false);
    ~QuickBLAST();
    void PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream);
    void SetThreadCount(int num_threads)
    {
        this->num_threads = num_threads;
        #ifdef _OPENMP
            omp_set_num_threads(num_threads);
        #endif
        arrow_wrapper->SetThreadCount(num_threads);
    }
    int GetHitCount()
    {
        return hit_count;
    }
    void AddHitCount(int val = 1)
    {
        #ifdef _OPENMP
            omp_set_lock(&hit_countLock);
        #endif
        hit_count += val;
        #ifdef _OPENMP
            omp_unset_lock(&hit_countLock);
        #endif
    }
    ncbi::blast::CBlastOptionsHandle &GetQuickBLASTOptions()
    {
        return *opts;
    }
    void ResetHitCount() { hit_count = 0; }
    template <typename T1>
    std::shared_ptr<arrow::RecordBatchVector> StreamFile(const std::string_view &filename, const char *delim = "\n", const int &num_threads = 1, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback = {});
    template <typename OptionsType>
    ncbi::blast::CBlastOptionsHandle *SetQuickBLASTOptions(const std::string &program_name, OptionsType &options);

    Rcpp::List BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress = true);
    std::shared_ptr<arrow::RecordBatchVector> BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, int blast_sequence_limit, int num_threads, const bool show_progress = true, const bool return_values = false);
    std::shared_ptr<arrow::RecordBatch> BLAST_seqs(const std::string &query, const std::string &subject);
    SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb);
    SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector);

private:
    std::vector<std::pair<std::string, std::string>> BLASTOptionsFromString(const std::string &input);
    int CountCharacter(std::string filename, char character, int num_threads);
    void PrintProgressBar(int current, int total, int barWidth = 50);
    SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name);
    int GetFrame(int start, int length, ncbi::objects::ENa_strand strand);
    std::string GetSSeqLocSequence(const SSeqLoc &seq_loc);
    std::shared_ptr<arrow::RecordBatch> ExtractFASTA(const FastaSequenceData &fasta_data);
    template <typename T>
    std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope);
    template <typename T>
    SSeqLoc *CreateSSeqLocFromType(T _t, CRef<ncbi::CScope> parent_scope);
    SSeqLoc *CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope);
};
