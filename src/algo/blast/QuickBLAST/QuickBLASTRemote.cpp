#include <algo/blast/QuickBLAST/commons.hpp>

// remote_blast_integration.cpp
#include <memory>
#include <string>
#include <thread>
#include <chrono>
#include <iostream>
#include <cstdio>

// typical toolkit includes (paths may vary by version)
#include <algo/blast/api/remote_blast.hpp>  
#include <algo/blast/api/blast_options_handle.hpp>
#include <objtools/blast/services/blast_services.hpp>
#include <objects/seqalign/Seq_align.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <corelib/ncbistd.hpp>

using namespace ncbi;
using namespace objects;
using namespace blast;
using namespace std::chrono_literals;
using namespace ncbi::objects;

#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_remote(
    const std::string &program,         // "blastn", "blastp", "tblastx", etc.
    const std::string &database,        // e.g. "nr" (or leave empty for defaults)
    const Rcpp::List &query_input,  // FASTA header+sequence or just sequence per your wrapper's expectations
    const QuickBLAST::EInputType input_type = EInputType::eSequenceString,
    std::string outFile = "",
    std::string outFormat = "parquet",
    const bool return_values = true,
    const unsigned int max_poll_seconds = 360,
    const unsigned int poll_interval_ms = 4000,
    bool verbose = true
){
  return pImpl->BLAST_remote(program, database, query_input, input_type, outFile, outFormat, return_values, max_poll_seconds, poll_interval_ms, verbose);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_remote(
    const std::string &program,         // "blastn", "blastp", "tblastx", etc.
    const std::string &database,        // e.g. "nr" (or leave empty for defaults)
    const Rcpp::List &query_input,  // FASTA header+sequence or just sequence per your wrapper's expectations
    const QuickBLAST::EInputType input_type = EInputType::eSequenceString,
    std::string outFile = "",
    std::string outFormat = "parquet",
    const bool return_values = true,
    const unsigned int max_poll_seconds = 360,
    const unsigned int poll_interval_ms = 4000,
    bool verbose = true
){

  if(outFile.empty() && return_values == false){
    Rcpp::Rcerr << "[BLAST_remote()] Error: Both outFile cannot be empty and return_values == FALSE." << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  if(max_poll_seconds <= 0){
    Rcpp::Rcerr << "[BLAST_remote()] Error: max_poll_seconds must be > 0s." << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  if(poll_interval_ms < 4000)
    Rcpp::Rcerr << "Warning: poll_interval < 4 seconds might not respect rate limits." << std::endl << std::flush;
  
  // if(outFile.empty()){
  //   outFile = std::tmpnam(nullptr); 
  // }
  
  arrow_wrapper->AddFASTAMetadata("BLAST locality", "CBlastOptions::eRemote");
  arrow_wrapper->AddFASTAMetadata("Input source", "sequence");
  arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);

  CBlastServices blast_service;
  if(!blast_service.IsValidBlastDb(database, (seq_type == ESeqType::eProtein))){
    Rcpp::Rcerr << "BLAST_remote: Not a valid NCBI database." << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  CRef<objects::CBioseq_set> bss(new objects::CBioseq_set());
  std::vector<CSeq_entry_Handle> subject_ent_vec;
  bss->SetClass(objects::CBioseq_set::eClass_nuc_prot); //EClass - https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/include/objects/seqset/Bioseq_set_.hpp
  
  CRef<CObjectManager> objMgr(CObjectManager::GetInstance());
  if (!objMgr) {
    Rcpp::Rcerr << "BLAST_remote: CObjectManager::GetInstance() returned NULL." << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  CRef<ncbi::CScope> scope(new ncbi::CScope(*objMgr));
  scope->AddDefaults();
  
  for(const std::string &query : query_input){
    Rcpp::checkUserInterrupt();
    auto q_type = this->arrow_wrapper->CastToType(query);
    if (q_type.header.empty() || q_type.seq.empty()) {
      Rcpp::Rcerr << "BLAST_remote: query header/sequence is empty." << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    int rec_no = q_type.rec_no;
    std::string fastaID(q_type.header.data());
    std::string fastaSequence(q_type.seq.data());
    
    const TSeqPos seqlen = fastaSequence.length();
    
    if(seqlen >= std::numeric_limits<TSeqPos>::max()){
      Rcpp::stop("[BLAST_remote()] seqlen >= std::numeric_limits<TSeqPos>::max().");
    }
  
    CRef<CSeq_id> id(new CSeq_id(fastaID, CSeq_id::fParse_RawText | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal));
    
    CRef<CSeq_loc>
      cseq_loc_obj(new CSeq_loc());
    cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
    cseq_loc_obj->SetId(*id);
    if (seq_type == ESeqType::eProtein)
    {
      cseq_loc_obj->SetStrand(eNa_strand_unknown);
    }
    else
    {
      switch (strand)
      {
      case EStrand::ePlus:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_plus);
        break;
      case EStrand::eMinus:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_minus);
        break;
      case EStrand::eUnknown:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_unknown);
        break;
      case EStrand::eBoth:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both);
        break;
      case EStrand::eBoth_rev:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both_rev);
        break;
      case EStrand::eOther:
        cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_other);
        break;
      }
    }
    
    CRef<CSeq_data> seq_data(new CSeq_data());
    seq_data->Select(seq_type == ESeqType::eProtein ? CSeq_data_Base::E_Choice::e_Iupacaa : CSeq_data_Base::E_Choice::e_Iupacna);
    switch (seq_type)
    {
    case ESeqType::eProtein:
    {
      seq_data->SetIupacaa(CIUPACaa(fastaSequence));
      break;
    }
    case ESeqType::eNucleotide:
    {
      seq_data->SetIupacna(CIUPACna(fastaSequence));
      break;
    }
    }
    
    CRef<CSeq_inst> seq_inst(new CSeq_inst());
    seq_inst->SetSeq_data(*seq_data);
    seq_inst->SetLength(fastaSequence.length());
    seq_inst->SetRepr(CSeq_inst::eRepr_raw);
    seq_inst->SetTopology(CSeq_inst::eTopology_linear);
    
    if (seq_type == ESeqType::eProtein)
    {
      seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_aa);
    }
    else
    {
      seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_na);
    }
    
    seq_inst->SetLength(seqlen);
    
    CRef<CBioseq> bioseq(new CBioseq());
    bioseq->SetId().push_back(id);      // ID that matches cseq_loc above
    bioseq->SetInst(*seq_inst);
    
    // optionally add description/title:
    // CRef<GBSeq> descr(new GBSeq()); ... bioseq->SetDescr(...);
    
    // -- put into an entry and add to scope (so the scope can resolve lcl|rec_no)
    CRef<CSeq_entry> ret_entry(new CSeq_entry());
    ret_entry->SetSeq(*bioseq);
    
    CSeq_entry_Handle tse_handle = scope->AddTopLevelSeqEntry(*ret_entry);
    subject_ent_vec.emplace_back(tse_handle);
    //DEBUG
    // for (auto &id0 : ret_entry->GetSeq().GetId()) {
    //   Rcpp::Rcout << "Added seq id: " << id0->AsFastaString() << std::endl;
    // }
    // try {
    //   CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(*ret_entry->GetSeq().GetId().front());
    //   CBioseq_Handle bh = scope->GetBioseqHandle(idh);
    //   Rcpp::Rcout << "Scope resolves id: " << idh.AsString() << std::endl;
    // } catch (const CException &e) {
    //   Rcpp::Rcerr << "Scope resolution failed after add: " << e.GetMsg() << std::endl;
    // }
    
    bss->SetSeq_set().push_back(ret_entry);
  }
  
  // // OPTIONAL: set number of results/hits etc
  // if(max_results > 0)
  //   remote.SetMaxTargetSequences(max_results); // API-dependent
  
  // Build CRemoteBlast
  CRemoteBlast remote(SetQuickBLASTOptions(program, GetQuickBLASTOptionString(), CBlastOptions::eRemote, verbose));
  remote.SetQueries(bss);
  remote.SetDatabase(database);
  // Submit synchronously, wait, then get results
  // remote.SubmitSync();
  try {
    if(verbose)
      Rcpp::Rcout << "Max wait time: " << max_poll_seconds << " seconds" << std::endl << std::flush;
    remote.Submit();                 // may return void or RID;
  }
  catch (const CBlastException &e) {
    // treat as network / submission failure
    throw std::runtime_error(std::string("Remote submit failed: ") + e.what());
  }
  
  Progress progress_bar(max_poll_seconds, verbose);
  
  // Use ms accumulator to avoid integer-division artifacts
  unsigned long long waited_ms = 0ULL;
  unsigned int reported_seconds = 0; // how many seconds we've reported to progress
  
  auto status = remote.CheckStatus(); // initial status
  
  while (status == CRemoteBlast::ESearchStatus::eStatus_Pending || status != CRemoteBlast::ESearchStatus::eStatus_Unknown) {
    Rcpp::checkUserInterrupt();
    
    // If remote changed between checks, break early
    if (status == CRemoteBlast::ESearchStatus::eStatus_Done) break;
    
    if (status == CRemoteBlast::ESearchStatus::eStatus_Failed) {
      std::vector<std::string> remoteErrors = remote.GetErrorVector();
      for (const std::string &error : remoteErrors) {
        Rcpp::Rcerr << "[BLAST_remote()] Error: Remote BLAST reported an error while processing the job : " << error << std::endl << std::flush;
      }
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    // Sleep for the requested interval
    std::this_thread::sleep_for(std::chrono::milliseconds(poll_interval_ms));
    waited_ms += static_cast<unsigned long long>(poll_interval_ms);
    
    // Compute whole seconds elapsed
    unsigned int seconds_elapsed = static_cast<unsigned int>(waited_ms / 1000ULL);
    
    // Increment progress for each whole second that passed since last report
    while (reported_seconds < seconds_elapsed && reported_seconds < static_cast<unsigned int>(max_poll_seconds)) {
      progress_bar.increment();
      ++reported_seconds;
    }
    
    // Time-out check (use seconds_elapsed to compare to max_poll_seconds)
    if (seconds_elapsed > static_cast<unsigned int>(max_poll_seconds)) {
      // optionally: remote.Cancel() if available
      throw std::runtime_error("Remote BLAST timed out");
    }
    
    // Re-check remote status
    status = remote.CheckStatus();
  }
  
  // If we exit before the bar reached max, you can optionally finish it:
  while (reported_seconds < static_cast<unsigned int>(max_poll_seconds)) {
    progress_bar.increment(); ++reported_seconds;
  }
  
  CRef<CSearchResultSet> results;
  auto rid = remote.GetRID();
  
  if(verbose)
    Rcpp::Rcout << "(Success) Remote BLAST run ID: " << rid << std::endl << std::flush; //DEBUG
  
  try {
    results = remote.GetResultSet();
  } catch (const ncbi::CException& e) {
    // Handle cases where the server might return an error
    Rcpp::Rcerr << std::string("[BLAST_remote()]: 1. NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  if (results.IsNull() || results->GetNumResults() == 0){
    Rcpp::Rcerr << "[BLAST_remote()] Remote BLAST search completed: No hits found." << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }
  // 5) Fetch result as Seq-align(s)
  // Many APIs provide a method like: remote.GetSeqAlignSet() or remote.GetSeqAligns()
  // We want to produce a TSeqAlignVector (same type your ExtractHits expects).
  TSeqAlignVector alignments;
  
  try {
    // The toolkit often returns a CRef<CSeq_align_set> or vector<CRef<CSeq_align_set>>
    CRef<CSeq_align_set> align_set = remote.GetAlignments();
    if (align_set.NotEmpty() && !align_set->IsEmpty()) { 
      alignments.emplace_back(align_set);
    } else {
      Rcpp::Rcerr << "[BLAST_remote()] No alignments found for this RID." << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    // AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);
    
    CBlastServices::TSeqIdVector seqids;
    seqids.reserve(256);
    std::unordered_set<std::string> seen; // de-dup by string
    
    for (const auto &align_set_ref : alignments) {
      if (!align_set_ref) continue;
      for (const auto &align_ref : align_set_ref->Get()) {
        if (!align_ref) continue;
        // find seg type and extract the subject row id(s).
        if (align_ref->GetSegs().IsDenseg()) {
          const CDense_seg &dseg = align_ref->GetSegs().GetDenseg();
          // usually row 0 = query, row 1 = subject for pairwise BLAST
          if (dseg.CanGetIds() && dseg.GetIds().size() >= 2) {
            // pick the subject id (row 1)
            const CRef<CSeq_id> &sid = dseg.GetIds()[1];
            if (sid) {
              // use a string key—GetSequence will accept CRef<CSeq_id> but we
              // want to de-duplicate identical ones.
              std::string key = sid->GetSeqIdString(true); // canonical string
              if (seen.emplace(key).second) {
                // make a *copy* of the seq-id for the request
                // CRef<CSeq_id> sid_copy(new CSeq_id(*sid));
                seqids.push_back(sid); //align_ref->GetSeq_id(1) //sid_copy
              }
            }
          }
        }
        
      }
    }
    
    if (seqids.empty()) {
      Rcpp::Rcerr << "RemoteBlast() - Query aligned to no subject IDs. (Subject IDs.size() == 0)" << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    // call the remote get-sequences service
    CBlastServices::TBioseqVector bioseqs;
    std::string errors, warnings;
    try {
      Rcpp::checkUserInterrupt();
      CBlastServices::GetSequences(seqids, database, seq_type == ESeqType::eProtein ? 'p' : 'n' , bioseqs, errors, warnings, /*verbose=*/false);
    }
    catch (const CException &e) {
      // toolkit exception — handle/log
      Rcpp::Rcerr << std::string("[RemoteBlast()] - GetSequences failed: ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    if (!errors.empty()) {
      Rcpp::Rcerr << std::string("RemoteBlast() - GetSequences errors: ") + errors << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    if (!warnings.empty()) {
      Rcpp::Rcerr << std::string("RemoteBlast() - GetSequences warnings: ") + warnings << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    // add returned Bioseqs to the scope so CSeqVector/CBioseq_Handle lookups work
    for (auto &bioseq_ref : bioseqs) {
      if (!bioseq_ref || bioseq_ref.IsNull()) continue;
      Rcpp::checkUserInterrupt();
      CRef<CSeq_entry> entry(new CSeq_entry());
      entry->SetSeq(*bioseq_ref);                 // copy CBioseq into Seq-entry
      scope->AddTopLevelSeqEntry(*entry);         // now scope can resolve those seq-ids
      // // Debug: print the IDs added
      // for (const auto &id : bioseq_ref->GetId()) {
      //   LOG_POST("Loaded subject seq id: " << id->GetSeqIdString(true));
      // }
    }
    
    return ExtractHitsRemote(alignments, subject_ent_vec, scope, return_values);

  }catch(const ncbi::CException& e) {
    Rcpp::Rcerr << std::string("[Blast_remote()]: 2. NCBI CException : Failed fetching remote results: ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[Blast_remote()]: C++ Runtime Error : Failed fetching remote results: ") + e.what() << std::endl << std::flush; 
    return std::make_shared<arrow::RecordBatchVector>();
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[Blast_remote()]: Rcpp Runtime Error : Failed fetching remote results: ") + e.what() << std::endl << std::flush; 
    return std::make_shared<arrow::RecordBatchVector>();
  }catch (const std::exception &e) {
    throw std::runtime_error(std::string("[Blast_remote()]: C++ Exception : Failed fetching remote results: ") + e.what());
  }catch (...) {
    throw std::runtime_error(std::string("[Blast_remote()]: Unknown Error."));
  }
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values){
  return pImpl->ExtractHitsRemote(alignments, sseq_entry_vec, scope, return_values);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values) 
{
  try{
    std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
    
    // quick sanity
    if (alignments.empty()) {
      Rcpp::Rcerr << "TSeqAlignVector - No alignments could be computed." << std::endl << std::flush;
      ret_val->emplace_back(empty_rb);
      return ret_val;
    }
    
    // Arrow Builders
    arrow::Int64Builder hsp_offset_builder, length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    arrow::StringBuilder frames_builder, strand_builder, qseqid_builder, sseqid_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder, qhsp_builder, shsp_builder;
    arrow::Int64Builder qlen_builder, slen_builder, num_alignments_builder;
    
    CScoreBuilder scorer;
    // if (effective_search_space > 0.0) scorer.SetEffectiveSearchSpace(effective_search_space);
    
    // Compute batch scores (AddScore has an overload for list)
    std::vector<CSeq_align::EScoreType> score_types = {
      CSeq_align::EScoreType::eScore_AlignLength,
      CSeq_align::EScoreType::eScore_BitScore,
      CSeq_align::EScoreType::eScore_Blast,
      CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped,
      CSeq_align::EScoreType::eScore_PercentIdentity,
      CSeq_align::EScoreType::eScore_GapCount,
      CSeq_align::EScoreType::eScore_EValue,
      CSeq_align::EScoreType::eScore_IdentityCount,
      CSeq_align::EScoreType::eScore_MismatchCount,
      CSeq_align::EScoreType::eScore_PercentCoverage,
      CSeq_align::EScoreType::eScore_Score,
      CSeq_align::EScoreType::eScore_PositiveCount,
      CSeq_align::EScoreType::eScore_Splices,
      CSeq_align::EScoreType::eScore_SumEValue,
      CSeq_align::EScoreType::eScore_ProductCoverage,
      CSeq_align::EScoreType::eScore_OverallIdentity,
      CSeq_align::EScoreType::eScore_NegativeCount,
      CSeq_align::EScoreType::eScore_Matches,
      CSeq_align::EScoreType::eScore_HighQualityPercentCoverage,
      CSeq_align::EScoreType::eScore_ExonIdentity,
      CSeq_align::EScoreType::eScore_ConsensusSplices,
      CSeq_align::EScoreType::eScore_CompAdjMethod
    };
    
    int64_t estimated_rows = 0;
    for (const auto &align_set_ref : alignments) {
      if (align_set_ref && align_set_ref->IsSet() && align_set_ref->CanGet()) {
        estimated_rows += align_set_ref->Get().size();
      }
    }

    if (estimated_rows == 0) {
      ret_val->emplace_back(empty_rb);
      return ret_val;
    }
    
    // reserve arrow containers for estimated row count
    static_cast<void>(hsp_offset_builder.Reserve(estimated_rows));
    static_cast<void>(length_builder.Reserve(estimated_rows));
    static_cast<void>(mismatch_builder.Reserve(estimated_rows));
    static_cast<void>(gapopen_builder.Reserve(estimated_rows));
    static_cast<void>(qstart_builder.Reserve(estimated_rows));
    static_cast<void>(qend_builder.Reserve(estimated_rows));
    static_cast<void>(sstart_builder.Reserve(estimated_rows));
    static_cast<void>(send_builder.Reserve(estimated_rows));
    static_cast<void>(gaps_builder.Reserve(estimated_rows));
    static_cast<void>(nident_builder.Reserve(estimated_rows));
    static_cast<void>(positive_builder.Reserve(estimated_rows));
    static_cast<void>(n_splices_builder.Reserve(estimated_rows));
    static_cast<void>(hsp_cnt_builder.Reserve(estimated_rows));
    static_cast<void>(negative_count_builder.Reserve(estimated_rows));
    static_cast<void>(pident_builder.Reserve(estimated_rows));
    static_cast<void>(pident_gap_builder.Reserve(estimated_rows));
    static_cast<void>(evalue_builder.Reserve(estimated_rows));
    static_cast<void>(bitscore_builder.Reserve(estimated_rows));
    static_cast<void>(score_builder.Reserve(estimated_rows));
    static_cast<void>(qcovhsp_builder.Reserve(estimated_rows));
    static_cast<void>(blast_score_builder.Reserve(estimated_rows));
    static_cast<void>(aln_len01_builder.Reserve(estimated_rows));
    static_cast<void>(sum_evalue_builder.Reserve(estimated_rows));
    static_cast<void>(product_coverage_builder.Reserve(estimated_rows));
    static_cast<void>(overall_identity_builder.Reserve(estimated_rows));
    static_cast<void>(matches_builder.Reserve(estimated_rows));
    static_cast<void>(high_quality_percent_coverage_builder.Reserve(estimated_rows));
    static_cast<void>(exon_identity_builder.Reserve(estimated_rows));
    static_cast<void>(consensus_splices_builder.Reserve(estimated_rows));
    static_cast<void>(comp_adj_method_builder.Reserve(estimated_rows));
    static_cast<void>(frames_builder.Reserve(estimated_rows));
    static_cast<void>(strand_builder.Reserve(estimated_rows));
    static_cast<void>(qseqid_builder.Reserve(estimated_rows));
    static_cast<void>(sseqid_builder.Reserve(estimated_rows));
    static_cast<void>(qseq_builder.Reserve(estimated_rows));
    static_cast<void>(sseq_builder.Reserve(estimated_rows));
    static_cast<void>(qhsp_builder.Reserve(estimated_rows));
    static_cast<void>(shsp_builder.Reserve(estimated_rows));
    static_cast<void>(qlen_builder.Reserve(estimated_rows));
    static_cast<void>(slen_builder.Reserve(estimated_rows));
    static_cast<void>(num_alignments_builder.Reserve(estimated_rows));
    
    std::string q_full, s_full, q_hsp, s_hsp, q_aligned, s_aligned;
    std::string qseq_id, sseq_id, strand_str, frames;
    
    // Pre-reserve common string capacities to avoid growth during appending
    q_full.reserve(4096); s_full.reserve(4096);
    q_hsp.reserve(1024); s_hsp.reserve(1024);
    
    int num_rows = 0;
    
    for (const auto &align_set_ref : alignments) {
      if (!align_set_ref || !align_set_ref->IsSet()) continue;
      RcppThread::checkUserInterrupt();
      
      auto &seq_align_list = align_set_ref->Get(); 
      int64_t parent_list_size = seq_align_list.size();
      
      for (const auto &seq_align : seq_align_list) {
        if (!seq_align || seq_align.IsNull() || !seq_align.NotEmpty()) continue;
        seq_align->Validate(true);
        RcppThread::checkUserInterrupt();
        
        // Clear hoisted strings for buffer reuse
        q_full.clear(); s_full.clear();
        q_hsp.clear(); s_hsp.clear(); 
        q_aligned.clear(); s_aligned.clear();
        qseq_id.clear(); sseq_id.clear();
        
        try {
          qseq_id = seq_align->GetSeq_id(0).GetSeqIdString(true);
        } catch (...) { qseq_id = "(unknown_query)"; }
        try {
          sseq_id = seq_align->GetSeq_id(1).GetSeqIdString(true);
        } catch (...) { sseq_id = "(unknown_subject)"; }
        
        // Branchless/Direct Strand Construction ---
        ncbi::objects::ENa_strand q_strand = seq_align->GetSeqStrand(0);
        ncbi::objects::ENa_strand s_strand = seq_align->GetSeqStrand(1);
        
        char q_strand_char = (q_strand == ncbi::objects::eNa_strand_plus) ? '+' : 
          (q_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
        char s_strand_char = (s_strand == ncbi::objects::eNa_strand_plus) ? '+' : 
          (s_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
        
        strand_str = std::string(1, q_strand_char) + "/" + s_strand_char;
        
        if (seq_align->GetSegs().IsDenseg()) {
          const auto& dseg = seq_align->GetSegs().GetDenseg();
          if (save_sequences && dseg.CanGetIds()) {
            if (dseg.GetIds().size() > 0) 
              GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[0]), scope, q_full);
            if (dseg.GetIds().size() > 1) 
              GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[1]), scope, s_full);
          }
          if (save_hsp_sequences) {
            GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned);
          }
        }
        
        // Variable grouping (Primitive types are cheap to initialize locally)
        double score = 0, n_splices = 0, num_ident = 0, aln_len = 0, gaps = 0, mismatches = 0, positive = 0, negative_count = 0;
        double bits = 0, evalue = 0, blast_score = 0, pident = 0, aln_len01 = 0, pident_gap = 0, qcovhsp = 0;
        double sum_evalue = 0, product_coverage = 0, overall_identity = 0, matches = 0, high_quality_percent_coverage = 0, exon_identity = 0, consensus_splices = 0, comp_adj_method = 0;
        
        // NCBI Score Extractions 
        if(!seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len)) {
          aln_len = seq_align->GetAlignLength(true);
        }
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); 
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
        bool hasid = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
        
        if (!seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident) && hasid) {
          pident = 100.0 * num_ident / seq_align->GetAlignLength(false); 
        }
        if (!seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap) && hasid) {
          pident_gap = 100.0 * num_ident / seq_align->GetAlignLength(true); 
        }
        if(!seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps)) {
          // gaps = seq_align->GetTotalGapCount(-1);
          // Proactively check the segment type to prevent NCBI_THROW from logging errors
          const auto& segs = seq_align->GetSegs();
          
          if (segs.IsDenseg() || segs.IsPacked()) {
            // Standard alignments (blastn, blastp, blastx) support this natively
            gaps = seq_align->GetTotalGapCount(-1);
          } else {
            // tblastx uses Std-seg or Dendiag which represent ungapped blocks.
            // Since there are no gaps in these specific alignment formats, it is exactly 0.
            gaps = 0;
          }
        }
        
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); 
        if(!seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches)) {
          mismatches = seq_align->GetAlignLength(true) - num_ident - gaps;
        }
        
        qcovhsp = q_full.length() > 0 ? (static_cast<double>(seq_align->GetAlignLength(false)) / static_cast<double>(q_full.length())) : 0.0;
        
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score); 
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices); 
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices);
        seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method);
        
        aln_len01 = seq_align->AlignLengthRatio();
        int qstart = seq_align->GetSeqStart(0) + 1; 
        int qend   = seq_align->GetSeqStop(0) + 1;
        int sstart = seq_align->GetSeqStart(1) + 1;
        int send   = seq_align->GetSeqStop(1) + 1;
        
        frames = std::to_string(GetFrame(qstart, aln_len, q_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, s_strand));
        
        // Append to builders
        static_cast<void>(qhsp_builder.Append(q_hsp));
        static_cast<void>(shsp_builder.Append(s_hsp));
        static_cast<void>(frames_builder.Append(frames));
        static_cast<void>(strand_builder.Append(strand_str));
        static_cast<void>(qseq_builder.Append(save_sequences ? q_full : ""));
        static_cast<void>(sseq_builder.Append(save_sequences ? s_full : ""));
        static_cast<void>(qseqid_builder.Append(qseq_id));
        static_cast<void>(sseqid_builder.Append(sseq_id));
        
        // Use UnsafeAppend because we already Reserved capacity!
        static_cast<void>(qstart_builder.UnsafeAppend(qstart));
        static_cast<void>(qend_builder.UnsafeAppend(qend));
        static_cast<void>(sstart_builder.UnsafeAppend(sstart));
        static_cast<void>(send_builder.UnsafeAppend(send));
        static_cast<void>(pident_builder.UnsafeAppend(pident));
        static_cast<void>(evalue_builder.UnsafeAppend(evalue));
        static_cast<void>(length_builder.UnsafeAppend(aln_len));
        static_cast<void>(aln_len01_builder.UnsafeAppend(aln_len01));
        static_cast<void>(bitscore_builder.UnsafeAppend(bits));
        static_cast<void>(score_builder.UnsafeAppend(score));
        static_cast<void>(qcovhsp_builder.UnsafeAppend(qcovhsp));
        static_cast<void>(blast_score_builder.UnsafeAppend(blast_score));
        static_cast<void>(pident_gap_builder.UnsafeAppend(pident_gap));
        static_cast<void>(gaps_builder.UnsafeAppend(gaps));
        static_cast<void>(nident_builder.UnsafeAppend(num_ident));
        static_cast<void>(mismatch_builder.UnsafeAppend(mismatches));
        static_cast<void>(positive_builder.UnsafeAppend(positive));
        static_cast<void>(n_splices_builder.UnsafeAppend(n_splices));
        static_cast<void>(hsp_cnt_builder.UnsafeAppend(num_rows + 1));
        static_cast<void>(sum_evalue_builder.UnsafeAppend(sum_evalue));
        static_cast<void>(product_coverage_builder.UnsafeAppend(product_coverage));
        static_cast<void>(overall_identity_builder.UnsafeAppend(overall_identity));
        static_cast<void>(negative_count_builder.UnsafeAppend(negative_count));
        static_cast<void>(matches_builder.UnsafeAppend(matches));
        static_cast<void>(high_quality_percent_coverage_builder.UnsafeAppend(high_quality_percent_coverage));
        static_cast<void>(exon_identity_builder.UnsafeAppend(exon_identity));
        static_cast<void>(consensus_splices_builder.UnsafeAppend(consensus_splices));
        static_cast<void>(comp_adj_method_builder.UnsafeAppend(comp_adj_method));
        
        TSeqPos qlen = scope->GetBioseqHandle(seq_align->GetSeq_id(0)).GetBioseqLength();
        TSeqPos slen = scope->GetBioseqHandle(seq_align->GetSeq_id(1)).GetBioseqLength();
        static_cast<void>(qlen_builder.UnsafeAppend(qlen));
        static_cast<void>(slen_builder.UnsafeAppend(slen));
        static_cast<void>(num_alignments_builder.UnsafeAppend(parent_list_size));
        static_cast<void>(hsp_offset_builder.UnsafeAppend(1));
        
        num_rows++;
        
      } 
    }
    
    // pour the builder contents into arrays
    std::shared_ptr<arrow::Array> qhsp_array;
    static_cast<void>(qhsp_builder.Finish(&qhsp_array));
    std::shared_ptr<arrow::Array> shsp_array;
    static_cast<void>(shsp_builder.Finish(&shsp_array));
    std::shared_ptr<arrow::Array> frames_array;
    static_cast<void>(frames_builder.Finish(&frames_array));
    std::shared_ptr<arrow::Array> pident_array;
    static_cast<void>(pident_builder.Finish(&pident_array));
    std::shared_ptr<arrow::Array> pident_gap_array;
    static_cast<void>(pident_gap_builder.Finish(&pident_gap_array));
    std::shared_ptr<arrow::Array> evalue_array;
    static_cast<void>(evalue_builder.Finish(&evalue_array));
    std::shared_ptr<arrow::Array> length_array;
    static_cast<void>(length_builder.Finish(&length_array));
    std::shared_ptr<arrow::Array> qstart_array;
    static_cast<void>(qstart_builder.Finish(&qstart_array));
    std::shared_ptr<arrow::Array> qend_array;
    static_cast<void>(qend_builder.Finish(&qend_array));
    std::shared_ptr<arrow::Array> sstart_array;
    static_cast<void>(sstart_builder.Finish(&sstart_array));
    std::shared_ptr<arrow::Array> send_array;
    static_cast<void>(send_builder.Finish(&send_array));
    std::shared_ptr<arrow::Array> aln_len01_array;
    static_cast<void>(aln_len01_builder.Finish(&aln_len01_array));
    std::shared_ptr<arrow::Array> bitscore_array;
    static_cast<void>(bitscore_builder.Finish(&bitscore_array));
    std::shared_ptr<arrow::Array> score_array;
    static_cast<void>(score_builder.Finish(&score_array));
    std::shared_ptr<arrow::Array> qcovhsp_array;
    static_cast<void>(qcovhsp_builder.Finish(&qcovhsp_array));
    std::shared_ptr<arrow::Array> blast_score_array;
    static_cast<void>(blast_score_builder.Finish(&blast_score_array));
    std::shared_ptr<arrow::Array> gaps_array;
    static_cast<void>(gaps_builder.Finish(&gaps_array));
    std::shared_ptr<arrow::Array> nident_array;
    static_cast<void>(nident_builder.Finish(&nident_array));
    std::shared_ptr<arrow::Array> mismatch_array;
    static_cast<void>(mismatch_builder.Finish(&mismatch_array));
    std::shared_ptr<arrow::Array> positive_array;
    static_cast<void>(positive_builder.Finish(&positive_array));
    std::shared_ptr<arrow::Array> n_splices_array;
    static_cast<void>(n_splices_builder.Finish(&n_splices_array));
    std::shared_ptr<arrow::Array> hsp_cnt_array;
    static_cast<void>(hsp_cnt_builder.Finish(&hsp_cnt_array));
    std::shared_ptr<arrow::Array> sum_evalue_array;
    static_cast<void>(sum_evalue_builder.Finish(&sum_evalue_array));
    std::shared_ptr<arrow::Array> product_coverage_array;
    static_cast<void>(product_coverage_builder.Finish(&product_coverage_array));
    std::shared_ptr<arrow::Array> overall_identity_array;
    static_cast<void>(overall_identity_builder.Finish(&overall_identity_array));
    std::shared_ptr<arrow::Array> negative_count_array;
    static_cast<void>(negative_count_builder.Finish(&negative_count_array));
    std::shared_ptr<arrow::Array> matches_array;
    static_cast<void>(matches_builder.Finish(&matches_array));
    std::shared_ptr<arrow::Array> high_quality_percent_coverage_array;
    static_cast<void>(high_quality_percent_coverage_builder.Finish(&high_quality_percent_coverage_array));
    std::shared_ptr<arrow::Array> exon_identity_array;
    static_cast<void>(exon_identity_builder.Finish(&exon_identity_array));
    std::shared_ptr<arrow::Array> consensus_splices_array;
    static_cast<void>(consensus_splices_builder.Finish(&consensus_splices_array));
    std::shared_ptr<arrow::Array> comp_adj_method_array;
    static_cast<void>(comp_adj_method_builder.Finish(&comp_adj_method_array));
    
    // package the arrays
    arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({
      qhsp_array,
      shsp_array,
      pident_array,
      pident_gap_array,
      frames_array,
      evalue_array,
      length_array,
      aln_len01_array,
      qstart_array,
      qend_array,
      sstart_array,
      send_array,
      bitscore_array,
      score_array,
      qcovhsp_array,
      blast_score_array,
      gaps_array,
      nident_array,
      mismatch_array,
      positive_array,
      n_splices_array,
      hsp_cnt_array,
      sum_evalue_array,
      product_coverage_array,
      overall_identity_array,
      negative_count_array,
      matches_array,
      high_quality_percent_coverage_array,
      exon_identity_array,
      consensus_splices_array,
      comp_adj_method_array},
      {"qhsp", "shsp", "pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});
    
    if (!aln_struct_array.ok()) {
      throw std::runtime_error("[ExtractHitsRemote()] 1. Failed to build StructArray: " + aln_struct_array.status().ToString());
    }
    
    std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
    
    // pour the seq-meta builder contents into arrays
    std::shared_ptr<arrow::Array> qseqid_array;
    static_cast<void>(qseqid_builder.Finish(&qseqid_array));
    std::shared_ptr<arrow::Array> sseqid_array;
    static_cast<void>(sseqid_builder.Finish(&sseqid_array));
    std::shared_ptr<arrow::Array> qseq_array;
    static_cast<void>(qseq_builder.Finish(&qseq_array));
    std::shared_ptr<arrow::Array> sseq_array;
    static_cast<void>(sseq_builder.Finish(&sseq_array));
    std::shared_ptr<arrow::Array> qlen_array;
    static_cast<void>(qlen_builder.Finish(&qlen_array));
    std::shared_ptr<arrow::Array> slen_array;
    static_cast<void>(slen_builder.Finish(&slen_array));
    std::shared_ptr<arrow::Array> strand_array;
    static_cast<void>(strand_builder.Finish(&strand_array));
    std::shared_ptr<arrow::Array> num_alignment_array;
    static_cast<void>(num_alignments_builder.Finish(&num_alignment_array));
    
    // Create the seq_info struct array and populate with the arrays
    std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
                                                                                                 seqids_struct_array,
                                                                                                 seqs_struct_array,
                                                                                                 strand_array,                                                                                               lengths_struct_array},
                                                                                                 {"num_alignments", "seqids", "seqs", "strands", "lengths"});
    
    if(!seq_info_array.ok()){
      std::runtime_error("[ExtractHitsRemote()] 2. Failed to build StructArray: " + seq_info_array.status().ToString());
    }
    
    std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
    
    std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                                num_rows,
                                                                                {seq_info_array_, aln_struct_array_});
    
    for(auto s_ent: sseq_entry_vec){
      scope->RemoveTopLevelSeqEntry(s_ent);
      // s_ent.Reset();
    }
    sseq_entry_vec.clear();
    sseq_entry_vec.shrink_to_fit();
    
    if(alignment_rb->num_rows() <= 0){
      Rcpp::Rcerr << "ExtractHitsRemote() - arrow::RecordBatch() - No alignments could be computed." << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    arrow::Status align_sts = alignment_rb->ValidateFull();
    if(!align_sts.ok()){
      // Rcpp::Rcout << align_sts.message()  << std::endl << align_sts.ToString() << std::endl << "rows:" << alignment_rb->num_rows() << "\ncols:" << alignment_rb->num_columns()  << std::endl << std::flush; //DEBUG
      Rcpp::Rcerr << "ExtractHitsRemote() - arrow::RecordBatch() - Alignments failed validation." << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    if (alignment_rb)
    {
      // Rcpp::Rcout << "DEBUG: RecordBatch::" << std::endl <<  alignment_rb->ToString() << std::endl << std::flush; //DEBUG
      if(return_values){
        ret_val->emplace_back(alignment_rb);
      }else{
        const auto &wrt_sts = arrow_wrapper->AddRB2Batch(alignment_rb);
        if (!wrt_sts.ok())
        {
          Rcpp::Rcerr << "ExtractHitsRemote() - Error writing RecordBatch..." << std::endl << std::flush; //DEBUG
        }
        ret_val->emplace_back(empty_rb); 
      }
    }else{
      // Rcpp::Rcerr << "[ExtractHitsRemote()] - Empty alignment_rb..." << std::endl << std::flush; //DEBUG
      ret_val->emplace_back(empty_rb);
    }
    
    return ret_val;
  }catch(const ncbi::CException& e) {
    Rcpp::Rcerr << std::string("[ExtractHitsRemote()]: 1. NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
    return std::make_shared<arrow::RecordBatchVector>();
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[ExtractHitsRemote()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush; 
    return std::make_shared<arrow::RecordBatchVector>();
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[ExtractHitsRemote()]: Rcpp Error : ") + e.what() << std::endl << std::flush; 
    return std::make_shared<arrow::RecordBatchVector>();
  }catch(const std::exception &e){
    throw std::runtime_error(std::string("[ExtractHitsRemote()]: C++ Exception : ") + e.what());
  }catch(...){
    throw std::runtime_error("[ExtractHitsRemote()]: Unknown Exception");
  }
}
