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

#include <algo/blast/QuickBLAST/commons.hpp>
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
    const unsigned int poll_interval_ms = 4000
){
  return pImpl->BLAST_remote(program, database, query_input, input_type, outFile, outFormat, return_values, max_poll_seconds, poll_interval_ms);
}

// wrapper that submits query/subject pair remotely and returns a TSeqAlignVector
// Note: function names used here (CRemoteBlast et al.) may be slightly different in your toolkit version.
// If compile-time names differ grep the toolkit include dir for "Remote" or "remote_blast".
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_remote(
    const std::string &program,         // "blastn", "blastp", "tblastx", etc.
    const std::string &database,        // e.g. "nr" (or leave empty for defaults)
    const Rcpp::List &query_input,  // FASTA header+sequence or just sequence per your wrapper's expectations
    const QuickBLAST::EInputType input_type = EInputType::eSequenceString,
    std::string outFile = "",
    std::string outFormat = "parquet",
    const bool return_values = true,
    const unsigned int max_poll_seconds = 360,
    const unsigned int poll_interval_ms = 4000
){
  assert(out_file.empty() || return_values == true);
  assert(!out_file.empty() || return_values == false);
  
  assert(max_poll_seconds > 0);
  assert(poll_interval_ms > 0);
  
  if(poll_interval_ms < 4000)
    Rcpp::Rcerr << "Warning: poll_interval < 4 seconds might not respect rate limits.";
  
  // if(outFile.empty()){
  //   outFile = std::tmpnam(nullptr); 
  // }
  
  arrow_wrapper->AddFASTAMetadata("BLAST locality", "CBlastOptions::eRemote");
  arrow_wrapper->AddFASTAMetadata("Input source", "sequence");
  arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);
  
  // 1) Create options handle (use the same options factory you use for local BLAST)
  // CRef<CBlastOptionsHandle> opts = CBlastOptionsFactory::Create(ProgramNameToEnum(program),  CBlastOptions::eRemote);
  
  // Optionally set program-specific remote options
  // opts->SetSomeParameter(...);
  
  // 2) Construct the remote object
  //   The actual class and constructor can differ across toolkit versions;
  //   search your toolkit include dir for "RemoteBlast" or "remote_blast".
  // CRemoteBlast remote(CBlastOptionsFactory::Create(ProgramNameToEnum(program)));             // adapt if class name is different
  
  // Provide sequences as either FASTA strings or as SSeqLoc objects depending on the API.
  // Many remote helpers accept a simple query string.
  // if (!database.empty()) remote.SetDatabase(database);
  
  CBlastServices blast_service;
  if(!blast_service.IsValidBlastDb(database, (seq_type == ESeqType::eProtein))){
    Rcpp::stop("BLAST_reomte: Not a valid NCBI database.");
  }
  
  // list<CRef<ncbi::objects::CBioseq>> query_input_list = {}; //list<CRef<CSeq_loc>> 
  CRef<objects::CBioseq_set> bss(new objects::CBioseq_set());
  std::vector<CSeq_entry_Handle> subject_ent_vec;
  bss->SetClass(objects::CBioseq_set::eClass_nuc_prot); //EClass - https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/include/objects/seqset/Bioseq_set_.hpp
  
  CRef<CObjectManager> objMgr(CObjectManager::GetInstance());
  if (!objMgr) {
    Rcpp::stop("BLAST_remote: CObjectManager::GetInstance() returned NULL.");
  }
  CRef<ncbi::CScope> scope(new ncbi::CScope(*objMgr));
  scope->AddDefaults();
  
  for(const std::string &query : query_input){
    Rcpp::checkUserInterrupt();
    auto q_type = this->arrow_wrapper->CastToType(query);
    if (q_type.header.empty() || q_type.seq.empty()) {
      Rcpp::stop("BLAST_remote: query header/sequence is empty.");
    }
    
    assert(!q_type.header.empty());
    assert(!q_type.seq.empty());
    
    int rec_no = q_type.rec_no;
    std::string fastaID(q_type.header.data());
    std::string fastaSequence(q_type.seq.data());
    
    const TSeqPos seqlen = fastaSequence.length();
    
    _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
    
    // CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
    // id->Select(CSeq_id_Base::E_Choice::e_Local);
    // id->SetLocal().SetId(rec_no);
    // // id->SetLocal().SetStr(fastaID);
    // std::string cleanID = fastaID;
    // std::replace(cleanID.begin(), cleanID.end(), ' ', '_');
    // id->SetLocal().SetStr(cleanID);
    
    // ENa_strand seqStrand;
    // 
    // if (seq_type == ESeqType::eProtein)
    // {
    //   seqStrand = eNa_strand_unknown;
    // }
    // else
    // {
    //   switch (strand)
    //   {
    //   case EStrand::ePlus:
    //     seqStrand = ENa_strand::eNa_strand_plus;
    //     break;
    //   case EStrand::eMinus:
    //     seqStrand = ENa_strand::eNa_strand_minus;
    //     break;
    //   case EStrand::eUnknown:
    //     seqStrand = ENa_strand::eNa_strand_unknown;
    //     break;
    //   case EStrand::eBoth:
    //     seqStrand = ENa_strand::eNa_strand_both;
    //     break;
    //   case EStrand::eBoth_rev:
    //     seqStrand = ENa_strand::eNa_strand_both_rev;
    //     break;
    //   case EStrand::eOther:
    //     seqStrand = ENa_strand::eNa_strand_other;
    //     break;
    //   }
    // }
    
    // TSeqPos from = static_cast<TSeqPos>(0);
    // TSeqPos to   = static_cast<TSeqPos>(seqlen > 0 ? seqlen - 1 : 0);
    // // const CRef<CSeq_loc> query_seqloc(id, from, to, seqStrand);
    // CRef<CSeq_loc> query_seqloc(new CSeq_loc());
    // 
    // // Create an interval and attach id, from/to, strand
    // CSeq_interval &interval = query_seqloc->SetInt();   // creates/returns interval
    // interval.SetId().Assign(*id);                        // copy the seq-id into the interval
    // interval.SetFrom(from);
    // interval.SetTo(to);
    // interval.SetStrand(seqStrand);
    
    // std::string cleanID = fastaID;
    // std::replace(cleanID.begin(), cleanID.end(), ' ', '_');
    // if(cleanID.length() > 50)
    //   cleanID.resize(50); //respecting the 50 char limit for FASTA headers by BLAST
    
    // create a seqdesc containing a title string (cleanID)
    // CRef<ncbi::objects::CSeqdesc> title_desc(new ncbi::objects::CSeqdesc());
    // title_desc->SetTitle(fastaID); //cleanID
    // 
    // // put it into a CSeq_descr and attach to the bioseq
    // CRef<ncbi::objects::CSeq_descr> descr(new ncbi::objects::CSeq_descr());
    // descr->Set().push_back(title_desc);
    
    
    CRef<CSeq_id> id(new CSeq_id(fastaID, CSeq_id::fParse_RawText | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal));
    // // id->SetLocal().SetId(rec_no);
    // // id->SetLocal().SetStr(cleanID);
    // id->SetLocal().SetStr(fastaID);
    
    CRef<CSeq_loc>
      cseq_loc_obj(new CSeq_loc());
    cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
    cseq_loc_obj->SetId(*id);
    // cseq_loc_obj->SetWhole()
    //   .SetLocal()
    //   .SetStr(fastaID); //cleanID
    // cseq_loc_obj->SetWhole()
    //   .SetLocal()
    //   .SetId(rec_no);
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
    // other useful defaults:
    seq_inst->SetTopology(CSeq_inst::eTopology_linear);
    
    if (seq_type == ESeqType::eProtein)
    {
      // seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
      seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_aa);
    }
    else
    {
      // seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
      seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_na);
    }
    
    // 
    // // enum EStrand {
    // //   eStrand_not_set =   0,
    // //   eStrand_ss      =   1,  ///< single strand
    // //   eStrand_ds      =   2,  ///< double strand
    // //   eStrand_mixed   =   3,
    // //   eStrand_other   = 255  ///< default ds for DNA, ss for RNA, pept
    // // };
    // 
    // seq_inst->SetRepr(CSeq_inst_Base::ERepr::eRepr_raw);
    // seq_inst->SetTopology(CSeq_inst_Base::ETopology::eTopology_linear);
    seq_inst->SetLength(seqlen);
    // seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
    // 
    // // { CSeq_inst::eMol_not_set, " " } ,
    // // { CSeq_inst::eMol_dna,     "DNA" } ,
    // // { CSeq_inst::eMol_rna,     "RNA" } ,
    // // { CSeq_inst::eMol_aa,      "protein" } ,
    // // { CSeq_inst::eMol_na,      "nucleotide" } ,
    // // { CSeq_inst::eMol_other,   "other" }
    
    // CRef<ncbi::objects::CBioseq> bioseq(new CBioseq(*cseq_loc_obj, cleanID)); //fastaID
    // bioseq->SetInst(*seq_inst);
    // 
    // CRef<CSeq_entry>
    //   ret_entry(new CSeq_entry());
    // ret_entry->SetSeq(*bioseq);
    
    
    CRef<CBioseq> bioseq(new CBioseq());
    bioseq->SetId().push_back(id);      // ID that matches cseq_loc above
    bioseq->SetInst(*seq_inst);
    // bioseq->SetDescr(*descr);
    
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
    
    // // query_input_list.push_back(query_seqloc);
    // query_input_list.push_back(bioseq);
    bss->SetSeq_set().push_back(ret_entry);
  }
  
  // remote.SetQueries(query_input_list);
  // // remote.SetSubject(subject_sequence);
  
  // // OPTIONAL: set number of results/hits etc
  // if(max_results > 0)
  //   remote.SetMaxTargetSequences(max_results); // API-dependent
  
  // 3) Submit job and get RID (request id)
  // try {
  //   remote.Submit();                 // may return void or RID; check your version
  // }
  // catch (const CBlastException &e) {
  //   // treat as network / submission failure
  //   throw std::runtime_error(std::string("Remote submit failed: ") + e.what());
  // }
  // 
  // // 4) Poll for completion or timeout
  // unsigned int waited = 0;
  // while (true) {
  //   // Check status - example API: remote.GetStatus()
  //   // Many implementations: remote.CheckStatus() or remote.GetRIDStatus()
  //   CRemoteBlast::ESearchStatus status = remote.CheckStatus(); // adapt to your header
  //   if (status == CRemoteBlast::ESearchStatus::eStatus_Done) break;
  //   if (status == CRemoteBlast::ESearchStatus::eStatus_Failed) {
  //     throw std::runtime_error("Remote BLAST reported an error while processing the job");
  //   }
  //   std::this_thread::sleep_for(std::chrono::milliseconds(poll_interval_ms));
  //   waited += poll_interval_ms/1000;
  //   if (waited > max_poll_seconds) {
  //     // you might want to call remote.Cancel() if API provides it
  //     throw std::runtime_error("Remote BLAST timed out");
  //   }
  // }
  
  // 3) Build CRemoteBlast
  // CRemoteBlast remote(CBlastOptionsFactory::Create(ProgramNameToEnum(program),  CBlastOptions::eRemote));
  CRemoteBlast remote(SetQuickBLASTOptions(program, GetQuickBLASTOptionString(), CBlastOptions::eRemote));
  remote.SetQueries(bss); //query_input_list
  remote.SetDatabase(database);
  // 4) Submit synchronously, wait, then get results
  // remote.SubmitSync();
  try {
    Rcpp::Rcout << "Max wait time: " << max_poll_seconds << std::endl << std::flush;
    remote.Submit();                 // may return void or RID; check your version
  }
  catch (const CBlastException &e) {
    // treat as network / submission failure
    throw std::runtime_error(std::string("Remote submit failed: ") + e.what());
  }
  // CRemoteBlast::ESearchStatus status = remote.CheckStatus();
  // Progress progress_bar(max_poll_seconds, true);
  // unsigned int waited = 0;
  // while (status == CRemoteBlast::ESearchStatus::eStatus_Pending) {
  //   Rcpp::checkUserInterrupt();
  //   // Check status - example API: remote.GetStatus()
  //   // Many implementations: remote.CheckStatus() or remote.GetRIDStatus()
  //   if (status == CRemoteBlast::ESearchStatus::eStatus_Done) break;
  //   if (status == CRemoteBlast::ESearchStatus::eStatus_Failed) {
  //     vector<std::string> remoteErrors = remote.GetErrorVector();
  //     for(const std::string &error : remoteErrors){
  //       Rcpp::Rcerr << error << std::endl << std::flush; 
  //     }
  //     throw std::runtime_error("Remote BLAST reported an error while processing the job");
  //   }
  //   std::this_thread::sleep_for(std::chrono::milliseconds(poll_interval_ms));
  //   waited += poll_interval_ms/1000;
  //   progress_bar.increment();
  //   if (waited > max_poll_seconds) {
  //     // you might want to call remote.Cancel() if API provides it
  //     throw std::runtime_error("Remote BLAST timed out");
  //   }
  //   status = remote.CheckStatus(); // adapt to your header
  // }
  
  Progress progress_bar(max_poll_seconds, true);
  
  // Use ms accumulator to avoid integer-division artifacts
  unsigned long long waited_ms = 0ULL;
  unsigned int reported_seconds = 0; // how many seconds we've reported to progress
  
  auto status = remote.CheckStatus(); // initial status
  
  while (status == CRemoteBlast::ESearchStatus::eStatus_Pending) {
    Rcpp::checkUserInterrupt();
    
    // If remote changed between checks, break early
    if (status == CRemoteBlast::ESearchStatus::eStatus_Done) break;
    
    if (status == CRemoteBlast::ESearchStatus::eStatus_Failed) {
      std::vector<std::string> remoteErrors = remote.GetErrorVector();
      for (const std::string &error : remoteErrors) {
        Rcpp::Rcerr << error << std::endl << std::flush;
      }
      throw std::runtime_error("Remote BLAST reported an error while processing the job");
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
  
  CRef<CSearchResultSet> results = remote.GetResultSet();
  auto rid = remote.GetRID();
  
  Rcpp::Rcout << "(Success) Remote BLAST run ID: " << rid << std::endl << std::flush; //DEBUG
  
  // 5) Fetch result as Seq-align(s)
  // Many APIs provide a method like: remote.GetSeqAlignSet() or remote.GetSeqAligns()
  // We want to produce a TSeqAlignVector (same type your ExtractHits expects).
  TSeqAlignVector alignments;
  
  try {
    // The toolkit often returns a CRef<CSeq_align_set> or vector<CRef<CSeq_align_set>>
    // Example:
    CRef<CSeq_align_set> align_set = remote.GetAlignments();
    if (!align_set->IsEmpty()) {
      // convert the returned alignment(s) into TSeqAlignVector (toolkit-specific)
      // Some helper exists already in the toolkit; otherwise wrap the single align_set inside vector:
      alignments.emplace_back(align_set);
    }
    
    // CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
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
        // else if (align_ref->GetSegs().IsStd()) {
        //   // handle Std-seg similarly if needed (left as exercise)
        // }
        // // (handle other seg types if your alignments use them)
      }
    }
    
    if (seqids.empty()) {
      Rcpp::stop("RemoteBlast() - Query aligned to no subject IDs. (Subject IDs.size() == 0)");
    }
    
    // 2) call the remote get-sequences service
    CBlastServices::TBioseqVector bioseqs;
    std::string errors, warnings;
    try {
      Rcpp::checkUserInterrupt();
      CBlastServices::GetSequences(seqids, database, seq_type == ESeqType::eProtein ? 'p' : 'n' , bioseqs, errors, warnings, /*verbose=*/false);
    }
    catch (const CException &e) {
      // toolkit exception — handle/log
      // ERR_POST("GetSequences failed: " << e.GetMsg());
      Rcpp::stop(std::string("RemoteBlast() - GetSequences failed: ") + e.GetMsg());
    }
    
    if (!errors.empty()) {
      Rcpp::stop(std::string("RemoteBlast() - GetSequences errors: ") + errors);
    }
    if (!warnings.empty()) {
      Rcpp::stop(std::string("RemoteBlast() - GetSequences warnings: ") + warnings);
    }
    
    // 3) add returned Bioseqs to the scope so CSeqVector/CBioseq_Handle lookups work
    for (auto &bioseq_ref : bioseqs) {
      if (!bioseq_ref) continue;
      Rcpp::checkUserInterrupt();
      CRef<CSeq_entry> entry(new CSeq_entry());
      entry->SetSeq(*bioseq_ref);                 // copy CBioseq into Seq-entry
      scope->AddTopLevelSeqEntry(*entry);         // now scope can resolve those seq-ids
      // // Debug: print the IDs added
      // for (const auto &id : bioseq_ref->GetId()) {
      //   LOG_POST("Loaded subject seq id: " << id->GetSeqIdString(true));
      // }
    }
    
    return ExtractHitsRemote(alignments, subject_ent_vec, *scope, return_values);
    // Rcpp::stop("CSeq_align_set - No alignments could be computed.");
  }
  catch (const std::exception &e) {
    throw std::runtime_error(std::string("Blast_remote(): C++ Exception : Failed fetching remote results: ") + e.what());
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("Blast_remote(): C++ Runtime Error : Failed fetching remote results: ") + e.what() << std::endl << std::flush; 
    return std::make_shared<arrow::RecordBatchVector>();
  }
  
  // return ExtractHitsRemote(alignments);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, CScope &scope, const bool &return_values){
  return pImpl->ExtractHitsRemote(alignments, sseq_entry_vec, scope, return_values);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::ExtractHitsRemote(const TSeqAlignVector &alignments, std::vector<CSeq_entry_Handle>& sseq_entry_vec, CScope &scope, const bool &return_values) 
{
 try{
   std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
    
    // quick sanity
    if (alignments.empty()) {
      Rcpp::stop("TSeqAlignVector - No alignments could be computed.");
      ret_val->emplace_back(empty_rb);
      return ret_val;
    }
    
    CScoreBuilder scorer;
    // if (effective_search_space > 0.0) scorer.SetEffectiveSearchSpace(effective_search_space);
    
    // Compute batch scores (AddScore has an overload for list)
    // We'll ask for a set of scores in a loop to leverage internal batching
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
    
    // Arrow builders for the columns (match columns in your schema)
    std::string qseq = "", sseq = "", frame = "*/*", strand, qseq_id, sseq_id; 
    arrow::StringBuilder qseqid_builder, sseqid_builder, strand_builder; // qseq_title_builder, sseq_title_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder;
    arrow::LargeStringBuilder qhsp_builder, shsp_builder;
    arrow::Int64Builder qlen_builder, slen_builder, num_alignments_builder;
    
    arrow::Int64Builder length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    arrow::StringBuilder frames_builder;
    
    arrow::Int64Builder hsp_offset_builder;
    
    int num_rows = 0;
    
    // CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
    
    for (const auto &align_set_ref : alignments) {
      if (!align_set_ref) continue;
      Rcpp::checkUserInterrupt();
      // seq_aligns (list) inside seq_align_set
      auto &seq_align_list = align_set_ref->Get(); //const
      // for (auto st : score_types) {
      //   try {
      //     scorer.ComputeScore(scope, seq_aligns, st); //scorer.AddScore(scope, seq_align_list, st);
      //   } catch (const CException& e) {
      //     // non-fatal; continue with others
      //     ERR_POST(Warning << "AddScore for type " << static_cast<int>(st) << " failed: " << e.GetMsg());
      //   }
      // }
      for (auto &seq_align : seq_align_list) { //const
        try
        {
        if (!seq_align) continue;
        
        assert(!seq_align.IsNull());
        if (!seq_align.NotEmpty())
        {
          continue;
        }
        
        assert(seq_align->IsSet());
        assert(seq_align->CanGet());
        seq_align->Validate(true);
        RcppThread::checkUserInterrupt();
        // Get seq ids of the two sequences involved in the alignment
        std::string qid = "(unknown)";
        std::string sid = "(unknown)";
        try {
          qid = seq_align->GetSeq_id(0).GetSeqIdString(true);
        } catch (...) { /* ignore — fallback */ }
        try {
          sid = seq_align->GetSeq_id(1).GetSeqIdString(true);
        } catch (...) { /* ignore */ }
        
        // CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(seq_align->GetSeq_id(0));
        // CBioseq_Handle q_bh = scope->GetBioseqHandle(q_idh);
        // const auto q_b = q_bh.GetCompleteObject();
        // const auto qdesc = q_b->GetDescr().Get();
        // std::string qseq_title = qid;
        // for (auto &d : qdesc) {
        //   if (d->IsTitle() && !d->GetTitle().empty()) 
        //     qseq_title = d->GetTitle();
        // }
        // 
        // CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(seq_align->GetSeq_id(1));
        // CBioseq_Handle s_bh = scope->GetBioseqHandle(s_idh);
        // const auto s_b = s_bh.GetCompleteObject();
        // const auto sdesc = s_b->GetDescr().Get();
        // std::string sseq_title = sid;
        // for (auto &d : sdesc) {
        //   if (d->IsTitle() && !d->GetTitle().empty()) 
        //     sseq_title = d->GetTitle();
        // }
        
  
        // const auto& seq_titles = GetTitlesFromSeqAlign(it, scope);
        // std::string qseq_title = seq_titles.first();
        // std::string sseq_title = seq_titles.second();
        
        // scorer.AddSplignScores(seq_align);
        
        ENa_strand q_strand = seq_align->GetSeqStrand(0); // query row
        ENa_strand s_strand = seq_align->GetSeqStrand(1); // subject row
        strand = q_strand + "/" + s_strand;
        
        // assert(!seq_aligns.empty());
        
        // if (seq_aligns.size() > 0) // FILL UP THE ARRAYS
        // {
        
            
        
        // iterate HSPs in this CSeq_align (the CSeq_align may represent one alignment/hsp or have segments)
        // Many toolkit objects treat each CSeq_align as an "hsp", so usually one align_ref is one hit/HSP.
          
        std::string q_full = "", s_full = "";
        std::string q_hsp = "", s_hsp = "", q_aligned = "", s_aligned = "";
        // handle Denseg case
        if (seq_align->GetSegs().IsDenseg()) {
          const CDense_seg& dseg = seq_align->GetSegs().GetDenseg();
          
          // // Get sequence ids (rows)
          // if (dseg.CanGetIds()) {
          //   const auto &ids = dseg.GetIds();
          //   // print/inspect id strings:
          //   for (size_t r = 0; r < ids.size(); ++r) {
          //     if (ids[r]) {
          //       NcbiCout << "Row " << r << " id: " << ids[r]->GetSeqIdString(true) << NcbiEndl;
          //     }
          //   }
          // }
          
          
          switch (save_sequences)
          {
          case true:
            // Full sequences for the two first rows (query, subject)
            if (dseg.CanGetIds()) {
              // try to fetch full sequences for rows 0 and 1
              if (dseg.GetIds().size() > 0) {
                GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[0]), scope, q_full);
              }
              if (dseg.GetIds().size() > 1) {
                GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[1]), scope, s_full);
              }
            }
            qseq = q_full;
            sseq = s_full;
            break;
          }
          
          if(save_hsp_sequences){
            // HSP sequences
            bool ok = GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned);
          }
          // NcbiCout << "Full query length: " << q_full.size() << " HSP ungapped length: " << q_hsp.size() << NcbiEndl;
          // NcbiCout << "Full subject length: " << s_full.size() << " HSP ungapped length: " << s_hsp.size() << NcbiEndl;
          // NcbiCout << "Aligned strings length: " << q_aligned.size() << " / " << s_aligned.size() << NcbiEndl;
          // NcbiCout << "Query HSP: " << q_hsp.substr(0, 200) << NcbiEndl;   // print only prefix
          // NcbiCout << "Subject HSP: " << s_hsp.substr(0, 200) << NcbiEndl;
        }
        // // handle Std-seg (a sequence of local 'loc' entries)
        // else if (seq_align->GetSegs().IsStd()) {
        //   // const CStd_seg &stdseg = seq_align->GetSegs().GetStd();
        //   // stdseg has a list of segments; each segment has a list of locs for each row
        //   // iterate and extract using the loc's intervals
        //   // For brevity, here's a simple approach that attempts to extract by using GetSeqStart/GetSeqStop
        //   int qstart = seq_align->GetSeqStart(0);
        //   int qstop  = seq_align->GetSeqStop(0);
        //   int sstart = seq_align->GetSeqStart(1);
        //   int sstop  = seq_align->GetSeqStop(1);
        //   // fetch sequences by slicing the bioseq handles (if available)
        //   // (You may prefer to iterate stdseg.Get() entries to get exact block-level offsets)
        //   NcbiCout << "Std-seg: q[" << qstart << "," << qstop << "] s[" << sstart << "," << sstop << "]" << NcbiEndl;
        //   // you can reuse GetFullSequenceString + substringing with CSeqVector for exact subrange
        // }
        // else {
        //   // Other seg types: disc, spliced, packed-int, etc.
        //   NcbiCout << "Unhandled seg type; implement specialized extraction if needed" << NcbiEndl;
        // }
                
                assert(seq_align->CanGetScore());
                double score = 0, n_splices = 0, num_ident = 0, aln_len = 0, gaps = 0, mismatches = 0, positive = 0, negative_count = 0;
                double bits = 0, evalue = 0, blast_score = 0, pident = 0, aln_len01 = 0, pident_gap = 0, qcovhsp = 0, sum_evalue = 0, product_coverage = 0, overall_identity = 0, high_quality_percent_coverage = 0, exon_identity = 0, consensus_splices = 0, comp_adj_method = 0, matches = 0;
                std::string frames = "*/*";
                
                bool ok;
                bool haslen = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
                if(!haslen){
                  aln_len = seq_align->GetAlignLength(/*include_gaps*/ true);
                  haslen = true;
                }
                
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
                bool hasid = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
                bool hasp = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident); 
                  
                // compute percent identity fallback per alignment if missing
                if (!hasp && hasid) {
                  double computed = 100.0 * double(num_ident) / seq_align->GetAlignLength(/*include_gaps*/ false); //double(aln_len);
                  // a->SetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, computed);
                  pident = computed;
                  hasp = true;
                }
                  
                bool hasp_gap = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap); 
                if (!hasp_gap && hasid) {
                  double computed = 100.0 * double(num_ident) / seq_align->GetAlignLength(/*include_gaps*/ true); //double(aln_len);
                  // a->SetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, computed);
                  pident_gap = computed;
                  hasp_gap = true;
                }
                
                bool hasgaps = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps); 
                if(!hasgaps){
                  gaps = seq_align->GetTotalGapCount(-1); //seq_align->GetTotalGapCount(0) + seq_align->GetTotalGapCount(1);
                  hasgaps = true;
                }
                
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); 
              
                bool hasmismatches = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
                if(!hasmismatches){
                  mismatches = seq_align->GetAlignLength(/*include_gaps*/ true) - num_ident - gaps;
                  hasmismatches = true;
                }
                
                bool hasqcovhsp = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp); 
                if(!hasqcovhsp){
                  qcovhsp = double(seq_align->GetAlignLength(/*include_gaps*/ false) / q_full.length()); //* 100;
                  hasqcovhsp = true;
                }
                
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices); 
                ok = seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method); 
                
                aln_len01 = seq_align->AlignLengthRatio();
                
                int qstart = seq_align->GetSeqStart(0);
                int qend = seq_align->GetSeqStop(0);
                int sstart = seq_align->GetSeqStart(1);
                int send = seq_align->GetSeqStop(1);
                
                frames = std::to_string(GetFrame(qstart, aln_len, q_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, s_strand));
                
                static_cast<void>(qhsp_builder.Append(q_hsp));
                static_cast<void>(shsp_builder.Append(s_hsp));
                static_cast<void>(frames_builder.Append(frames));
                static_cast<void>(qstart_builder.Append(qstart));
                static_cast<void>(qend_builder.Append(qend));
                static_cast<void>(sstart_builder.Append(sstart));
                static_cast<void>(send_builder.Append(send));
                static_cast<void>(pident_builder.Append(pident)); //pident
                static_cast<void>(evalue_builder.Append(evalue));
                static_cast<void>(length_builder.Append(aln_len));
                static_cast<void>(aln_len01_builder.Append(aln_len01));
                static_cast<void>(bitscore_builder.Append(bits));
                static_cast<void>(score_builder.Append(score));
                static_cast<void>(qcovhsp_builder.Append(qcovhsp));
                static_cast<void>(blast_score_builder.Append(blast_score));
                static_cast<void>(pident_gap_builder.Append(pident_gap));
                static_cast<void>(gaps_builder.Append(gaps));
                static_cast<void>(nident_builder.Append(num_ident));
                static_cast<void>(mismatch_builder.Append(mismatches));
                static_cast<void>(positive_builder.Append(positive));
                static_cast<void>(n_splices_builder.Append(n_splices));
                static_cast<void>(hsp_cnt_builder.Append(num_rows + 1));
                static_cast<void>(sum_evalue_builder.Append(sum_evalue));
                static_cast<void>(product_coverage_builder.Append(product_coverage));
                static_cast<void>(overall_identity_builder.Append(overall_identity));
                static_cast<void>(negative_count_builder.Append(negative_count));
                static_cast<void>(matches_builder.Append(matches));
                static_cast<void>(high_quality_percent_coverage_builder.Append(high_quality_percent_coverage));
                static_cast<void>(exon_identity_builder.Append(exon_identity));
                static_cast<void>(consensus_splices_builder.Append(consensus_splices));
                static_cast<void>(comp_adj_method_builder.Append(comp_adj_method));
                
                /// SEQ INFO
                static_cast<void>(qseqid_builder.Append(qid));
                static_cast<void>(sseqid_builder.Append(sid));
                static_cast<void>(qseq_builder.Append(qseq));
                static_cast<void>(sseq_builder.Append(sseq));
                static_cast<void>(qlen_builder.Append(q_full.length()));
                static_cast<void>(slen_builder.Append(s_full.length()));
                static_cast<void>(num_alignments_builder.Append(seq_align_list.size()));
                
                static_cast<void>(strand_builder.Append(strand));
                static_cast<void>(hsp_offset_builder.Append(1));
    
                // static_cast<void>(qseq_title_builder.Append(qseq_title));
                // static_cast<void>(sseq_title_builder.Append(sseq_title));
                
                num_rows++;
            
          } catch (const std::exception &e) {
          // best effort: continue to next alignment
          std::cerr << "ExtractHitsRemote(): Warning: exception while processing alignment: " << e.what() << std::endl;
          continue;
        }
      } // end each CSeq_align in set
    } // end each align_set
    
    if (num_rows == 0) {
      ret_val->emplace_back(empty_rb);
      return ret_val;
    }
    
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
    
    assert(aln_struct_array.ok());
    
    std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
    
    // std::shared_ptr<arrow::Array> qseq_title_array;
    // static_cast<void>(qseq_title_builder.Finish(&qseq_title_array));
    // 
    // std::shared_ptr<arrow::Array> sseq_title_array;
    // static_cast<void>(sseq_title_builder.Finish(&sseq_title_array));
    
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
    // std::shared_ptr<arrow::StructArray> seqtitle_struct_array = *arrow::StructArray::Make({qseq_title_array, sseq_title_array}, {arrow::field("qseq_title", arrow::utf8()), arrow::field("sseq_title", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
  seqids_struct_array,
  seqs_struct_array,
  strand_array,                                                                                               lengths_struct_array},
  {"num_alignments", "seqids", "seqs", "strands", "lengths"});
    
    assert(seq_info_array.ok());
    
    std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
    
    // Rprintf("\n%d\n", num_rows); //DEBUG
    std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                                num_rows,
                                                                                {seq_info_array_, aln_struct_array_});
    
    for(auto s_ent: sseq_entry_vec){
      scope.RemoveTopLevelSeqEntry(s_ent);
      // s_ent.Reset();
    }
    sseq_entry_vec.clear();
    sseq_entry_vec.shrink_to_fit();
    
    if(alignment_rb->num_rows() <= 0){
      Rcpp::stop("ExtractHitsRemote() - arrow::RecordBatch() - No alignments could be computed.");
    }
    
    arrow::Status align_sts = alignment_rb->ValidateFull();
    if(!align_sts.ok()){
      // Rcpp::Rcout << align_sts.message()  << std::endl << align_sts.ToString() << std::endl << "rows:" << alignment_rb->num_rows() << "\ncols:" << alignment_rb->num_columns()  << std::endl << std::flush; //DEBUG
      Rcpp::stop("ExtractHitsRemote() - arrow::RecordBatch() - Alignments failed validation.");
    }
    
    if (alignment_rb)
    {
      // Rcpp::Rcout << "DEBUG: RecordBatch::" << std::endl <<  alignment_rb->ToString() << std::endl << std::flush; //DEBUG
      // if(save_sequences){
        if(return_values){
          ret_val->emplace_back(alignment_rb);
        }else{
          const auto &wrt_sts = arrow_wrapper->AddRB2Batch(alignment_rb);
          if (!wrt_sts.ok())
          {
            // ret_val->emplace_back(alignment_rb);
            Rcpp::Rcerr << "ExtractHitsRemote() - Error writing RecordBatch..." << std::endl << std::flush; //DEBUG 
          }
          ret_val->emplace_back(empty_rb); 
        // }else{
        //   ret_val->emplace_back(alignment_rb);
        // }
        // ret_val->emplace_back(alignment_rb);
        }
    }else{
      Rcpp::Rcerr << "ExtractHitsRemote() - Empty alignment_rb..." << std::endl << std::flush; //DEBUG
      ret_val->emplace_back(empty_rb);
    }
    
    return ret_val;
 }
 catch(const std::exception &e){
   Rcpp::stop(std::string("ExtractHitsRemote(): C++ Exception : ") + e.what());
 }
 catch(const std::runtime_error &e){
   Rcpp::stop(std::string("ExtractHitsRemote(): C++ Runtime Error : ") + e.what());
 }
 catch(...){
   Rcpp::stop("ExtractHitsRemote(): Unknown Exception");
 }
}
