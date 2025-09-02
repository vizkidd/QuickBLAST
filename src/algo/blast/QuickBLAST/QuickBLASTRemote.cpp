// remote_blast_integration.cpp
#include <memory>
#include <string>
#include <thread>
#include <chrono>
#include <iostream>
#include <cstdio>

// typical toolkit includes (paths may vary by version)
#include <algo/blast/api/remote_blast.hpp>     // may be remote_blast.hpp or remote.hpp in your tree
#include <algo/blast/api/blast_options_handle.hpp>
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
    const bool return_values = true,
    const unsigned int max_poll_seconds = 120,
    const unsigned int poll_interval_ms = 4000
){
  return pImpl->BLAST_remote(program, database, query_input, input_type, outFile, return_values, max_poll_seconds, poll_interval_ms);
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
    const bool return_values = true,
    const unsigned int max_poll_seconds = 120,
    const unsigned int poll_interval_ms = 4000
){
  assert(out_file.empty() || return_values == true);
  assert(!out_file.empty() || return_values == false);
  
  if(outFile.empty()){
    outFile = std::tmpnam(nullptr); 
  }
  
  // 1) Create options handle (use the same options factory you use for local BLAST)
  // CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(ProgramNameToEnum(program)));
  
  // Optionally set program-specific remote options
  // opts->SetSomeParameter(...);
  
  // 2) Construct the remote object
  //   The actual class and constructor can differ across toolkit versions;
  //   search your toolkit include dir for "RemoteBlast" or "remote_blast".
  CRemoteBlast remote(CBlastOptionsFactory::Create(ProgramNameToEnum(program)));             // adapt if class name is different
  
  // Provide sequences as either FASTA strings or as SSeqLoc objects depending on the API.
  // Many remote helpers accept a simple query string.
  if (!database.empty()) remote.SetDatabase(database);
  
  list<CRef<CSeq_loc>> query_input_list = {};
  
  for(const std::string &query : query_input){
    
    auto q_type = this->arrow_wrapper->CastToType(query);
    if (q_type.header.empty() || q_type.seq.empty()) {
      Rcpp::stop("BLAST_seqs: query header/sequence is empty.");
    }
    
    int rec_no = q_type.rec_no;
    std::string fastaID(q_type.header.data());
    std::string fastaSequence(q_type.seq.data());
    
    const TSeqPos seqlen = fastaSequence.length();
    
    _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
    
    CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
    id->Select(CSeq_id_Base::E_Choice::e_Local);
    id->SetLocal().SetId(rec_no);
    id->SetLocal().SetStr(fastaID);
    
    ENa_strand seqStrand;
    
    if (seq_type == ESeqType::eProtein)
    {
      seqStrand = eNa_strand_unknown;
    }
    else
    {
      switch (strand)
      {
      case EStrand::ePlus:
        seqStrand = ENa_strand::eNa_strand_plus;
        break;
      case EStrand::eMinus:
        seqStrand = ENa_strand::eNa_strand_minus;
        break;
      case EStrand::eUnknown:
        seqStrand = ENa_strand::eNa_strand_unknown;
        break;
      case EStrand::eBoth:
        seqStrand = ENa_strand::eNa_strand_both;
        break;
      case EStrand::eBoth_rev:
        seqStrand = ENa_strand::eNa_strand_both_rev;
        break;
      case EStrand::eOther:
        seqStrand = ENa_strand::eNa_strand_other;
        break;
      }
    }
    
    TSeqPos from = static_cast<TSeqPos>(0);
    TSeqPos to   = static_cast<TSeqPos>(seqlen > 0 ? seqlen - 1 : 0);
    // const CRef<CSeq_loc> query_seqloc(id, from, to, seqStrand);
    CRef<CSeq_loc> query_seqloc(new CSeq_loc());
    
    // Create an interval and attach id, from/to, strand
    CSeq_interval &interval = query_seqloc->SetInt();   // creates/returns interval
    interval.SetId().Assign(*id);                        // copy the seq-id into the interval
    interval.SetFrom(from);
    interval.SetTo(to);
    interval.SetStrand(seqStrand);
    
    query_input_list.push_back(query_seqloc);
  }
  remote.SetQueries(query_input_list);
  // remote.SetSubject(subject_sequence);
  
  // // OPTIONAL: set number of results/hits etc
  // if(max_results > 0)
  //   remote.SetMaxTargetSequences(max_results); // API-dependent
  
  // 3) Submit job and get RID (request id)
  try {
    remote.Submit();                 // may return void or RID; check your version
  }
  catch (const CBlastException &e) {
    // treat as network / submission failure
    throw std::runtime_error(std::string("Remote submit failed: ") + e.what());
  }
  
  // 4) Poll for completion or timeout
  unsigned int waited = 0;
  while (true) {
    // Check status - example API: remote.GetStatus()
    // Many implementations: remote.CheckStatus() or remote.GetRIDStatus()
    CRemoteBlast::ESearchStatus status = remote.CheckStatus(); // adapt to your header
    if (status == CRemoteBlast::ESearchStatus::eStatus_Done) break;
    if (status == CRemoteBlast::ESearchStatus::eStatus_Failed) {
      throw std::runtime_error("Remote BLAST reported an error while processing the job");
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(poll_interval_ms));
    waited += poll_interval_ms/1000;
    if (waited > max_poll_seconds) {
      // you might want to call remote.Cancel() if API provides it
      throw std::runtime_error("Remote BLAST timed out");
    }
  }
  
  // 5) Fetch result as Seq-align(s)
  // Many APIs provide a method like: remote.GetSeqAlignSet() or remote.GetSeqAligns()
  // We want to produce a TSeqAlignVector (same type your ExtractHits expects).
  TSeqAlignVector alignments;
  
  try {
    // The toolkit often returns a CRef<CSeq_align_set> or vector<CRef<CSeq_align_set>>
    // Example:
    CRef<CSeq_align_set> align_set = remote.GetAlignments();
    if (align_set) {
      // empty
      return ExtractHitsRemote(alignments);
    }
    // convert the returned alignment(s) into TSeqAlignVector (toolkit-specific)
    // Some helper exists already in the toolkit; otherwise wrap the single align_set inside vector:
    alignments.emplace_back(align_set);
  }
  catch (const std::exception &e) {
    throw std::runtime_error(std::string("Failed fetching remote results: ") + e.what());
  }
  
  return ExtractHitsRemote(alignments);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::ExtractHitsRemote(const TSeqAlignVector &alignments){
  return pImpl->ExtractHitsRemote(alignments);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::ExtractHitsRemote(const TSeqAlignVector &alignments) 
{
 
 std::shared_ptr<arrow::RecordBatchVector> ret_val = std::make_shared<arrow::RecordBatchVector>();
  
  // quick sanity
  if (alignments.empty()) {
    ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie());
    return ret_val;
  }
 
 
  
  // Arrow builders for the columns (match columns in your schema)
  arrow::StringBuilder qseqid_builder, sseqid_builder;
  arrow::LargeStringBuilder qseq_builder, sseq_builder; // sequence text (likely empty)
  arrow::Int32Builder qstart_builder, qend_builder, sstart_builder, send_builder;
  arrow::Int32Builder aln_len_builder, hsp_num_builder;

  arrow::Int32Builder hsp_offset_builder;
  arrow::Int8Builder length_builder, mismatch_builder, gapopen_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
  arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
  arrow::StringBuilder frames_builder, strand_builder;
  std::string qseq, sseq, frame, strand, qseq_id, sseq_id;
  arrow::Int8Builder qlen_builder, slen_builder, num_alignments_builder;
  
  int num_rows = 0;
  
  for (const auto &align_set_ref : alignments) {
    if (!align_set_ref) continue;
    
    // seq_aligns (list) inside seq_align_set
    const auto &seq_align_list = align_set_ref->Get();
    for (const auto &seq_align : seq_align_list) {
      try
      {
      if (!seq_align) continue;
      
      // coordinates
      // if (seq_align->IsEmpty())
      // {
      //   break;
      // }
      assert(seq_align->IsSet());
      assert(seq_align->CanGet());
      
      // Get seq ids of the two sequences involved in the alignment
      std::string qid = "(unknown)";
      std::string sid = "(unknown)";
      try {
        qid = seq_align->GetSeq_id(0).GetSeqIdString(true);
      } catch (...) { /* ignore — fallback */ }
      try {
        sid = seq_align->GetSeq_id(1).GetSeqIdString(true);
      } catch (...) { /* ignore */ }
      
      ENa_strand q_strand = seq_align->GetSeqStrand(0); // query row
      ENa_strand s_strand = seq_align->GetSeqStrand(1); // subject row
      strand = q_strand + "/" + s_strand;
      
      // assert(!seq_aligns.empty());
      
      // if (seq_aligns.size() > 0) // FILL UP THE ARRAYS
      // {
      
          
      
      // iterate HSPs in this CSeq_align (the CSeq_align may represent one alignment/hsp or have segments)
      // Many toolkit objects treat each CSeq_align as an "hsp", so usually one align_ref is one hit/HSP.
        
        
              
              assert(!seq_align.IsNull());
              if (!seq_align.NotEmpty())
              {
                break;
              }
              
              assert(seq_align->CanGetScore());
              int score, n_splices, num_ident, aln_len, aln_len01, gaps, mismatches, positive, qstart, qend, sstart, send, negative_count;
              double bits, evalue, blast_score, pident, pident_gap, qcovhsp, sum_evalue, product_coverage, overall_identity, high_quality_percent_coverage, exon_identity, consensus_splices, comp_adj_method, matches;
              std::string frames = "*/*";
              
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); 
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident); 
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap); 
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps); 
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); 
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
              seq_align->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp); 
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
              
              qstart = seq_align->GetSeqStart(0);
              qend = seq_align->GetSeqStop(0);
              sstart = seq_align->GetSeqStart(1);
              send = seq_align->GetSeqStop(1);
              
              frames = std::to_string(GetFrame(qstart, aln_len, q_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, s_strand));
              
              static_cast<void>(frames_builder.Append(frames));
              static_cast<void>(qstart_builder.Append(qstart));
              static_cast<void>(qend_builder.Append(qend));
              static_cast<void>(sstart_builder.Append(sstart));
              static_cast<void>(send_builder.Append(send));
              static_cast<void>(pident_builder.Append(pident));
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
              static_cast<void>(qseq_builder.Append(std::string()));
              static_cast<void>(sseq_builder.Append(std::string()));
              static_cast<void>(qlen_builder.Append(qseq.length()));
              static_cast<void>(slen_builder.Append(sseq.length()));
              static_cast<void>(num_alignments_builder.Append(seq_align_list.size()));
              
              static_cast<void>(strand_builder.Append(strand));
              static_cast<void>(hsp_offset_builder.Append(1));
              
              num_rows++;
          
        } catch (const std::exception &e) {
        // best effort: continue to next alignment
        std::cerr << "Warning: exception while processing alignment: " << e.what() << std::endl;
        continue;
      }
    } // end each CSeq_align in set
  } // end each align_set
  
  if (num_rows == 0) {
    ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie());
    return ret_val;
  }
  
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
  
  arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({pident_array,
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
                                                                                                 {"pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});
  
  assert(aln_struct_array.ok());
  
  std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
  
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
  std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int8()), arrow::field("slen", arrow::int8())});
  
  arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
                                                                                               seqids_struct_array,
                                                                                               seqs_struct_array,
                                                                                               strand_array,
                                                                                               lengths_struct_array},
                                                                                               {"num_alignments", "seqids", "seqs", "strands", "lengths"});
  
  assert(seq_info_array.ok());
  
  std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
  
  // Rprintf("\n%d\n", num_rows); //DEBUG
  std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                              num_rows,
                                                                              {seq_info_array_, aln_struct_array_});
  
  if (alignment_rb)
  {
    ret_val->emplace_back(alignment_rb);
  }else{
    ret_val->emplace_back(arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie());
  }
  


  return ret_val;
  
}
