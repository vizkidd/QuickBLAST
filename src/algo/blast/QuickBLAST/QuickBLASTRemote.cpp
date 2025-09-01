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

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_remote(
    const std::string &program,         // "blastn", "blastp", "tblastx", etc.
    const std::string &database,        // e.g. "nr" (or leave empty for defaults)
    const std::string &query_input,  // FASTA header+sequence or just sequence per your wrapper's expectations
    std::string &outFile = "",
    const bool return_values = true,
    const std::string &APIKey = "",
    const unsigned int max_results = 0,
    const unsigned int max_poll_seconds = 120,
    const unsigned int poll_interval_ms = 4000){
  return pImpl->(program, database, query_input, input_type, outFile, return_values, APIKey, max_results, max_poll_seconds, poll_interval_ms);
}

// wrapper that submits query/subject pair remotely and returns a TSeqAlignVector
// Note: function names used here (CRemoteBlast et al.) may be slightly different in your toolkit version.
// If compile-time names differ grep the toolkit include dir for "Remote" or "remote_blast".
std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::BLAST_remote(
    const std::string &program,         // "blastn", "blastp", "tblastx", etc.
    const std::string &database,        // e.g. "nr" (or leave empty for defaults)
    const std::string &query_input,  // FASTA header+sequence or just sequence per your wrapper's expectations
    std::string &outFile = "",
    const bool return_values = true,
    const std::string &APIKey = "",
    const unsigned int max_results = 0,
    const unsigned int max_poll_seconds = 120,
    const unsigned int poll_interval_ms = 4000)
{
  assert(out_file.empty() || return_values == true);
  assert(!out_file.empty() || return_values == false);
  
  if(outFile.empty()){
    outFile = std::tmpnam(nullptr); 
  }
  
  // 1) Create options handle (use the same options factory you use for local BLAST)
  CRef<CBlastOptionsHandle> opts = CBlastOptionsFactory::Create(program);
  
  // Optionally set program-specific remote options
  // opts->SetSomeParameter(...);
  
  // 2) Construct the remote object
  //   The actual class and constructor can differ across toolkit versions;
  //   search your toolkit include dir for "RemoteBlast" or "remote_blast".
  CRemoteBlast remote(*opts);             // adapt if class name is different
  
  // If you want to set Entrez API key (raises rate limits):
  if(!APIKey.empty())
    remote.SetEntrezKey(APIKey); 
  
  // Provide sequences as either FASTA strings or as SSeqLoc objects depending on the API.
  // Many remote helpers accept a simple query string.
  remote.SetProgram(program);
  if (!database.empty()) remote.SetDatabase(database);
  remote.SetQuery(query_input);
  // remote.SetSubject(subject_sequence);
  
  // OPTIONAL: set number of results/hits etc
  if(max_results > 0)
    remote.SetMaxTargetSequences(max_results); // API-dependent
  
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
    ERemoteBlastStatus status = remote.GetStatus(); // adapt to your header
    if (status == eRemoteBlastStatus_Done) break;
    if (status == eRemoteBlastStatus_Error) {
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
    CRef<CSeq_align_set> align_set = remote.GetSeqAlignSet();
    if (align_set && !align_set->Get()) {
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

std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractHitsRemote(const TSeqAlignVector &alignments){
  return pImpl->ExtractHitsRemote(alignments);
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractHitsRemote(const TSeqAlignVector &alignments)
{
  // quick checks
  if (!arrow_wrapper) {
    Rcpp::Rcerr << "ExtractHitsRemote: arrow_wrapper not initialized\n";
    return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
  }
  if (alignments.empty()) {
    // return an empty but typed record batch
    return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
  }
  
  try {
    // Builders for columns (match the schema used elsewhere)
    arrow::StringBuilder qseqid_builder, sseqid_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder; // left empty in remote mode
    arrow::Int32Builder qlen_builder, slen_builder, num_alignments_builder;
    
    arrow::StringBuilder frames_builder;
    arrow::Int32Builder qstart_builder, qend_builder, sstart_builder, send_builder, hsp_offset_builder;
    arrow::Int32Builder length_builder, hsp_cnt_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder;
    arrow::DoubleBuilder qcovhsp_builder, blast_score_builder;
    arrow::Int32Builder gaps_builder, nident_builder, mismatch_builder, positive_builder, n_splices_builder;
    arrow::Int32Builder negative_count_builder, matches_builder;
    arrow::DoubleBuilder sum_evalue_builder, product_coverage_builder, overall_identity_builder, high_quality_percent_coverage_builder;
    arrow::DoubleBuilder exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    // Add other builders you need (names taken from your previous schema)
    
    int total_rows = 0;
    
    // iterate over alignment sets
    for (const auto &align_set_ref : alignments) {
      if (!align_set_ref) continue;
      if (!align_set_ref->IsSet()) continue;
      
      const auto &vec = align_set_ref->Get();
      if (vec.empty()) continue;
      
      // each element in vec is a CRef<CSeq_align>
      for (const auto &aln_ref : vec) {
        if (!aln_ref) continue; //|| aln_ref->IsEmpty()
        
        // extract common fields; many of these calls exist across toolkit versions;
        // if your toolkit uses other method names, adapt accordingly.
        try {
          // seq ids (row 0 = query, row 1 = subject)
          std::string qid = "NA", sid = "NA";
          try {
            if (aln_ref->GetSegs().Which() == CSeq_align::TSegs::e_Denseg) {
              // Typically GetSeq_id is available: get id for row 0 and 1 when possible
            }
          } catch(...) { /* ignore */ }
          
          // toolkit has method GetSeq_id(row) in many versions - attempt and fallback
          try {
            const CSeq_id_Handle qidh = aln_ref->GetSeq_id(0); // might be different type in some versions
            qid = qidh.GetSeqIdString(true);
          } catch (...) {
            // fallback: try via GetSeq_id() on the alignment segments/locs if available
            try {
              qid = aln_ref->GetSeq_id(0)->GetSeqIdString(true);
            } catch (...) { qid = "NA"; }
          }
          
          try {
            const CSeq_id_Handle sidh = aln_ref->GetSeq_id(1);
            sid = sidh.GetSeqIdString(true);
          } catch (...) {
            try {
              sid = aln_ref->GetSeq_id(1)->GetSeqIdString(true);
            } catch (...) { sid = "NA"; }
          }
          
          // alignment coordinates and scores
          int qstart = 0, qend = 0, sstart = 0, send = 0;
          double bits = 0.0, evalue = 0.0, blast_score = 0.0, pident = 0.0, pident_gap = 0.0;
          int aln_len = 0, gaps = 0, num_ident = 0, mismatches = 0, positive = 0, n_splices = 0;
          double sum_evalue = 0.0, product_coverage = 0.0, overall_identity = 0.0;
          int negative_count = 0, matches = 0;
          double high_quality_percent_coverage = 0.0, exon_identity = 0.0, consensus_splices = 0.0, comp_adj_method = 0.0;
          
          // many toolkits provide GetNamedScore; wrap calls in try/catch to avoid exceptions
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices); } catch(...) {}
          try { aln_ref->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method); } catch(...) {}
          
          // get coordinates (these calls should be available)
          try { qstart = aln_ref->GetSeqStart(0); } catch(...) { qstart = 0; }
          try { qend   = aln_ref->GetSeqStop(0);  } catch(...) { qend = 0; }
          try { sstart = aln_ref->GetSeqStart(1); } catch(...) { sstart = 0; }
          try { send   = aln_ref->GetSeqStop(1);  } catch(...) { send = 0; }
          
          // compute alignment length
          try {
            if (aln_ref->CanGetSegs() && aln_ref->GetSegs().Which() == CSeq_align::TSegs::e_Denseg) {
              // denseg has length info; fallback to using qend-qstart+1
              aln_len = (qend - qstart) + 1;
            } else {
              aln_len = (qend - qstart) + 1;
            }
          } catch(...) {
            aln_len = (qend - qstart) + 1;
          }
          
          // frame: remote queries may not carry explicit strand info. create a placeholder
          std::string frames = "*/*";
          // if you have frame logic use qstart/aln_len/etc to compute as you did before
          
          // append to builders
          static_cast<void>(qseqid_builder.Append(qid));
          static_cast<void>(sseqid_builder.Append(sid));
          static_cast<void>(qseq_builder.Append("")); // no sequence text in remote mode
          static_cast<void>(sseq_builder.Append("")); // no sequence text in remote mode
          static_cast<void>(qlen_builder.Append((int) (qend - qstart + 1)));
          static_cast<void>(slen_builder.Append((int) (send - sstart + 1)));
          static_cast<void>(num_alignments_builder.Append(1)); // treat each alignment as 1
          static_cast<void>(frames_builder.Append(frames));
          static_cast<void>(qstart_builder.Append(qstart));
          static_cast<void>(qend_builder.Append(qend));
          static_cast<void>(sstart_builder.Append(sstart));
          static_cast<void>(send_builder.Append(send));
          static_cast<void>(length_builder.Append(aln_len));
          static_cast<void>(pident_builder.Append(pident));
          static_cast<void>(pident_gap_builder.Append(pident_gap));
          static_cast<void>(evalue_builder.Append(evalue));
          static_cast<void>(bitscore_builder.Append(bits));
          static_cast<void>(score_builder.Append(0.0)); // if you can't get score, use 0
          static_cast<void>(qcovhsp_builder.Append(0.0)); // remote may not contain
          static_cast<void>(blast_score_builder.Append(blast_score));
          static_cast<void>(gaps_builder.Append(gaps));
          static_cast<void>(nident_builder.Append(num_ident));
          static_cast<void>(mismatch_builder.Append(mismatches));
          static_cast<void>(positive_builder.Append(positive));
          static_cast<void>(n_splices_builder.Append(n_splices));
          static_cast<void>(hsp_cnt_builder.Append(total_rows + 1));
          static_cast<void>(sum_evalue_builder.Append(sum_evalue));
          static_cast<void>(product_coverage_builder.Append(product_coverage));
          static_cast<void>(overall_identity_builder.Append(overall_identity));
          static_cast<void>(negative_count_builder.Append(negative_count));
          static_cast<void>(matches_builder.Append(matches));
          static_cast<void>(high_quality_percent_coverage_builder.Append(high_quality_percent_coverage));
          static_cast<void>(exon_identity_builder.Append(exon_identity));
          static_cast<void>(consensus_splices_builder.Append(consensus_splices));
          static_cast<void>(comp_adj_method_builder.Append(comp_adj_method));
          
          // increment row counter
          ++total_rows;
        } catch (const std::exception &e) {
          // log and continue with next alignment
          Rcpp::Rcerr << "Warning ExtractHitsRemote: exception while parsing alignment: " << e.what() << "\n";
          continue;
        }
      } // each aln_ref
    } // each align_set_ref
    
    // If no rows created, return empty RB
    if (total_rows == 0) {
      return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
    }
    
    // Finish builders into arrays
    std::shared_ptr<arrow::Array> qseqid_array, sseqid_array, qseq_array, sseq_array;
    std::shared_ptr<arrow::Array> qlen_array, slen_array, num_alignments_array;
    std::shared_ptr<arrow::Array> frames_array;
    std::shared_ptr<arrow::Array> qstart_array, qend_array, sstart_array, send_array, length_array;
    std::shared_ptr<arrow::Array> pident_array, pident_gap_array, evalue_array, bitscore_array, score_array;
    std::shared_ptr<arrow::Array> qcovhsp_array, blast_score_array, gaps_array, nident_array, mismatch_array;
    std::shared_ptr<arrow::Array> positive_array, n_splices_array, hsp_cnt_array;
    std::shared_ptr<arrow::Array> sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array, matches_array;
    std::shared_ptr<arrow::Array> high_quality_percent_coverage_array, exon_identity_array, consensus_splices_array, comp_adj_method_array;
    
    // finish (ignore status details for brevity; you can check them and handle errors)
    frames_builder.Finish(&frames_array);
    qseqid_builder.Finish(&qseqid_array);
    sseqid_builder.Finish(&sseqid_array);
    qseq_builder.Finish(&qseq_array);
    sseq_builder.Finish(&sseq_array);
    qlen_builder.Finish(&qlen_array);
    slen_builder.Finish(&slen_array);
    num_alignments_builder.Finish(&num_alignments_array);
    qstart_builder.Finish(&qstart_array);
    qend_builder.Finish(&qend_array);
    sstart_builder.Finish(&sstart_array);
    send_builder.Finish(&send_array);
    length_builder.Finish(&length_array);
    pident_builder.Finish(&pident_array);
    pident_gap_builder.Finish(&pident_gap_array);
    evalue_builder.Finish(&evalue_array);
    bitscore_builder.Finish(&bitscore_array);
    score_builder.Finish(&score_array);
    qcovhsp_builder.Finish(&qcovhsp_array);
    blast_score_builder.Finish(&blast_score_array);
    gaps_builder.Finish(&gaps_array);
    nident_builder.Finish(&nident_array);
    mismatch_builder.Finish(&mismatch_array);
    positive_builder.Finish(&positive_array);
    n_splices_builder.Finish(&n_splices_array);
    hsp_cnt_builder.Finish(&hsp_cnt_array);
    sum_evalue_builder.Finish(&sum_evalue_array);
    product_coverage_builder.Finish(&product_coverage_array);
    overall_identity_builder.Finish(&overall_identity_array);
    negative_count_builder.Finish(&negative_count_array);
    matches_builder.Finish(&matches_array);
    high_quality_percent_coverage_builder.Finish(&high_quality_percent_coverage_array);
    exon_identity_builder.Finish(&exon_identity_array);
    consensus_splices_builder.Finish(&consensus_splices_array);
    comp_adj_method_builder.Finish(&comp_adj_method_array);
    
    // Construct struct arrays if your schema expects struct columns (seqs, seqids),
    // otherwise pass arrays in the same order your schema requires.
    // Here we create two struct arrays similar to your previous code:
    auto seqids_struct = *arrow::StructArray::Make({qseqid_array, sseqid_array},
    {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    auto seqs_struct = *arrow::StructArray::Make({qseq_array, sseq_array},
    {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    
    // Create the final RecordBatch using your BLAST schema and columns in order.
    std::shared_ptr<arrow::RecordBatch> rb =
      arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                               total_rows,
                               { seqids_struct,
                                 // the aln struct - create using previously finished arrays
                                 *arrow::StructArray::Make({
                                   pident_array, pident_gap_array, frames_array, evalue_array, length_array,
                                   /* length01 */ length_array, qstart_array, qend_array, sstart_array, send_array,
                                   bitscore_array, score_array, qcovhsp_array, blast_score_array, gaps_array,
                                   nident_array, mismatch_array, positive_array, n_splices_array, hsp_cnt_array,
                                   sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array,
                                   matches_array, high_quality_percent_coverage_array, exon_identity_array,
                                   consensus_splices_array, comp_adj_method_array },
                                   { arrow::field("pident", arrow::float64()), /* ... remaining field definitions ... */ } )
                                 // you may need to pass more fields/structs in the same order as your schema
                               });
    
    return rb;
  }
  catch (const std::exception &e) {
    Rcpp::Rcerr << "ExtractHitsRemote: caught exception: " << e.what() << "\n";
    return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
  }
}

