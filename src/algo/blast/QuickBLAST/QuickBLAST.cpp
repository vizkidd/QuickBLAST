// #ifdef WIN32

// //#include <algo/blast/QuickBLAST/quick_blast.hpp>

// /* template <typename OptionsType>
// ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
// {
//   this->blast_options = options;
//   this->program = program_name;
//   if constexpr (std::is_same_v<OptionsType, std::string> || std::is_same_v<OptionsType, Rcpp::String>)
//   {
//     // this->blast_options_str = options;
//     return SetQuickBLASTOptions<std::string>(program_name, options);
//   }
//   else if constexpr (std::is_same_v<OptionsType, Rcpp::List>)
//   {
//     // this->blast_options_list = options;
//     return SetQuickBLASTOptions<Rcpp::List>(program_name, options);
//   }
//   else
//   {
//     static_assert(std::is_same_v<OptionsType, OptionsType>, "Unsupported type");
//   }
// }*/

// void QuickBLAST::SetThreadCount(int num_threads)
// {
//     this->num_threads = num_threads;
// #if defined(_OPENMP) || defined(WIN32)
//     omp_set_num_threads(num_threads);
// #endif
//     arrow_wrapper->SetThreadCount(num_threads);
// }
// int QuickBLAST::GetHitCount()
// {
//     return hit_count;
// }
// void QuickBLAST::AddHitCount(int val)
// {
// #if defined(_OPENMP) || defined(WIN32)
//     omp_set_lock(&hit_countLock);
// #endif
//     hit_count += val;
// #if defined(_OPENMP) || defined(WIN32)
//     omp_unset_lock(&hit_countLock);
// #endif
// }
// ncbi::blast::CBlastOptionsHandle& QuickBLAST::GetQuickBLASTOptions()
// {
//     return *opts;
// }

// // template <>
// ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
// {
//   assert(!program_name.empty());

//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options.empty())
//   {
//     cout << "Using " << program_name << " Defaults..." << std::endl;
//     return opts;
//   }

//   std::vector<std::pair<std::string, std::string>> keyValuePairs = BLASTOptionsFromString(options);

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (const auto &pair : keyValuePairs)
//   {
//     std::string key_str = pair.second;
//     std::size_t key = std::hash<std::string>{}(pair.first);

//     if (key == hashMap["evalue"])
//     {
//       double val = std::stod(key_str);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = std::stod(key_str);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = (key_str == "TRUE" || key_str == "True" || key_str == "true" || key_str == "1");
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = key_str;
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = std::stod(key_str);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = std::stod(key_str);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = std::stod(key_str);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = std::stoi(key_str);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// }

// /* template <>
// ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
// {
//   assert(!program_name.empty());
//   Rcpp::List options_(options);
//   for (int i = 0; i < options_.size(); ++i)
//   {
//     // Check if the element is empty
//     if (Rf_isNull(options_[i]))
//     {
//       // Remove the empty element
//       options_.erase(i);
//       --i; // Decrement the index since the list size has changed
//     }
//   }
//   ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
//   // Create a CBlastOptionsHandle object
//   ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
//   opts->SetDefaults();
//   // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
//   // Example: Extracting and setting the BLAST database

//   if (options_.size() == 0 || options_.isNULL())
//   {
//     Rcpp::Rcout << "Using " << program_name << " Defaults..." << std::endl;
//     return opts;
//   }

//   Rcpp::CharacterVector keys = options_.names();

//   std::unordered_map<std::string, std::size_t> hashMap;

//   hashMap["evalue"] = std::hash<std::string>{}("evalue");
//   hashMap["pident"] = std::hash<std::string>{}("pident");
//   hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
//   hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
//   hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
//   hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
//   hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
//   hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
//   hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
//   hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
//   hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
//   hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
//   hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
//   hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
//   hashMap["window_size"] = std::hash<std::string>{}("window_size");

//   for (int i = 0; i < keys.size(); ++i)
//   {
//     std::string key_str = Rcpp::as<std::string>(keys[i]);
//     std::size_t key = std::hash<std::string>{}(key_str);

//     if (key == hashMap["evalue"])
//     {
//       double val = Rcpp::as<double>(options_["evalue"]);
//       opts->SetEvalueThreshold(val);
//     }
//     else if (key == hashMap["pident"])
//     {
//       double val = Rcpp::as<double>(options_["pident"]);
//       opts->SetPercentIdentity(val);
//     }
//     else if (key == hashMap["gapped_mode"])
//     {
//       bool val = Rcpp::as<bool>(options_["gapped_mode"]);
//       opts->SetGappedMode(val);
//     }
//     else if (key == hashMap["filter_string"])
//     {
//       std::string val = Rcpp::as<std::string>(options_["filter_string"]);
//       opts->SetFilterString(val.c_str());
//     }
//     else if (key == hashMap["effective_search_space"])
//     {
//       int val = Rcpp::as<int>(options_["effective_search_space"]);
//       opts->SetEffectiveSearchSpace(val);
//     }
//     else if (key == hashMap["cutoff_score"])
//     {
//       int val = Rcpp::as<int>(options_["cutoff_score"]);
//       opts->SetCutoffScore(val);
//     }
//     else if (key == hashMap["gap_trigger"])
//     {
//       double val = Rcpp::as<double>(options_["gap_trigger"]);
//       opts->SetGapTrigger(val);
//     }
//     else if (key == hashMap["gap_x_dropoff"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff"]);
//       opts->SetGapXDropoff(val);
//     }
//     else if (key == hashMap["gap_x_dropoff_final"])
//     {
//       double val = Rcpp::as<double>(options_["gap_x_dropoff_final"]);
//       opts->SetGapXDropoffFinal(val);
//     }
//     else if (key == hashMap["hit_list_size"])
//     {
//       int val = Rcpp::as<int>(options_["hit_list_size"]);
//       opts->SetHitlistSize(val);
//     }
//     else if (key == hashMap["low_score_percentage"])
//     {
//       double val = Rcpp::as<double>(options_["low_score_percentage"]);
//       opts->SetLowScorePerc(val);
//     }
//     else if (key == hashMap["max_hsp_per_subject"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_subject"]);
//       opts->SetMaxHspsPerSubject(val);
//     }
//     else if (key == hashMap["max_hsp_per_sequence"])
//     {
//       int val = Rcpp::as<int>(options_["max_hsp_per_sequence"]);
//       opts->SetMaxNumHspPerSequence(val);
//     }
//     else if (key == hashMap["qcovhsp_perc"])
//     {
//       double val = Rcpp::as<double>(options_["qcovhsp_perc"]);
//       opts->SetQueryCovHspPerc(val);
//     }
//     else if (key == hashMap["window_size"])
//     {
//       int val = Rcpp::as<int>(options_["window_size"]);
//       opts->SetWindowSize(val);
//     }
//   }

//   opts->Validate();

//   return opts;
// } */

// QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences)
// {

// #if defined(_OPENMP) || defined(WIN32)
//   this->num_threads = omp_get_num_threads();
// #else
//   this->num_threads = 1;
// #endif
//   arrow_wrapper = std::make_shared<ArrowWrapper>();
//   this->save_sequences = save_sequences;
//   this->program = program;
//   //this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<std::string>(program, options));
//   this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
//   this->strand = strand;
//   this->seq_type = seq_type;
//   ok_promise.set_value(arrow::Status::OK());
// #if defined(_OPENMP) || defined(WIN32)
//   omp_init_lock(&hit_countLock);
// #endif

//   this->obj_id = (int)obj_list.size();

// }

// QuickBLAST::~QuickBLAST()
// {
// #if defined(_OPENMP) || defined(WIN32)
//   omp_destroy_lock(&hit_countLock);
// #endif

//   // DO NOT DELETE NCBI C++ OBJECTs or PTRs or face Corruption
//   //  delete self;
//   //  opts->ReleaseReference();
//   // delete opts;

//   // delete arrow_wrapper;

//   cout << "~QuickBLAST " << std::endl;
// }

// // Function to process a single FASTA block
// void QuickBLAST::PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream)
// {
//   if (outputStream != nullptr)
//   {

//     // Print FastaSequenceData object
//     (*outputStream) << "No: " << data->rec_no << std::endl;
//     (*outputStream) << "Header: " << data->header << std::endl;
//     (*outputStream) << "Sequence: " << data->seq << std::endl;
//     (*outputStream) << std::endl;
//     outputStream->flush();
//   }
// }

// template <typename T>
// SSeqLoc *QuickBLAST::CreateSSeqLocFromType(T _t, CRef<ncbi::CScope> parent_scope)
// {
//   return CreateSSeqLocFromType<FastaSequenceData>(arrow_wrapper->CastToType<FastaSequenceData>(_t), parent_scope);
// }

// template <>
// SSeqLoc *QuickBLAST::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope)
// {
//   int rec_no = fasta_data.rec_no;
//   std::string fastaID(fasta_data.header.data());
//   std::string fastaSequence(fasta_data.seq.data());

//   const TSeqPos seqlen = fastaSequence.length();

//   _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
//   ncbi::CRef<ncbi::objects::CSeq_interval> interval(new ncbi::objects::CSeq_interval());
//   interval->SetFrom(0);
//   interval->SetTo(seqlen - 1);

//   CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
//   id->Select(CSeq_id_Base::E_Choice::e_Local);
//   id->SetLocal().SetId(rec_no);
//   id->SetLocal().SetStr(fastaID);

//   CRef<CSeq_loc>
//       cseq_loc_obj(new CSeq_loc());
//   cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
//   cseq_loc_obj->SetWhole()
//       .SetLocal()
//       .SetStr(fastaID);
//   if (seq_type == ESeqType::eProtein)
//   {
//     cseq_loc_obj->SetStrand(eNa_strand_unknown);
//   }
//   else
//   {
//     switch (strand)
//     {
//     case EStrand::ePlus:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_plus);
//       break;
//     case EStrand::eMinus:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_minus);
//       break;
//     case EStrand::eUnknown:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_unknown);
//       break;
//     case EStrand::eBoth:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both);
//       break;
//     case EStrand::eBoth_rev:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_both_rev);
//       break;
//     case EStrand::eOther:
//       cseq_loc_obj->SetStrand(ENa_strand::eNa_strand_other);
//       break;
//     }
//   }

//   CRef<CSeq_data> seq_data(new CSeq_data());
//   seq_data->Select(seq_type == ESeqType::eProtein ? CSeq_data_Base::E_Choice::e_Iupacaa : CSeq_data_Base::E_Choice::e_Iupacna);
//   switch (seq_type)
//   {
//   case ESeqType::eProtein:
//   {
//     seq_data->SetIupacaa(CIUPACaa(fastaSequence));
//     break;
//   }
//   case ESeqType::eNucleotide:
//   {
//     seq_data->SetIupacna(CIUPACna(fastaSequence));
//     break;
//   }
//   }

//   CRef<CSeq_inst> seq_inst(new CSeq_inst());
//   seq_inst->SetSeq_data(*seq_data);
//   seq_inst->SetLength(fastaSequence.length());
//   seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
//   seq_inst->SetTopology(CSeq_inst::eTopology_linear);
//   seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
//   seq_inst->SetRepr(CSeq_inst_Base::ERepr::eRepr_raw);
//   seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
//   seq_inst->SetTopology(CSeq_inst_Base::ETopology::eTopology_linear);
//   seq_inst->SetStrand(CSeq_inst_Base::EStrand::eStrand_ss);
//   seq_inst->SetLength(seqlen);

//   CRef<ncbi::objects::CBioseq> bioseq(new CBioseq(*cseq_loc_obj, fastaID));
//   bioseq->SetInst(*seq_inst);

//   CRef<CSeq_entry>
//       ret_entry(new CSeq_entry());
//   ret_entry->SetSeq(*bioseq);

//   parent_scope->AddTopLevelSeqEntry(*ret_entry);

//   return new SSeqLoc(cseq_loc_obj.GetObject(), parent_scope.GetObject());
// }

// int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
// {
//   int frame = 0;
//   if (strand == eNa_strand_plus)
//   {
//     frame = (start % 3) + 1;
//   }
//   else if (strand == eNa_strand_minus)
//   {
//     frame = -(((int)length - start - 1) % 3 + 1);
//   }
//   return frame;
// }

// template <typename T>
// std::conditional_t<std::is_same_v<T, TSeqLocVector>, std::shared_ptr<arrow::RecordBatchVector>, std::shared_ptr<arrow::RecordBatch>> QuickBLAST::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const T &sloc, const CScope &scope)
// {
//   if constexpr (std::is_same_v<T, TSeqLocVector>)
//   {
//     arrow::RecordBatchVector recBth_vec;

//     for (const auto &s_it : sloc)
//     {
//       // Rcpp::checkUserInterrupt();

//       std::shared_ptr<arrow::RecordBatch> rb = ExtractHits<SSeqLoc>(alignments, qloc, s_it, scope);

//       if (rb)
//       {
//         recBth_vec.emplace_back(std::move(rb));
//       }
//     }

//     const auto &wrt_sts = arrow_wrapper->AddRBV2Batch(recBth_vec);
//     if (wrt_sts.ok())
//     {
//       return wrt_sts.ValueOrDie();
//     }
//     else
//     {
//       cerr << "Warn : Invalid Alignment RBV (Returning Empty) : " << wrt_sts.status().detail() << std::endl
//                   << wrt_sts.status().message() << std::endl;
//     }
//     return std::make_shared<arrow::RecordBatchVector>();
//   }
//   // For SSeqLoc
//   else if constexpr (std::is_same_v<T, SSeqLoc>)
//   {
//     std::shared_ptr<arrow::RecordBatch> rb = ExtractHits<SSeqLoc>(alignments, qloc, sloc, scope);

//     const auto &wrt_sts = arrow_wrapper->AddRB2Batch(rb);
//     if (wrt_sts.ok())
//     {
//       return wrt_sts.ValueOrDie();
//     }
//     else
//     {
//       cerr << "Warn : Invalid Alignment RB (Returning Empty) : " << wrt_sts.status().detail() << std::endl
//                   << wrt_sts.status().message() << std::endl;
//     }

//     return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie();
//   }
//   else
//   {
//     static_assert(std::is_same_v<T, T>, "Unsupported type, only ncbi::blast::TSeqLocVector & ncbi::blast::SSeqLoc are supported");
//   }
// }

// template <>
// std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, const CScope &scope)
// {

//   // Rcpp::checkUserInterrupt();
//   assert(!alignments.empty());

//   std::string strand;
//   arrow::StringBuilder strand_builder;

//   auto query_strand = qloc.seqloc->GetStrand();
//   auto subject_strand = sloc.seqloc->GetStrand();

//   switch (query_strand)
//   {
//   case eNa_strand_minus:
//     strand = strand + "-";
//     break;
//   case eNa_strand_plus:
//     strand = strand + "+";
//     break;
//   default:
//     strand = strand + "*";
//     break;
//   }

//   switch (subject_strand)
//   {
//   case eNa_strand_minus:
//     strand = strand + "/-";
//     break;
//   case eNa_strand_plus:
//     strand = strand + "/+";
//     break;
//   default:
//     strand = strand + "/*";
//     break;
//   }

//   std::string qseq, sseq, frame, qseq_id, sseq_id;
//   arrow::StringBuilder qseqid_builder,
//       sseqid_builder;
//   arrow::LargeStringBuilder qseq_builder, sseq_builder;
//   arrow::Int8Builder qlen_builder, slen_builder, num_alignments_builder;

//   qseq_id = qloc.seqloc->GetId()->GetSeqIdString(true);
//   sseq_id = sloc.seqloc->GetId()->GetSeqIdString(true);

//   qseq = "";
//   sseq = "";
//   switch (save_sequences)
//   {
//   case true:
//     qseq = GetSSeqLocSequence(qloc);
//     sseq = GetSSeqLocSequence(sloc);
//     break;
//   }

//   arrow::Int32Builder hsp_offset_builder;

//   int num_rows = 0;

//   arrow::Int8Builder length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
//   arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
//   arrow::StringBuilder frames_builder;

//   for (const auto &seq_align_set : alignments)
//   {

//     if (seq_align_set->IsEmpty())
//     {
//       break;
//     }
//     assert(seq_align_set->IsSet());
//     assert(seq_align_set->CanGet());
//     const auto &seq_aligns = seq_align_set->Get();
//     assert(!seq_aligns.empty());

//     if (seq_aligns.size() > 0) // FILL UP THE ARRAYS
//     {

//       try
//       {

//         for (const auto &it : seq_aligns)
//         {

//           assert(!it.IsNull());
//           if (!it.NotEmpty())
//           {
//             break;
//           }

//           assert(it->CanGetScore());
//           int score, n_splices, num_ident, aln_len, aln_len01, gaps, mismatches, positive, qstart, qend, sstart, send, negative_count;
//           double bits, evalue, blast_score, pident, pident_gap, qcovhsp, sum_evalue, product_coverage, overall_identity, high_quality_percent_coverage, exon_identity, consensus_splices, comp_adj_method, matches;
//           std::string frames;

//           it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices);

//           it->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices);
//           it->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method);

//           aln_len01 = it->AlignLengthRatio();

//           qstart = it->GetSeqStart(0);
//           qend = it->GetSeqStop(0);
//           sstart = it->GetSeqStart(1);
//           send = it->GetSeqStop(1);

//           frames = std::to_string(GetFrame(qstart, aln_len, query_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, subject_strand));

//           static_cast<void>(frames_builder.Append(frames));
//           static_cast<void>(qstart_builder.Append(qstart));
//           static_cast<void>(qend_builder.Append(qend));
//           static_cast<void>(sstart_builder.Append(sstart));
//           static_cast<void>(send_builder.Append(send));
//           static_cast<void>(pident_builder.Append(pident));
//           static_cast<void>(evalue_builder.Append(evalue));
//           static_cast<void>(length_builder.Append(aln_len));
//           static_cast<void>(aln_len01_builder.Append(aln_len01));
//           static_cast<void>(bitscore_builder.Append(bits));
//           static_cast<void>(score_builder.Append(score));
//           static_cast<void>(qcovhsp_builder.Append(qcovhsp));
//           static_cast<void>(blast_score_builder.Append(blast_score));
//           static_cast<void>(pident_gap_builder.Append(pident_gap));
//           static_cast<void>(gaps_builder.Append(gaps));
//           static_cast<void>(nident_builder.Append(num_ident));
//           static_cast<void>(mismatch_builder.Append(mismatches));
//           static_cast<void>(positive_builder.Append(positive));
//           static_cast<void>(n_splices_builder.Append(n_splices));
//           static_cast<void>(hsp_cnt_builder.Append(num_rows + 1));
//           static_cast<void>(sum_evalue_builder.Append(sum_evalue));
//           static_cast<void>(product_coverage_builder.Append(product_coverage));
//           static_cast<void>(overall_identity_builder.Append(overall_identity));
//           static_cast<void>(negative_count_builder.Append(negative_count));
//           static_cast<void>(matches_builder.Append(matches));
//           static_cast<void>(high_quality_percent_coverage_builder.Append(high_quality_percent_coverage));
//           static_cast<void>(exon_identity_builder.Append(exon_identity));
//           static_cast<void>(consensus_splices_builder.Append(consensus_splices));
//           static_cast<void>(comp_adj_method_builder.Append(comp_adj_method));

//           /// SEQ INFO
//           static_cast<void>(qseqid_builder.Append(qseq_id));
//           static_cast<void>(sseqid_builder.Append(sseq_id));
//           static_cast<void>(qseq_builder.Append(qseq));
//           static_cast<void>(sseq_builder.Append(sseq));
//           static_cast<void>(qlen_builder.Append(qseq.length()));
//           static_cast<void>(slen_builder.Append(sseq.length()));
//           static_cast<void>(num_alignments_builder.Append(seq_aligns.size()));

//           static_cast<void>(strand_builder.Append(strand));
//           static_cast<void>(hsp_offset_builder.Append(1));

//           num_rows++;
//         }
//       }
//       catch (const std::exception &e)
//       {
//         cout << e.what() << std::endl
//                     << std::flush;
//       }
//     }
//     else
//     {
//       return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
//     }
//   }

//   if (num_rows == 0)
//   {
//     return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
//   }

//   std::shared_ptr<arrow::Array> frames_array;
//   static_cast<void>(frames_builder.Finish(&frames_array));
//   std::shared_ptr<arrow::Array> pident_array;
//   static_cast<void>(pident_builder.Finish(&pident_array));
//   std::shared_ptr<arrow::Array> pident_gap_array;
//   static_cast<void>(pident_gap_builder.Finish(&pident_gap_array));
//   std::shared_ptr<arrow::Array> evalue_array;
//   static_cast<void>(evalue_builder.Finish(&evalue_array));
//   std::shared_ptr<arrow::Array> length_array;
//   static_cast<void>(length_builder.Finish(&length_array));
//   std::shared_ptr<arrow::Array> qstart_array;
//   static_cast<void>(qstart_builder.Finish(&qstart_array));
//   std::shared_ptr<arrow::Array> qend_array;
//   static_cast<void>(qend_builder.Finish(&qend_array));
//   std::shared_ptr<arrow::Array> sstart_array;
//   static_cast<void>(sstart_builder.Finish(&sstart_array));
//   std::shared_ptr<arrow::Array> send_array;
//   static_cast<void>(send_builder.Finish(&send_array));
//   std::shared_ptr<arrow::Array> aln_len01_array;
//   static_cast<void>(aln_len01_builder.Finish(&aln_len01_array));
//   std::shared_ptr<arrow::Array> bitscore_array;
//   static_cast<void>(bitscore_builder.Finish(&bitscore_array));
//   std::shared_ptr<arrow::Array> score_array;
//   static_cast<void>(score_builder.Finish(&score_array));
//   std::shared_ptr<arrow::Array> qcovhsp_array;
//   static_cast<void>(qcovhsp_builder.Finish(&qcovhsp_array));
//   std::shared_ptr<arrow::Array> blast_score_array;
//   static_cast<void>(blast_score_builder.Finish(&blast_score_array));
//   std::shared_ptr<arrow::Array> gaps_array;
//   static_cast<void>(gaps_builder.Finish(&gaps_array));
//   std::shared_ptr<arrow::Array> nident_array;
//   static_cast<void>(nident_builder.Finish(&nident_array));
//   std::shared_ptr<arrow::Array> mismatch_array;
//   static_cast<void>(mismatch_builder.Finish(&mismatch_array));
//   std::shared_ptr<arrow::Array> positive_array;
//   static_cast<void>(positive_builder.Finish(&positive_array));
//   std::shared_ptr<arrow::Array> n_splices_array;
//   static_cast<void>(n_splices_builder.Finish(&n_splices_array));
//   std::shared_ptr<arrow::Array> hsp_cnt_array;
//   static_cast<void>(hsp_cnt_builder.Finish(&hsp_cnt_array));
//   std::shared_ptr<arrow::Array> sum_evalue_array;
//   static_cast<void>(sum_evalue_builder.Finish(&sum_evalue_array));
//   std::shared_ptr<arrow::Array> product_coverage_array;
//   static_cast<void>(product_coverage_builder.Finish(&product_coverage_array));
//   std::shared_ptr<arrow::Array> overall_identity_array;
//   static_cast<void>(overall_identity_builder.Finish(&overall_identity_array));
//   std::shared_ptr<arrow::Array> negative_count_array;
//   static_cast<void>(negative_count_builder.Finish(&negative_count_array));
//   std::shared_ptr<arrow::Array> matches_array;
//   static_cast<void>(matches_builder.Finish(&matches_array));
//   std::shared_ptr<arrow::Array> high_quality_percent_coverage_array;
//   static_cast<void>(high_quality_percent_coverage_builder.Finish(&high_quality_percent_coverage_array));
//   std::shared_ptr<arrow::Array> exon_identity_array;
//   static_cast<void>(exon_identity_builder.Finish(&exon_identity_array));
//   std::shared_ptr<arrow::Array> consensus_splices_array;
//   static_cast<void>(consensus_splices_builder.Finish(&consensus_splices_array));
//   std::shared_ptr<arrow::Array> comp_adj_method_array;
//   static_cast<void>(comp_adj_method_builder.Finish(&comp_adj_method_array));

//   arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({pident_array,
//                                                                                                   pident_gap_array,
//                                                                                                   frames_array,
//                                                                                                   evalue_array,
//                                                                                                   length_array,
//                                                                                                   aln_len01_array,
//                                                                                                   qstart_array,
//                                                                                                   qend_array,
//                                                                                                   sstart_array,
//                                                                                                   send_array,
//                                                                                                   bitscore_array,
//                                                                                                   score_array,
//                                                                                                   qcovhsp_array,
//                                                                                                   blast_score_array,
//                                                                                                   gaps_array,
//                                                                                                   nident_array,
//                                                                                                   mismatch_array,
//                                                                                                   positive_array,
//                                                                                                   n_splices_array,
//                                                                                                   hsp_cnt_array,
//                                                                                                   sum_evalue_array,
//                                                                                                   product_coverage_array,
//                                                                                                   overall_identity_array,
//                                                                                                   negative_count_array,
//                                                                                                   matches_array,
//                                                                                                   high_quality_percent_coverage_array,
//                                                                                                   exon_identity_array,
//                                                                                                   consensus_splices_array,
//                                                                                                   comp_adj_method_array},
//                                                                                                  {"pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});

//   assert(aln_struct_array.ok());

//   std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();

//   std::shared_ptr<arrow::Array> qseqid_array;
//   static_cast<void>(qseqid_builder.Finish(&qseqid_array));

//   std::shared_ptr<arrow::Array> sseqid_array;
//   static_cast<void>(sseqid_builder.Finish(&sseqid_array));

//   std::shared_ptr<arrow::Array> qseq_array;
//   static_cast<void>(qseq_builder.Finish(&qseq_array));

//   std::shared_ptr<arrow::Array> sseq_array;
//   static_cast<void>(sseq_builder.Finish(&sseq_array));

//   std::shared_ptr<arrow::Array> qlen_array;
//   static_cast<void>(qlen_builder.Finish(&qlen_array));

//   std::shared_ptr<arrow::Array> slen_array;
//   static_cast<void>(slen_builder.Finish(&slen_array));

//   std::shared_ptr<arrow::Array> strand_array;
//   static_cast<void>(strand_builder.Finish(&strand_array));

//   std::shared_ptr<arrow::Array> num_alignment_array;
//   static_cast<void>(num_alignments_builder.Finish(&num_alignment_array));

//   // Create the seq_info struct array and populate with the arrays
//   std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
//   std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
//   std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int8()), arrow::field("slen", arrow::int8())});

//   arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
//                                                                                                 seqids_struct_array,
//                                                                                                 seqs_struct_array,
//                                                                                                 strand_array,
//                                                                                                 lengths_struct_array},
//                                                                                                {"num_alignments", "seqids", "seqs", "strands", "lengths"});

//   assert(seq_info_array.ok());

//   std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();

//   std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
//                                                                               num_rows,
//                                                                               {seq_info_array_, aln_struct_array_});
//   if (alignment_rb)
//   {
//     return alignment_rb;
//   }

//   return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // //ERROR RETURN, END
// }

// std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data)
// {
//   // Rcpp::checkUserInterrupt();
//   std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr;
//   std::shared_ptr<arrow::Int32Builder> rec_no_builder;
//   std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
//   rec_no_builder = std::make_shared<arrow::Int32Builder>();
//   fasta_seq_builder = std::make_shared<arrow::StringBuilder>();
//   fasta_h_builder = std::make_shared<arrow::StringBuilder>();
//   static_cast<void>(rec_no_builder->Append(fasta_data.rec_no));
//   static_cast<void>(fasta_h_builder->Append(fasta_data.header));
//   static_cast<void>(fasta_seq_builder->Append(fasta_data.seq));
//   static_cast<void>(fasta_seq_builder->Finish(&seqArr));
//   static_cast<void>(fasta_h_builder->Finish(&hArr));
//   static_cast<void>(rec_no_builder->Finish(&recnoArr));
//   return arrow::RecordBatch::Make(arrow_wrapper->GetFASTASchema(), 1, {recnoArr, hArr, seqArr});
// }

// std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc)
// {
//   const CSeq_id &id = *(seq_loc.seqloc->GetId());

//   // Get the Bioseq using the Seq-id.
//   CBioseq_Handle bioseq_handle = seq_loc.scope->GetBioseqHandle(id);

//   // Terminate the program if the GI cannot be resolved to a Bioseq.
//   if (!bioseq_handle)
//   {
//     ERR_POST(Fatal << "Bioseq not found");
//   }

//   // Get the sequence using CSeqVector.
//   // Use Iupac encoding: CSeq_data::e_Iupacna or CSeq_data::e_Iupacaa.
//   // const auto &length = bioseq_handle.GetBioseqLength();
//   const auto &seq_vect_begin = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).begin();
//   const auto &seq_vect_end = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).end();

//   std::string str(seq_vect_begin, seq_vect_end);

//   return NStr::PrintableString(str);
// }

// std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
// {
//   // assert(num_threads > 0);
//   /*   if (!arrow_wrapper || arrow_wrapper.get() == nullptr)
//     {
//       arrow_wrapper = std::make_shared<ArrowWrapper>();
//     } */
//   /* if (this->opts.Empty() || this->opts.IsNull())
//   {
//     // if (this->blast_options_list.size() > 0)
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_list);
//     // }
//     // else if (!this->blast_options_str.empty())
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_str);
//     // }
//     // else
//     // {
//     //   this->opts = SetQuickBLASTOptions(this->program, "");
//     // }
//     this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, this->blast_options));
//   } */

// #if defined(_OPENMP) || defined(WIN32)
//   int n_threads = num_threads > omp_get_num_threads() ? omp_get_num_threads() : num_threads;
// #else
//   int n_threads = 1;
// #endif

//   n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
//   arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile);
//   if (!outfile_sts.ok())
//   {
//     /* Rcpp::Rcerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
//                 << outfile_sts.message() << std::endl; */
//       cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
//           << outfile_sts.message() << std::endl;
//     return std::make_shared<arrow::RecordBatchVector>();
//   }

//   SetThreadCount(n_threads);

//   int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);

//   int s_seq_count = arrow_wrapper->CountCharacter(subjectFile, '>', n_threads);

//   const int totalIterations = q_seq_count * s_seq_count;
//   if (blast_sequence_limit > 0)
//   {
//     blast_sequence_limit = blast_sequence_limit > totalIterations ? totalIterations : blast_sequence_limit;
//     blast_sequence_limit = blast_sequence_limit > s_seq_count ? s_seq_count - 1 : blast_sequence_limit; // - 1;
//   }
//   else
//   {
//     blast_sequence_limit = s_seq_count - 1;
//   }
//   assert(totalIterations > 0);

//   //// int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
//   ////  batch_size = 32 * n_threads; // batch_size > 0 ? batch_size : 1024;
//   arrow_wrapper->SetBatchSize(batch_size);

//   // Progress progress_bar(totalIterations, show_progress);

//   std::shared_ptr<arrow::RecordBatchVector> final_ret = StreamFile<FastaSequenceData>(
//       queryFile, ">", n_threads, [this, n_threads, subjectFile, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_q) //&progress_bar
//       {
//                                                                                         CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
//                                                                                         const std::shared_ptr<SSeqLoc> query_loc(std::move(CreateSSeqLocFromType<FastaSequenceData>(*data_q, scope)));

//                                                                                         _ASSERT(query_loc->seqloc.NotEmpty());

//                                                                                         // Rcpp::checkUserInterrupt();

//                                                                                         std::shared_ptr<TSeqLocVector> subjects_loc_vec(new TSeqLocVector());

// #if defined(_OPENMP) || defined(WIN32)
//                                                                                           omp_lock_t query_locLock;
//                                                                                           omp_lock_t subjects_loc_vecLock;
//                                                                                           omp_init_lock(&query_locLock);
//                                                                                           omp_init_lock(&subjects_loc_vecLock);
// #endif

//                                                                                           std::shared_ptr<arrow::RecordBatchVector> ret_results = StreamFile<FastaSequenceData>(
//                                                                                               subjectFile, ">", n_threads, [this, query_loc,
// #if defined(_OPENMP) || defined(WIN32)
//                                                                                                                             &query_locLock, &subjects_loc_vecLock,
// #endif
//                                                                                                                             &scope, &subjects_loc_vec, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_s) // &progress_bar
//                                                                                               {
//                                                                                                 const std::unique_ptr<SSeqLoc> subject_loc(CreateSSeqLocFromType<FastaSequenceData>(*data_s, scope));
//                                                                                                 _ASSERT(subject_loc->seqloc.NotEmpty());

//                                                                                                 if (strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) != 0)
//                                                                                                 {
//                                                                                                   // Rcpp::checkUserInterrupt();

//                                                                                                   CBl2Seq *blaster;

//                                                                                                   try
//                                                                                                   {

//                                                                                                     switch (blast_sequence_limit)
//                                                                                                     {
//                                                                                                     case 0:
//                                                                                                     {
//                                                                                                       blaster = new CBl2Seq(*query_loc, *subject_loc, this->GetQuickBLASTOptions());
//                                                                                                       arrow::RecordBatchVector tmp_rbv = {ExtractHits<SSeqLoc>(blaster->Run(), *query_loc, *subject_loc, *scope)};
//                                                                                                       // progress_bar.increment();
//                                                                                                       if(return_values)
//                                                                                                       {
//                                                                                                         return std::make_shared<arrow::RecordBatchVector>(tmp_rbv);
//                                                                                                       }
//                                                                                                       else
//                                                                                                       {
//                                                                                                         tmp_rbv.clear();
//                                                                                                         return std::make_shared<arrow::RecordBatchVector>();
//                                                                                                       }
//                                                                                                     }
//                                                                                                     break;
//                                                                                                     default:
//                                                                                                     {
// #if defined(_OPENMP) || defined(WIN32)
//                                                                                                       omp_set_lock(&subjects_loc_vecLock);
// #endif
//                                                                                                       subjects_loc_vec->emplace_back(*subject_loc);
// #if defined(_OPENMP) || defined(WIN32)
//                                                                                                       omp_unset_lock(&subjects_loc_vecLock);
// #endif

//                                                                                                       // progress_bar.increment();
//                                                                                                       if (subjects_loc_vec->size() >= blast_sequence_limit)
//                                                                                                       {
//                                                                                                         TSeqLocVector subjects_buffer_vec;
// #if defined(_OPENMP) || defined(WIN32)
//                                                                                                         omp_set_lock(&subjects_loc_vecLock);

// #endif
//                                                                                                         subjects_buffer_vec.swap(*subjects_loc_vec);
//                                                                                                         subjects_loc_vec->clear();
// #if defined(_OPENMP) || defined(WIN32)
//                                                                                                         omp_unset_lock(&subjects_loc_vecLock);
// #endif
//                                                                                                         blaster = new CBl2Seq(*query_loc, subjects_buffer_vec, this->GetQuickBLASTOptions(), true);

//                                                                                                         AddHitCount(subjects_buffer_vec.size());
//                                                                                                         std::shared_ptr<arrow::RecordBatchVector> tmp_rbv = ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
//                                                                                                         //return ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
//                                                                                                         if(return_values)
//                                                                                                         {
//                                                                                                           return tmp_rbv;
//                                                                                                         }
//                                                                                                         else
//                                                                                                         {
//                                                                                                           tmp_rbv->clear();
//                                                                                                           return std::make_shared<arrow::RecordBatchVector>();
//                                                                                                         }
//                                                                                                       }
//                                                                                                     }
//                                                                                                     break;
//                                                                                                     }
//                                                                                                   }
//                                                                                                   catch (const CException &e)
//                                                                                                   {
//                                                                                                     // Handle exception ...
//                                                                                                     cout << e.GetFunction() << std::endl;
//                                                                                                     cout << e.GetErrCodeString() << std::endl;
//                                                                                                     cout << e.GetErrCode() << std::endl;
//                                                                                                     cout << e.GetModule() << std::endl;
//                                                                                                     cout << e.GetPredecessor() << std::endl;
//                                                                                                     cout << e.GetFile() << std::endl;
//                                                                                                     cout << e.GetLine() << std::endl;
//                                                                                                     cout << e.GetMsg() << std::endl;
//                                                                                                     cout << e.GetStackTrace() << std::endl;
//                                                                                                     cout << e.GetStackTraceLevel() << std::endl;
//                                                                                                     cout << e.GetClass() << std::endl
//                                                                                                                 << std::flush;
//                                                                                                   }
//                                                                                                   catch (const std::exception &e)
//                                                                                                   {
//                                                                                                     cout << e.what() << std::endl
//                                                                                                                 << std::flush;
//                                                                                                   }
//                                                                                                 }

//                                                                                                 return std::make_shared<arrow::RecordBatchVector>(); // EMPTY ERROR Return
//                                                                                               },
//                                                                                               return_values);

//                                                                                           if (subjects_loc_vec->size() > 0)
//                                                                                           {
//                                                                                             CBl2Seq blaster(*query_loc, *subjects_loc_vec, this->GetQuickBLASTOptions(), true);
//                                                                                             AddHitCount(subjects_loc_vec->size());
//                                                                                             std::shared_ptr<arrow::RecordBatchVector> ret_vec = ExtractHits(blaster.Run(), *query_loc, *subjects_loc_vec, *scope);

//                                                                                             if (return_values)
//                                                                                             {
//                                                                                               ret_results->insert(ret_results->end(), ret_vec->begin(), ret_vec->end());
//                                                                                             }
//                                                                                             else
//                                                                                             {
//                                                                                               ret_vec->clear();
//                                                                                             }
//                                                                                         }

// #if defined(_OPENMP) || defined(WIN32)
//                                                                                             omp_destroy_lock(&query_locLock);
//                                                                                             omp_destroy_lock(&subjects_loc_vecLock);
// #endif
//                                                             if(return_values){
//                                                                                         return ret_results;
//                                                                                         }else
//                                                                                         {
//                                                                                           ret_results->clear();
//                                                                                         return std::make_shared<arrow::RecordBatchVector>();
//                                                                                         } },
//       return_values);

// #if defined(_OPENMP) || defined(WIN32)
// #pragma omp barrier
// #endif

//   arrow_wrapper->FinishOutputStream();
//   if (return_values)
//   {
//     return final_ret;
//   }
//   else
//   {
//     final_ret->clear();
//     return std::make_shared<arrow::RecordBatchVector>();
//   }
// }

// std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject)
// {
//   // Rcpp::checkUserInterrupt();

//   CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));

//   std::unique_ptr<SSeqLoc>
//       query_seqloc(CreateSSeqLocFromType<std::string>(query, scope));
//   std::unique_ptr<SSeqLoc> subject_seqloc(CreateSSeqLocFromType<std::string>(subject, scope));

//   CBl2Seq blaster(*query_seqloc, *subject_seqloc, GetQuickBLASTOptions());

//   return ExtractHits<SSeqLoc>(blaster.Run(), *query_seqloc, *subject_seqloc, *scope);
// }

// inline std::map<int, std::shared_ptr<QuickBLAST>> QuickBLAST::GetInstanceList()
// {
//     return obj_list;
// }

// //' @name BLAST C++ Call
// //' @title BLAST C++ Call
// //'
// //' @description BLAST 2 Files/Seqs. This is for the QuickBLAST C++ object exposed in R
// //'
// //' @param query (string) Query FASTA File/Seq
// //' @param subject (string) Subject FASTA File/Seq
// //' @param outputFile (string) Output Filename (Arrow Feather/IPC Format)  - Not used for Sequence BLAST
// //' @param input_type - (QuickBLAST::EInputType) 0 - eFile, 1 - eSequenceString
// //' @param blast_sequence_limit (int) Batch Size to BLAST at a time  - Not used for Sequence BLAST
// //' @param show_progress (bool) TRUE - Show progress, Set FALSE for multiple instances - Not used for Sequence BLAST
// //' @return Nested List of BLAST Hits
// auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
// {

//   assert(std::filesystem::exists(query));
//   assert(std::filesystem::exists(subject));

//   switch (input_type)
//   {
//   case QuickBLAST::EInputType::eFile:
//   {
// #if defined(_OPENMP) || defined(WIN32)
//     int n_threads = omp_get_num_threads();
// #else
//     int n_threads = 1;
// #endif
//     return BLAST_files(query, subject, outputFile, blast_sequence_limit, n_threads, show_progress);
//     //return Hits2RList(*ret_val);
//   }
//   break;
//   case QuickBLAST::EInputType::eSequenceString:
//   {
//       arrow::RecordBatchVector ret_val;
//       ret_val.emplace_back(BLAST_seqs(query, subject));
//       return std::make_shared<arrow::RecordBatchVector>(ret_val);
//     //return Hits2RList(ret_val);
//   }
//   break;
//   default:
//   {
//     // Rcpp::Rcerr << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
//       cout << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
//       //return false; //Rcpp::wrap(false);
//       return std::make_shared<arrow::RecordBatchVector>();
//   }
//   break;
//   }
// }

// // #ifdef linux
// //     //SEXP QuickBLAST::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
// //     SEXP Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
// //     {

// //       // Dispatch based on the data type of the array
// //       if (type->id() == arrow::Type::STRUCT)
// //       {

// //         auto struct_array = std::static_pointer_cast<arrow::StructArray>(array);
// //         int num_fields = struct_array->num_fields();

// //         // Create an Rcpp list to hold the data frames representing each field of the struct
// //         Rcpp::List struct_list(num_fields);
// //         Rcpp::CharacterVector names(num_fields);

// //         for (int i = 0; i < num_fields; i++)
// //         {

// //           auto field_array = struct_array->field(i);
// //           auto field_type = type->field(i)->type();
// //           auto field_name = type->field(i)->name();
// //           names[i] = field_name;
// //           struct_list[i] = Hits2RList_internal(field_array, field_type, field_name);
// //         }

// //         struct_list.names() = names;

// //         return struct_list;
// //       }
// //       else if (type->id() == arrow::Type::LIST)
// //       {

// //         auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
// //         auto value_type = type->field(0)->type();

// //         // Convert the list array to an Rcpp list
// //         Rcpp::List list_values(list_array->length());
// //         Rcpp::CharacterVector names(list_array->length());

// //         for (int i = 0; i < list_array->length(); i++)
// //         {

// //           auto sublist_array = list_array->values()->Slice(list_array->value_offset(i), list_array->value_length(i));

// //           names[i] = field_name + "[" + std::to_string(i) + "]";
// //           auto sublist_name = field_name + "[" + std::to_string(i) + "]";
// //           list_values[i] = Hits2RList_internal(sublist_array, value_type, sublist_name);
// //         }

// //         list_values.names() = names;

// //         return list_values;
// //       }
// //       else if (type->id() == arrow::Type::STRING || type->id() == arrow::Type::LARGE_STRING)
// //       {

// //         auto string_array = std::static_pointer_cast<arrow::StringArray>(array);

// //         Rcpp::StringVector strings(string_array->length());

// //         for (int i = 0; i < string_array->length(); ++i)
// //         {

// //           if (string_array->IsValid(i))
// //           {
// //             strings[i] = Rcpp::String(string_array->GetString(i));
// //           }
// //           // else
// //           // {
// //           //   strings[i] = NA_STRING;
// //           // }
// //         }

// //         return strings;
// //       }
// //       else if (type->id() == arrow::Type::INT8)
// //       {

// //         auto int_array = std::static_pointer_cast<arrow::Int8Array>(array);

// //         Rcpp::IntegerVector ints(int_array->length());

// //         for (int i = 0; i < int_array->length(); ++i)
// //         {

// //           if (int_array->IsValid(i))
// //           {
// //             ints[i] = int_array->Value(i);
// //           }
// //           // else
// //           // {
// //           //   ints[i] = NA_INTEGER;
// //           // }
// //         }

// //         return ints;
// //       }
// //       else if (type->id() == arrow::Type::DOUBLE)
// //       { // Use arrow::Type::DOUBLE instead of FLOAT64

// //         auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
// //         Rcpp::NumericVector doubles(double_array->length());

// //         for (int i = 0; i < double_array->length(); ++i)
// //         {

// //           if (double_array->IsValid(i))
// //           {
// //             doubles[field_name] = double_array->Value(i);
// //           }
// //           // else
// //           // {
// //           //   doubles[i] = NA_REAL;
// //           // }
// //         }

// //         return doubles;
// //       }
// //       else
// //       {
// //         // For other data types that don't have a direct conversion, return R_NilValue (NA)
// //         return R_NilValue;
// //       }
// //     }

// //     //SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
// //     SEXP Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
// //     {
// //       // Assuming the schema of the RecordBatch is accessible here

// //       auto rb_schema = rb->schema();

// //       // Convert each column of the RecordBatch to R objects and store in a list
// //       Rcpp::List result_list(rb_schema->num_fields());

// //       for (int i = 0; i < rb_schema->num_fields(); ++i)
// //       {

// //         auto array = rb->column(i);
// //         auto field_type = rb_schema->field(i)->type();
// //         auto field_name = rb_schema->field(i)->name();
// //         result_list[i] = Hits2RList_internal(array, field_type, field_name);
// //       }

// //       return result_list;
// //     }

// //     //SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
// //     SEXP Hits2RList(const arrow::RecordBatchVector &rb_vector)
// //     {
// //       Rcpp::List result_list(rb_vector.size());

// //       // Traverse the vector of RecordBatches and convert each RecordBatch
// //       for (size_t i = 0; i < rb_vector.size(); ++i)
// //       {
// //         std::shared_ptr<arrow::RecordBatch> rb = rb_vector[i];
// //         result_list[i] = Hits2RList(rb);
// //       }

// //       return result_list;
// //     }
// // #endif

// std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input)
// {
//   std::vector<std::pair<std::string, std::string>> keyValuePairs;
//   std::istringstream iss(input);
//   std::string token;

//   while (iss >> token)
//   {
//     if (token[0] == '-')
//     {
//       // Extract key-value pair
//       std::string key = token.substr(1);
//       std::string value;

//       if (iss >> value)
//       {
//         keyValuePairs.emplace_back(key, value);
//       }
//       else
//       {
//         // Handle error: Missing value for key
//         cerr << "Error: Missing value for key '" << key << "'." << std::endl;
//         break;
//       }
//     }
//     else
//     {
//       // Handle error: Invalid token (not starting with '-')
//       cerr << "Error: Invalid token '" << token << "'." << std::endl;
//       break;
//     }
//   }

//   return keyValuePairs;
// }

// /* void QuickBLAST::PrintProgressBar(int current, int total, int barWidth)
// {
//   Rcpp::Rcout << "\033[2J\033[1;1H";
//   float progress = static_cast<float>(current) / total;
//   int filledWidth = static_cast<int>(progress * barWidth);

//   Rcpp::Rcout << "[" << std::setw(barWidth) << std::left << std::string(filledWidth, '=') << std::right << "] ";
//   Rcpp::Rcout << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
//   Rcpp::Rcout << "  (" << current << " / " << total << ")";
//   Rcpp::Rcout << std::flush;
// } */

// template <typename T1>
// std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
// {

//   if constexpr (!std::is_same_v<T1, std::string> || !std::is_same_v<T1, FastaSequenceData>)
//   {
//     static_assert(std::is_same_v<T1, T1>, "Unsupported type, only std::string & FastaSequenceData are supported");
//   }

//   return arrow_wrapper->SplitFilesIntoEntries<T1>(filename, delim, num_threads, Entry_callback, return_values);
// }

// #endif

// #if defined(linux) || defined(MINGW32)

#include <RcppCommon.h>
#include <Rcpp.h>
#include <chrono>
#include <iostream>
#include <string_view>
#include <map>
#include <tuple>
#include <future>
#include <iomanip>
#include <sys/mman.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unistd.h>
#include <cassert>
#include <thread>
#include <filesystem>

#include <algo/blast/QuickBLAST/commons.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
// #include <ArrowWrapper.cpp>
// #include <algo/blast/QuickBLAST/quick_blast-functions.hpp>
// #include <algo/blast/QuickBLAST/arrow_wrapper-functions.hpp>

template <typename OptionsType>
ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const OptionsType &options)
{
  this->blast_options = options;
  this->program = program_name;
  if constexpr (std::is_same_v<OptionsType, std::string> || std::is_same_v<OptionsType, Rcpp::String>)
  {
    // this->blast_options_str = options;
    return SetQuickBLASTOptions<std::string>(program_name, options);
  }
  else if constexpr (std::is_same_v<OptionsType, Rcpp::List>)
  {
    // this->blast_options_list = options;
    return SetQuickBLASTOptions<Rcpp::List>(program_name, options);
  }
  else
  {
    static_assert(std::is_same_v<OptionsType, OptionsType>, "Unsupported type");
  }
}

template <>
ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options)
{
  assert(!program_name.empty());

  ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
  // Create a CBlastOptionsHandle object
  ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
  opts->SetDefaults();
  // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
  // Example: Extracting and setting the BLAST database

  if (options.empty())
  {
    cout << "Using " << program_name << " Defaults..." << std::endl;
    return opts;
  }

  std::vector<std::pair<std::string, std::string>> keyValuePairs = BLASTOptionsFromString(options);

  std::unordered_map<std::string, std::size_t> hashMap;

  hashMap["evalue"] = std::hash<std::string>{}("evalue");
  hashMap["pident"] = std::hash<std::string>{}("pident");
  hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
  hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
  hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
  hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
  hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
  hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
  hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
  hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
  hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
  hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
  hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
  hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
  hashMap["window_size"] = std::hash<std::string>{}("window_size");

  for (const auto &pair : keyValuePairs)
  {
    std::string key_str = pair.second;
    std::size_t key = std::hash<std::string>{}(pair.first);

    if (key == hashMap["evalue"])
    {
      double val = std::stod(key_str);
      opts->SetEvalueThreshold(val);
    }
    else if (key == hashMap["pident"])
    {
      double val = std::stod(key_str);
      opts->SetPercentIdentity(val);
    }
    else if (key == hashMap["gapped_mode"])
    {
      bool val = (key_str == "TRUE" || key_str == "True" || key_str == "true" || key_str == "1");
      opts->SetGappedMode(val);
    }
    else if (key == hashMap["filter_string"])
    {
      std::string val = key_str;
      opts->SetFilterString(val.c_str());
    }
    else if (key == hashMap["effective_search_space"])
    {
      int val = std::stoi(key_str);
      opts->SetEffectiveSearchSpace(val);
    }
    else if (key == hashMap["cutoff_score"])
    {
      int val = std::stoi(key_str);
      opts->SetCutoffScore(val);
    }
    else if (key == hashMap["gap_trigger"])
    {
      double val = std::stod(key_str);
      opts->SetGapTrigger(val);
    }
    else if (key == hashMap["gap_x_dropoff"])
    {
      double val = std::stod(key_str);
      opts->SetGapXDropoff(val);
    }
    else if (key == hashMap["gap_x_dropoff_final"])
    {
      double val = std::stod(key_str);
      opts->SetGapXDropoffFinal(val);
    }
    else if (key == hashMap["hit_list_size"])
    {
      int val = std::stoi(key_str);
      opts->SetHitlistSize(val);
    }
    else if (key == hashMap["low_score_percentage"])
    {
      double val = std::stod(key_str);
      opts->SetLowScorePerc(val);
    }
    else if (key == hashMap["max_hsp_per_subject"])
    {
      int val = std::stoi(key_str);
      opts->SetMaxHspsPerSubject(val);
    }
    else if (key == hashMap["max_hsp_per_sequence"])
    {
      int val = std::stoi(key_str);
      opts->SetMaxNumHspPerSequence(val);
    }
    else if (key == hashMap["qcovhsp_perc"])
    {
      double val = std::stod(key_str);
      opts->SetQueryCovHspPerc(val);
    }
    else if (key == hashMap["window_size"])
    {
      int val = std::stoi(key_str);
      opts->SetWindowSize(val);
    }
  }

  opts->Validate();

  return opts;
}

template <>
ncbi::blast::CBlastOptionsHandle *QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const Rcpp::List &options)
{
  assert(!program_name.empty());
  Rcpp::List options_(options);
  for (int i = 0; i < options_.size(); ++i)
  {
    // Check if the element is empty
    if (Rf_isNull(options_[i]))
    {
      // Remove the empty element
      options_.erase(i);
      --i; // Decrement the index since the list size has changed
    }
  }
  ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
  // Create a CBlastOptionsHandle object
  ncbi::blast::CBlastOptionsHandle *opts = ncbi::blast::CBlastOptionsFactory::Create(program);
  opts->SetDefaults();
  // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
  // Example: Extracting and setting the BLAST database

  if (options_.size() == 0 || options_.isNULL())
  {
    Rcpp::Rcout << "Using " << program_name << " Defaults..." << std::endl;
    return opts;
  }

  Rcpp::CharacterVector keys = options_.names();

  std::unordered_map<std::string, std::size_t> hashMap;

  hashMap["evalue"] = std::hash<std::string>{}("evalue");
  hashMap["pident"] = std::hash<std::string>{}("pident");
  hashMap["gapped_mode"] = std::hash<std::string>{}("gapped_mode");
  hashMap["filter_string"] = std::hash<std::string>{}("filter_string");
  hashMap["effective_search_space"] = std::hash<std::string>{}("effective_search_space");
  hashMap["cutoff_score"] = std::hash<std::string>{}("cutoff_score");
  hashMap["gap_trigger"] = std::hash<std::string>{}("gap_trigger");
  hashMap["gap_x_dropoff"] = std::hash<std::string>{}("gap_x_dropoff");
  hashMap["gap_x_dropoff_final"] = std::hash<std::string>{}("gap_x_dropoff_final");
  hashMap["hit_list_size"] = std::hash<std::string>{}("hit_list_size");
  hashMap["low_score_percentage"] = std::hash<std::string>{}("low_score_percentage");
  hashMap["max_hsp_per_subject"] = std::hash<std::string>{}("max_hsp_per_subject");
  hashMap["max_hsp_per_sequence"] = std::hash<std::string>{}("max_hsp_per_sequence");
  hashMap["qcovhsp_perc"] = std::hash<std::string>{}("qcovhsp_perc");
  hashMap["window_size"] = std::hash<std::string>{}("window_size");

  for (int i = 0; i < keys.size(); ++i)
  {
    std::string key_str = Rcpp::as<std::string>(keys[i]);
    std::size_t key = std::hash<std::string>{}(key_str);

    if (key == hashMap["evalue"])
    {
      double val = Rcpp::as<double>(options_["evalue"]);
      opts->SetEvalueThreshold(val);
    }
    else if (key == hashMap["pident"])
    {
      double val = Rcpp::as<double>(options_["pident"]);
      opts->SetPercentIdentity(val);
    }
    else if (key == hashMap["gapped_mode"])
    {
      bool val = Rcpp::as<bool>(options_["gapped_mode"]);
      opts->SetGappedMode(val);
    }
    else if (key == hashMap["filter_string"])
    {
      std::string val = Rcpp::as<std::string>(options_["filter_string"]);
      opts->SetFilterString(val.c_str());
    }
    else if (key == hashMap["effective_search_space"])
    {
      int val = Rcpp::as<int>(options_["effective_search_space"]);
      opts->SetEffectiveSearchSpace(val);
    }
    else if (key == hashMap["cutoff_score"])
    {
      int val = Rcpp::as<int>(options_["cutoff_score"]);
      opts->SetCutoffScore(val);
    }
    else if (key == hashMap["gap_trigger"])
    {
      double val = Rcpp::as<double>(options_["gap_trigger"]);
      opts->SetGapTrigger(val);
    }
    else if (key == hashMap["gap_x_dropoff"])
    {
      double val = Rcpp::as<double>(options_["gap_x_dropoff"]);
      opts->SetGapXDropoff(val);
    }
    else if (key == hashMap["gap_x_dropoff_final"])
    {
      double val = Rcpp::as<double>(options_["gap_x_dropoff_final"]);
      opts->SetGapXDropoffFinal(val);
    }
    else if (key == hashMap["hit_list_size"])
    {
      int val = Rcpp::as<int>(options_["hit_list_size"]);
      opts->SetHitlistSize(val);
    }
    else if (key == hashMap["low_score_percentage"])
    {
      double val = Rcpp::as<double>(options_["low_score_percentage"]);
      opts->SetLowScorePerc(val);
    }
    else if (key == hashMap["max_hsp_per_subject"])
    {
      int val = Rcpp::as<int>(options_["max_hsp_per_subject"]);
      opts->SetMaxHspsPerSubject(val);
    }
    else if (key == hashMap["max_hsp_per_sequence"])
    {
      int val = Rcpp::as<int>(options_["max_hsp_per_sequence"]);
      opts->SetMaxNumHspPerSequence(val);
    }
    else if (key == hashMap["qcovhsp_perc"])
    {
      double val = Rcpp::as<double>(options_["qcovhsp_perc"]);
      opts->SetQueryCovHspPerc(val);
    }
    else if (key == hashMap["window_size"])
    {
      int val = Rcpp::as<int>(options_["window_size"]);
      opts->SetWindowSize(val);
    }
  }

  opts->Validate();

  return opts;
}

QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, Rcpp::List options, bool save_sequences)
{
#ifdef _OPENMP
  this->num_threads = omp_get_num_threads();
#else
  this->num_threads = 1;
#endif
  arrow_wrapper = std::make_shared<ArrowWrapper>();
  this->save_sequences = save_sequences;
  this->program = program;
  this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<Rcpp::List>(program, options));
  // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
  this->strand = strand;
  this->seq_type = seq_type;
  ok_promise.set_value(arrow::Status::OK());
#ifdef _OPENMP
  omp_init_lock(&hit_countLock);
#endif
}

QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences)
{

#if defined(_OPENMP) || defined(WIN32)
  this->num_threads = omp_get_num_threads();
#else
  this->num_threads = 1;
#endif
  arrow_wrapper = std::make_shared<ArrowWrapper>();
  this->save_sequences = save_sequences;
  this->program = program;
  this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions<std::string>(program, options));
  // this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, options));
  this->strand = strand;
  this->seq_type = seq_type;
  ok_promise.set_value(arrow::Status::OK());
#if defined(_OPENMP) || defined(WIN32)
  omp_init_lock(&hit_countLock);
#endif
}

QuickBLAST::~QuickBLAST()
{
#if defined(_OPENMP) || defined(WIN32)
  omp_destroy_lock(&hit_countLock);
#endif

  // DO NOT DELETE NCBI C++ OBJECTs or PTRs or face Corruption
  //  delete self;
  //  opts->ReleaseReference();
  // delete opts;

  // delete arrow_wrapper;

  cout << "~QuickBLAST " << std::endl;
}

// Function to process a single FASTA block
void QuickBLAST::PrintFastaBlock(FastaSequenceData *data, std::shared_ptr<std::ostringstream> outputStream)
{
  if (outputStream != nullptr)
  {

    // Print FastaSequenceData object
    (*outputStream) << "No: " << data->rec_no << std::endl;
    (*outputStream) << "Header: " << data->header << std::endl;
    (*outputStream) << "Sequence: " << data->seq << std::endl;
    (*outputStream) << std::endl;
    outputStream->flush();
  }
}

template <>
SSeqLoc *QuickBLAST::CreateSSeqLocFromType(FastaSequenceData fasta_data, CRef<ncbi::CScope> parent_scope)
{
  int rec_no = fasta_data.rec_no;
  std::string fastaID(fasta_data.header.data());
  std::string fastaSequence(fasta_data.seq.data());

  const TSeqPos seqlen = fastaSequence.length();

  _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
  ncbi::CRef<ncbi::objects::CSeq_interval> interval(new ncbi::objects::CSeq_interval());
  interval->SetFrom(0);
  interval->SetTo(seqlen - 1);

  CRef<CSeq_id> id(new CSeq_id(fastaID, (ncbi::objects::CSeq_id::fParse_RawText | ncbi::objects::CSeq_id::fParse_PartialOK | ncbi::objects::CSeq_id::fParse_ValidLocal)));
  id->Select(CSeq_id_Base::E_Choice::e_Local);
  id->SetLocal().SetId(rec_no);
  id->SetLocal().SetStr(fastaID);

  CRef<CSeq_loc>
      cseq_loc_obj(new CSeq_loc());
  cseq_loc_obj->Select(CSeq_loc_Base::E_Choice::e_Whole);
  cseq_loc_obj->SetWhole()
      .SetLocal()
      .SetStr(fastaID);
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
  seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
  seq_inst->SetTopology(CSeq_inst::eTopology_linear);
  seq_inst->SetStrand(CSeq_inst_Base::TStrand::eStrand_ss);
  seq_inst->SetRepr(CSeq_inst_Base::ERepr::eRepr_raw);
  seq_inst->SetMol(CSeq_inst_Base::EMol::eMol_dna);
  seq_inst->SetTopology(CSeq_inst_Base::ETopology::eTopology_linear);
  seq_inst->SetStrand(CSeq_inst_Base::EStrand::eStrand_ss);
  seq_inst->SetLength(seqlen);

  CRef<ncbi::objects::CBioseq> bioseq(new CBioseq(*cseq_loc_obj, fastaID));
  bioseq->SetInst(*seq_inst);

  CRef<CSeq_entry>
      ret_entry(new CSeq_entry());
  ret_entry->SetSeq(*bioseq);

  parent_scope->AddTopLevelSeqEntry(*ret_entry);

  return new SSeqLoc(cseq_loc_obj.GetObject(), parent_scope.GetObject());
}

int QuickBLAST::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
{
  int frame = 0;
  if (strand == eNa_strand_plus)
  {
    frame = (start % 3) + 1;
  }
  else if (strand == eNa_strand_minus)
  {
    frame = -(((int)length - start - 1) % 3 + 1);
  }
  return frame;
}

template <>
std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, const CScope &scope)
{

  // Rcpp::checkUserInterrupt();
  assert(!alignments.empty());

  std::string strand;
  arrow::StringBuilder strand_builder;

  auto query_strand = qloc.seqloc->GetStrand();
  auto subject_strand = sloc.seqloc->GetStrand();

  switch (query_strand)
  {
  case eNa_strand_minus:
    strand = strand + "-";
    break;
  case eNa_strand_plus:
    strand = strand + "+";
    break;
  default:
    strand = strand + "*";
    break;
  }

  switch (subject_strand)
  {
  case eNa_strand_minus:
    strand = strand + "/-";
    break;
  case eNa_strand_plus:
    strand = strand + "/+";
    break;
  default:
    strand = strand + "/*";
    break;
  }

  std::string qseq, sseq, frame, qseq_id, sseq_id;
  arrow::StringBuilder qseqid_builder,
      sseqid_builder;
  arrow::LargeStringBuilder qseq_builder, sseq_builder;
  arrow::Int8Builder qlen_builder, slen_builder, num_alignments_builder;

  qseq_id = qloc.seqloc->GetId()->GetSeqIdString(true);
  sseq_id = sloc.seqloc->GetId()->GetSeqIdString(true);

  qseq = "";
  sseq = "";
  switch (save_sequences)
  {
  case true:
    qseq = GetSSeqLocSequence(qloc);
    sseq = GetSSeqLocSequence(sloc);
    break;
  }

  arrow::Int32Builder hsp_offset_builder;

  int num_rows = 0;

  arrow::Int8Builder length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
  arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
  arrow::StringBuilder frames_builder;

  for (const auto &seq_align_set : alignments)
  {

    if (seq_align_set->IsEmpty())
    {
      break;
    }
    assert(seq_align_set->IsSet());
    assert(seq_align_set->CanGet());
    const auto &seq_aligns = seq_align_set->Get();
    assert(!seq_aligns.empty());

    if (seq_aligns.size() > 0) // FILL UP THE ARRAYS
    {

      try
      {

        for (const auto &it : seq_aligns)
        {

          assert(!it.IsNull());
          if (!it.NotEmpty())
          {
            break;
          }

          assert(it->CanGetScore());
          int score, n_splices, num_ident, aln_len, aln_len01, gaps, mismatches, positive, qstart, qend, sstart, send, negative_count;
          double bits, evalue, blast_score, pident, pident_gap, qcovhsp, sum_evalue, product_coverage, overall_identity, high_quality_percent_coverage, exon_identity, consensus_splices, comp_adj_method, matches;
          std::string frames;

          it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentCoverage, qcovhsp);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_Score, score);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_PositiveCount, positive);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_Splices, n_splices);

          it->GetNamedScore(CSeq_align::EScoreType::eScore_SumEValue, sum_evalue);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_ProductCoverage, product_coverage);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_OverallIdentity, overall_identity);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_NegativeCount, negative_count);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_Matches, matches);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_HighQualityPercentCoverage, high_quality_percent_coverage);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_ExonIdentity, exon_identity);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_ConsensusSplices, consensus_splices);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_CompAdjMethod, comp_adj_method);

          aln_len01 = it->AlignLengthRatio();

          qstart = it->GetSeqStart(0);
          qend = it->GetSeqStop(0);
          sstart = it->GetSeqStart(1);
          send = it->GetSeqStop(1);

          frames = std::to_string(GetFrame(qstart, aln_len, query_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, subject_strand));

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
          static_cast<void>(qseqid_builder.Append(qseq_id));
          static_cast<void>(sseqid_builder.Append(sseq_id));
          static_cast<void>(qseq_builder.Append(qseq));
          static_cast<void>(sseq_builder.Append(sseq));
          static_cast<void>(qlen_builder.Append(qseq.length()));
          static_cast<void>(slen_builder.Append(sseq.length()));
          static_cast<void>(num_alignments_builder.Append(seq_aligns.size()));

          static_cast<void>(strand_builder.Append(strand));
          static_cast<void>(hsp_offset_builder.Append(1));

          num_rows++;
        }
      }
      catch (const std::exception &e)
      {
        cout << e.what() << std::endl
             << std::flush;
      }
    }
    else
    {
      return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
    }
  }

  if (num_rows == 0)
  {
    return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // CORRECT RETURN, NO ALIGNMENTS
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

  std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                              num_rows,
                                                                              {seq_info_array_, aln_struct_array_});
  if (alignment_rb)
  {
    return alignment_rb;
  }

  return arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema()).ValueOrDie(); // //ERROR RETURN, END
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::ExtractFASTA(const FastaSequenceData &fasta_data)
{
  // Rcpp::checkUserInterrupt();
  std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr;
  std::shared_ptr<arrow::Int32Builder> rec_no_builder;
  std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
  rec_no_builder = std::make_shared<arrow::Int32Builder>();
  fasta_seq_builder = std::make_shared<arrow::StringBuilder>();
  fasta_h_builder = std::make_shared<arrow::StringBuilder>();
  static_cast<void>(rec_no_builder->Append(fasta_data.rec_no));
  static_cast<void>(fasta_h_builder->Append(fasta_data.header));
  static_cast<void>(fasta_seq_builder->Append(fasta_data.seq));
  static_cast<void>(fasta_seq_builder->Finish(&seqArr));
  static_cast<void>(fasta_h_builder->Finish(&hArr));
  static_cast<void>(rec_no_builder->Finish(&recnoArr));
  return arrow::RecordBatch::Make(arrow_wrapper->GetFASTASchema(), 1, {recnoArr, hArr, seqArr});
}

std::string QuickBLAST::GetSSeqLocSequence(const SSeqLoc &seq_loc)
{
  const CSeq_id &id = *(seq_loc.seqloc->GetId());

  // Get the Bioseq using the Seq-id.
  CBioseq_Handle bioseq_handle = seq_loc.scope->GetBioseqHandle(id);

  // Terminate the program if the GI cannot be resolved to a Bioseq.
  if (!bioseq_handle)
  {
    ERR_POST(Fatal << "Bioseq not found");
  }

  // Get the sequence using CSeqVector.
  // Use Iupac encoding: CSeq_data::e_Iupacna or CSeq_data::e_Iupacaa.
  // const auto &length = bioseq_handle.GetBioseqLength();
  const auto &seq_vect_begin = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).begin();
  const auto &seq_vect_end = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac, ncbi::objects::eNa_strand_plus).end();

  std::string str(seq_vect_begin, seq_vect_end);

  return NStr::PrintableString(str);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, const std::string &outFile, unsigned int blast_sequence_limit, int num_threads, const bool show_progress, const bool return_values, int batch_size)
{
  // assert(num_threads > 0);
  /*   if (!arrow_wrapper || arrow_wrapper.get() == nullptr)
    {
      arrow_wrapper = std::make_shared<ArrowWrapper>();
    } */
  /* if (this->opts.Empty() || this->opts.IsNull())
  {
    // if (this->blast_options_list.size() > 0)
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_list);
    // }
    // else if (!this->blast_options_str.empty())
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, this->blast_options_str);
    // }
    // else
    // {
    //   this->opts = SetQuickBLASTOptions(this->program, "");
    // }
    this->opts = CRef<ncbi::blast::CBlastOptionsHandle>(SetQuickBLASTOptions(program, this->blast_options));
  } */

#if defined(_OPENMP) || defined(WIN32)
  int n_threads = num_threads > omp_get_num_threads() ? omp_get_num_threads() : num_threads;
#else
  int n_threads = 1;
#endif

  n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
  arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile);
  if (!outfile_sts.ok())
  {
    /* Rcpp::Rcerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
                << outfile_sts.message() << std::endl; */
    cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
         << outfile_sts.message() << std::endl;
    return std::make_shared<arrow::RecordBatchVector>();
  }

  SetThreadCount(n_threads);

  unsigned int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);

  unsigned int s_seq_count = arrow_wrapper->CountCharacter(subjectFile, '>', n_threads);

  const unsigned int totalIterations = q_seq_count * s_seq_count;
  if (blast_sequence_limit > 0)
  {
    blast_sequence_limit = blast_sequence_limit > totalIterations ? totalIterations : blast_sequence_limit;
    blast_sequence_limit = blast_sequence_limit > s_seq_count ? s_seq_count - 1 : blast_sequence_limit; // - 1;
  }
  else
  {
    blast_sequence_limit = s_seq_count - 1;
  }
  assert(totalIterations > 0);

  //// int batch_size = 96 * num_threads; // int(ceil(totalIterations / pow(2, n_threads))); // int(ceil(sqrt(totalIterations) * (n_threads * 2)) / 2);
  ////  batch_size = 32 * n_threads; // batch_size > 0 ? batch_size : 1024;
  arrow_wrapper->SetBatchSize(batch_size);

  // Progress progress_bar(totalIterations, show_progress);

  std::shared_ptr<arrow::RecordBatchVector> final_ret = StreamFile<FastaSequenceData>(
      queryFile, ">", n_threads, [this, n_threads, subjectFile, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_q) //&progress_bar
      {
            CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
            const std::shared_ptr<SSeqLoc> query_loc(std::move(CreateSSeqLocFromType<FastaSequenceData>(*data_q, scope)));

            _ASSERT(query_loc->seqloc.NotEmpty());

            // Rcpp::checkUserInterrupt();

            std::shared_ptr<TSeqLocVector> subjects_loc_vec(new TSeqLocVector());

#if defined(_OPENMP) || defined(WIN32)
            omp_lock_t query_locLock;
            omp_lock_t subjects_loc_vecLock;
            omp_init_lock(&query_locLock);
            omp_init_lock(&subjects_loc_vecLock);
#endif

            std::shared_ptr<arrow::RecordBatchVector> ret_results = StreamFile<FastaSequenceData>(
                subjectFile, ">", n_threads, [this, query_loc,
#if defined(_OPENMP) || defined(WIN32)
                &query_locLock, &subjects_loc_vecLock,
#endif
                & scope, &subjects_loc_vec, blast_sequence_limit, return_values](std::shared_ptr<FastaSequenceData> data_s) // &progress_bar
                {
                    const std::unique_ptr<SSeqLoc> subject_loc(CreateSSeqLocFromType<FastaSequenceData>(*data_s, scope));
                    _ASSERT(subject_loc->seqloc.NotEmpty());

                    if (strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) != 0)
                    {
                        // Rcpp::checkUserInterrupt();

                        CBl2Seq* blaster;

                        try
                        {

                            switch (blast_sequence_limit)
                            {
                            case 0:
                            {
                                blaster = new CBl2Seq(*query_loc, *subject_loc, this->GetQuickBLASTOptions());
                                arrow::RecordBatchVector tmp_rbv = { ExtractHits<SSeqLoc>(blaster->Run(), *query_loc, *subject_loc, *scope) };
                                // progress_bar.increment();
                                if (return_values)
                                {
                                    return std::make_shared<arrow::RecordBatchVector>(tmp_rbv);
                                }
                                else
                                {
                                    tmp_rbv.clear();
                                    return std::make_shared<arrow::RecordBatchVector>();
                                }
                            }
                            break;
                            default:
                            {
#if defined(_OPENMP) || defined(WIN32)
                                omp_set_lock(&subjects_loc_vecLock);
#endif
                                subjects_loc_vec->emplace_back(*subject_loc);
#if defined(_OPENMP) || defined(WIN32)
                                omp_unset_lock(&subjects_loc_vecLock);
#endif

                                // progress_bar.increment();
                                if (subjects_loc_vec->size() >= blast_sequence_limit)
                                {
                                    TSeqLocVector subjects_buffer_vec;
#if defined(_OPENMP) || defined(WIN32)
                                    omp_set_lock(&subjects_loc_vecLock);

#endif
                                    subjects_buffer_vec.swap(*subjects_loc_vec);
                                    subjects_loc_vec->clear();
#if defined(_OPENMP) || defined(WIN32)
                                    omp_unset_lock(&subjects_loc_vecLock);
#endif
                                    blaster = new CBl2Seq(*query_loc, subjects_buffer_vec, this->GetQuickBLASTOptions(), true);

                                    AddHitCount(subjects_buffer_vec.size());
                                    std::shared_ptr<arrow::RecordBatchVector> tmp_rbv = ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
                                    //return ExtractHits(blaster->Run(), *query_loc, subjects_buffer_vec, *scope);
                                    if (return_values)
                                    {
                                        return tmp_rbv;
                                    }
                                    else
                                    {
                                        tmp_rbv->clear();
                                        return std::make_shared<arrow::RecordBatchVector>();
                                    }
                                }
                            }
                            break;
                            }
                        }
                        catch (const CException& e)
                        {
                            // Handle exception ...
                            cout << e.GetFunction() << std::endl;
                            cout << e.GetErrCodeString() << std::endl;
                            cout << e.GetErrCode() << std::endl;
                            cout << e.GetModule() << std::endl;
                            cout << e.GetPredecessor() << std::endl;
                            cout << e.GetFile() << std::endl;
                            cout << e.GetLine() << std::endl;
                            cout << e.GetMsg() << std::endl;
                            cout << e.GetStackTrace() << std::endl;
                            cout << e.GetStackTraceLevel() << std::endl;
                            cout << e.GetClass() << std::endl
                                << std::flush;
                        }
                        catch (const std::exception& e)
                        {
                            cout << e.what() << std::endl
                                << std::flush;
                        }
                    }

                    return std::make_shared<arrow::RecordBatchVector>(); // EMPTY ERROR Return
                },
                return_values);

            if (subjects_loc_vec->size() > 0)
            {
                CBl2Seq blaster(*query_loc, *subjects_loc_vec, this->GetQuickBLASTOptions(), true);
                AddHitCount(subjects_loc_vec->size());
                std::shared_ptr<arrow::RecordBatchVector> ret_vec = ExtractHits(blaster.Run(), *query_loc, *subjects_loc_vec, *scope);

                if (return_values)
                {
                    ret_results->insert(ret_results->end(), ret_vec->begin(), ret_vec->end());
                }
                else
                {
                    ret_vec->clear();
                }
            }

#if defined(_OPENMP) || defined(WIN32)
            omp_destroy_lock(&query_locLock);
            omp_destroy_lock(&subjects_loc_vecLock);
#endif
            if (return_values) {
                return ret_results;
            }
            else
            {
                ret_results->clear();
                return std::make_shared<arrow::RecordBatchVector>();
            } },
      return_values);

#if defined(_OPENMP) || defined(WIN32)
#pragma omp barrier
#endif

  arrow_wrapper->FinishOutputStream();
  if (return_values)
  {
    return final_ret;
  }
  else
  {
    final_ret->clear();
    return std::make_shared<arrow::RecordBatchVector>();
  }
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject)
{
  // Rcpp::checkUserInterrupt();

  CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));

  std::unique_ptr<SSeqLoc>
      query_seqloc(CreateSSeqLocFromType<std::string>(query, scope));
  std::unique_ptr<SSeqLoc> subject_seqloc(CreateSSeqLocFromType<std::string>(subject, scope));

  CBl2Seq blaster(*query_seqloc, *subject_seqloc, GetQuickBLASTOptions());

  return ExtractHits<SSeqLoc>(blaster.Run(), *query_seqloc, *subject_seqloc, *scope);
}

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
auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, const std::string &outputFile, QuickBLAST::EInputType input_type, int blast_sequence_limit, const bool show_progress)
{

  assert(std::filesystem::exists(query));
  assert(std::filesystem::exists(subject));

  switch (input_type)
  {
  case QuickBLAST::EInputType::eFile:
  {
#if defined(_OPENMP) || defined(WIN32)
    int n_threads = omp_get_num_threads();
#else
    int n_threads = 1;
#endif
    return BLAST_files(query, subject, outputFile, blast_sequence_limit, n_threads, show_progress);
    // return Hits2RList(*ret_val);
  }
  break;
  case QuickBLAST::EInputType::eSequenceString:
  {
    arrow::RecordBatchVector ret_val;
    ret_val.emplace_back(BLAST_seqs(query, subject));
    return std::make_shared<arrow::RecordBatchVector>(ret_val);
    // return Hits2RList(ret_val);
  }
  break;
  default:
  {
    // Rcpp::Rcerr << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
    cout << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !";
    // return false; //Rcpp::wrap(false);
    return std::make_shared<arrow::RecordBatchVector>();
  }
  break;
  }
}

SEXP QuickBLAST::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
{

  // Dispatch based on the data type of the array
  if (type->id() == arrow::Type::STRUCT)
  {

    auto struct_array = std::static_pointer_cast<arrow::StructArray>(array);
    int num_fields = struct_array->num_fields();

    // Create an Rcpp list to hold the data frames representing each field of the struct
    Rcpp::List struct_list(num_fields);
    Rcpp::CharacterVector names(num_fields);

    for (int i = 0; i < num_fields; i++)
    {

      auto field_array = struct_array->field(i);
      auto field_type = type->field(i)->type();
      auto field_name = type->field(i)->name();
      names[i] = field_name;
      struct_list[i] = Hits2RList_internal(field_array, field_type, field_name);
    }

    struct_list.names() = names;

    return struct_list;
  }
  else if (type->id() == arrow::Type::LIST)
  {

    auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
    auto value_type = type->field(0)->type();

    // Convert the list array to an Rcpp list
    Rcpp::List list_values(list_array->length());
    Rcpp::CharacterVector names(list_array->length());

    for (int i = 0; i < list_array->length(); i++)
    {

      auto sublist_array = list_array->values()->Slice(list_array->value_offset(i), list_array->value_length(i));

      names[i] = field_name + "[" + std::to_string(i) + "]";
      auto sublist_name = field_name + "[" + std::to_string(i) + "]";
      list_values[i] = Hits2RList_internal(sublist_array, value_type, sublist_name);
    }

    list_values.names() = names;

    return list_values;
  }
  else if (type->id() == arrow::Type::STRING || type->id() == arrow::Type::LARGE_STRING)
  {

    auto string_array = std::static_pointer_cast<arrow::StringArray>(array);

    Rcpp::StringVector strings(string_array->length());

    for (int i = 0; i < string_array->length(); ++i)
    {

      if (string_array->IsValid(i))
      {
        strings[i] = Rcpp::String(string_array->GetString(i));
      }
      // else
      // {
      //   strings[i] = NA_STRING;
      // }
    }

    return strings;
  }
  else if (type->id() == arrow::Type::INT8)
  {

    auto int_array = std::static_pointer_cast<arrow::Int8Array>(array);

    Rcpp::IntegerVector ints(int_array->length());

    for (int i = 0; i < int_array->length(); ++i)
    {

      if (int_array->IsValid(i))
      {
        ints[i] = int_array->Value(i);
      }
      // else
      // {
      //   ints[i] = NA_INTEGER;
      // }
    }

    return ints;
  }
  else if (type->id() == arrow::Type::DOUBLE)
  { // Use arrow::Type::DOUBLE instead of FLOAT64

    auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
    Rcpp::NumericVector doubles(double_array->length());

    for (int i = 0; i < double_array->length(); ++i)
    {

      if (double_array->IsValid(i))
      {
        doubles[field_name] = double_array->Value(i);
      }
      // else
      // {
      //   doubles[i] = NA_REAL;
      // }
    }

    return doubles;
  }
  else
  {
    // For other data types that don't have a direct conversion, return R_NilValue (NA)
    return R_NilValue;
  }
}

SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
{
  // Assuming the schema of the RecordBatch is accessible here

  auto rb_schema = rb->schema();

  // Convert each column of the RecordBatch to R objects and store in a list
  Rcpp::List result_list(rb_schema->num_fields());

  for (int i = 0; i < rb_schema->num_fields(); ++i)
  {

    auto array = rb->column(i);
    auto field_type = rb_schema->field(i)->type();
    auto field_name = rb_schema->field(i)->name();
    result_list[i] = Hits2RList_internal(array, field_type, field_name);
  }

  return result_list;
}

SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
{
  Rcpp::List result_list(rb_vector.size());

  // Traverse the vector of RecordBatches and convert each RecordBatch
  for (size_t i = 0; i < rb_vector.size(); ++i)
  {
    std::shared_ptr<arrow::RecordBatch> rb = rb_vector[i];
    result_list[i] = Hits2RList(rb);
  }

  return result_list;
}

std::vector<std::pair<std::string, std::string>> QuickBLAST::BLASTOptionsFromString(const std::string &input)
{
  std::vector<std::pair<std::string, std::string>> keyValuePairs;
  std::istringstream iss(input);
  std::string token;

  while (iss >> token)
  {
    if (token[0] == '-')
    {
      // Extract key-value pair
      std::string key = token.substr(1);
      std::string value;

      if (iss >> value)
      {
        keyValuePairs.emplace_back(key, value);
      }
      else
      {
        // Handle error: Missing value for key
        cerr << "Error: Missing value for key '" << key << "'." << std::endl;
        break;
      }
    }
    else
    {
      // Handle error: Invalid token (not starting with '-')
      cerr << "Error: Invalid token '" << token << "'." << std::endl;
      break;
    }
  }

  return keyValuePairs;
}

/* void QuickBLAST::PrintProgressBar(int current, int total, int barWidth)
{
  Rcpp::Rcout << "\033[2J\033[1;1H";
  float progress = static_cast<float>(current) / total;
  int filledWidth = static_cast<int>(progress * barWidth);

  Rcpp::Rcout << "[" << std::setw(barWidth) << std::left << std::string(filledWidth, '=') << std::right << "] ";
  Rcpp::Rcout << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
  Rcpp::Rcout << "  (" << current << " / " << total << ")";
  Rcpp::Rcout << std::flush;
} */

template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback, bool return_values)
{

  if constexpr (!std::is_same_v<T1, std::string> || !std::is_same_v<T1, FastaSequenceData>)
  {
    static_assert(std::is_same_v<T1, T1>, "Unsupported type, only std::string & FastaSequenceData are supported");
  }

  return arrow_wrapper->SplitFilesIntoEntries<T1>(filename, delim, num_threads, Entry_callback, return_values);
}

// #endif // linux
