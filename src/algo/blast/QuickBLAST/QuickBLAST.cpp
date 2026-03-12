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
#include <algo/blast/QuickBLAST/QuickBLAST.hpp>
#include <algo/blast/QuickBLAST/ArrowWrapper.hpp>

std::shared_ptr<ArrowWrapper> QuickBLAST::Impl::GetArrowWrapper() { return arrow_wrapper; }

std::shared_ptr<arrow::Schema> QuickBLAST::Impl::GetSchema() { 
  if (!arrow_wrapper) {
    // initialize default options or throw helpful error
    Rcpp::stop("QuickBLAST::GetSchema(): arrow_wrapper is not initialised");
  }
  return arrow_wrapper->GetSchema(); 
};

QuickBLAST::ESeqType QuickBLAST::GetSeqType(){
  return pImpl->GetSeqType();
}
QuickBLAST::ESeqType QuickBLAST::Impl::GetSeqType(){
  return seq_type;
}
void QuickBLAST::Impl::SetThreadCount(unsigned int num_threads)
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_set_num_threads(num_threads);
#endif
  this->num_threads = num_threads;
  arrow_wrapper->SetThreadCount(this->num_threads);
}
unsigned int QuickBLAST::Impl::GetThreadCount()
{
  return num_threads;
}

int QuickBLAST::Impl::GetHitCount()
{
  return hit_count;
}
void QuickBLAST::Impl::AddHitCount(int val)
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_set_lock(&hit_countLock);
#endif
  hit_count += val;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_unset_lock(&hit_countLock);
#endif
}

unsigned int QuickBLAST::Impl::GetProcRecordCount() { return arrow_wrapper->GetProcRecordCount(); }
unsigned int QuickBLAST::GetProcRecordCount() { return pImpl->GetProcRecordCount(); }

std::string QuickBLAST::GetQuickBLASTOptionString(){
  return pImpl->GetQuickBLASTOptionString();
}
std::string QuickBLAST::Impl::GetQuickBLASTOptionString(){
  return blast_options;
}

ncbi::blast::CBlastOptionsHandle &QuickBLAST::Impl::GetQuickBLASTOptions()
{
  if (!opts) {
    // initialize default options or throw helpful error
    Rcpp::stop("GetQuickBLASTOptions: options handle 'opts' is not initialised");
  }
  return *opts;
}
void QuickBLAST::Impl::ResetHitCount() { hit_count = 0; }
unsigned int QuickBLAST::Impl::GetObjectID()
{
  return obj_id;
}
void QuickBLAST::Impl::SetObjectID(unsigned int id)
{
  obj_id = id;
}

static std::vector<std::pair<std::string, std::string>> BLASTOptionsFromString(const std::string &input)
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

static CRef<ncbi::blast::CBlastOptionsHandle> MakeQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality = CBlastOptions::eLocal, const bool &verbose = true){
  if(program_name.empty()){
    Rcpp::stop("MakeQuickBLASTOptions(): program_name cannot be empty.");
  }
  ncbi::blast::EProgram program = ncbi::blast::ProgramNameToEnum(program_name);
  // Create a CBlastOptionsHandle object
  CRef<ncbi::blast::CBlastOptionsHandle> opts(ncbi::blast::CBlastOptionsFactory::Create(program, locality));
  opts->SetDefaults();
  // Extract the relevant options from the R list and set them in the CBlastOptionsHandle object
  // Example: Extracting and setting the BLAST database
  
  if (options.empty())
  {
    if(verbose)
      Rcpp::Rcout << "Using " << program_name << " Defaults..." << std::endl;
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
      opts->ClearFilterOptions();
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
    }else{
      Rcpp::Rcerr << "Unidentified Option: key: " << pair.first << " and value: " << key_str << std::endl <<std::flush; //DEBUG
      continue;
    }
    if(verbose)
      Rcpp::Rcout << "Option: " << pair.first << " set to : " << key_str << std::endl <<std::flush; //DEBUG
  }
  if(!opts->Validate()){
    Rcpp::stop("MakeQuickBLASTOptions(): Error : Input BLAST options failed validation.");
  }
  return opts;
}

CRef<ncbi::blast::CBlastOptionsHandle> QuickBLAST::Impl::SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality = CBlastOptions::eLocal, const bool& verbose = true)
{
  if(program_name.empty()){
    Rcpp::stop("SetQuickBLASTOptions(): program_name cannot be empty.");
  }
  
  this->opts = MakeQuickBLASTOptions(program_name, options, locality, verbose); 
  
  blast_options = options;
  
  return opts;
}

static Boolean BlastInterruptFn(SBlastProgress* progress) {
  // progress can be nullptr in some calls, so defensively check
  if (!progress || !progress->user_data) return (Boolean)0;
  InterruptContext* ctx = static_cast<InterruptContext*>(progress->user_data);
  // return non-zero (true) to INTERRUPT/stop the BLAST run
  return ctx->stop.load() ? (Boolean)1 : (Boolean)0;
}

QuickBLAST::Impl::Impl(ESeqType seq_type, EStrand strand, std::string program, std::string options, bool save_sequences, bool save_hsp_sequences)
  : seq_type(seq_type), strand(strand), program(program), blast_options(options), save_sequences(save_sequences), save_hsp_sequences(save_hsp_sequences)
{
  try
  {
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    this->num_threads = omp_get_num_threads();
#else
    this->num_threads = 1;
#endif
    arrow_wrapper = std::make_shared<ArrowWrapper>();
    arrow_wrapper->AddFASTAMetadata("program", program);
    arrow_wrapper->AddFASTAMetadata("options", options);
    
    const auto empty_rb_ = arrow::RecordBatch::MakeEmpty(arrow_wrapper->GetBLASTSchema());
    if(!empty_rb_.ok()){
      throw std::runtime_error("QuickBLAST::Impl(): Error creating empty record batch.");
    }
    empty_rb = empty_rb_.ValueOrDie();
    this->save_sequences = save_sequences;
    this->save_hsp_sequences = save_hsp_sequences;
    this->program = program;
    this->opts = MakeQuickBLASTOptions(program, options, CBlastOptions::eLocal, /*verbose*/ true); 
    this->strand = strand;
    this->seq_type = seq_type;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    omp_init_lock(&hit_countLock);
    // omp_init_lock(&cleaner_threadsLock);
#endif
  }catch(const ncbi::CException &e){
    Rcpp::Rcerr << std::string("[QuickBLAST::Impl()]: NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[QuickBLAST::Impl()]: C++ Runtime Error : ") + e.what() << std::endl<< std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[QuickBLAST::Impl()]: Rcpp::exception : ") + e.what() << std::endl<< std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[QuickBLAST::Impl()]: C++ Exception : ") + e.what() << std::endl<< std::flush;
  }catch(...){
    Rcpp::Rcerr << "[QuickBLAST::Impl()]: Unknown Exception" << std::endl<< std::flush;
  }
}

QuickBLAST::Impl::~Impl()
{
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
  omp_destroy_lock(&hit_countLock);
  // omp_destroy_lock(&cleaner_threadsLock);
#endif
  // arrow_wrapper->~ArrowWrapper();
  // DO NOT DELETE NCBI C++ OBJECTs or PTRs or face Corruption
  //  delete self;
  //  opts->ReleaseReference();
  // delete opts;
  // delete arrow_wrapper;
  // Rprintf("~QuickBLAST::Impl ");
}


int QuickBLAST::Impl::GetFrame(int start, int length, ncbi::objects::ENa_strand strand)
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

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractFASTA(const FastaSequenceData &fasta_data)
{
  RcppThread::checkUserInterrupt();
  std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr;
  std::shared_ptr<arrow::Int64Builder> rec_no_builder;
  std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
  rec_no_builder = std::make_shared<arrow::Int64Builder>();
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

static std::string sanitize_protein(const std::string &s_in, bool &had_bad) {
  static const std::string allowed = "ACDEFGHIKLMNPQRSTVWYBXZJUO"; // include common ambiguity letters
  std::string out;
  out.reserve(s_in.size());
  had_bad = false;
  for (unsigned char c0 : s_in) {
    if (c0 == '\0') { had_bad = true; continue; }
    unsigned char c = static_cast<unsigned char>(std::toupper(c0));
    if (c >= 'A' && c <= 'Z') {
      if (allowed.find(c) != std::string::npos) out.push_back(char(c));
      else { had_bad = true; /* optional: map unknown to 'X' */ out.push_back('X'); }
    } else {
      // skip digits, punctuation, whitespace
      had_bad = true;
    }
  }
  return out;
}

static std::string sanitize_nucleotide(const std::string &s_in, bool &had_bad) {
  static const std::string allowed = "ACGTUN*-"; // IUPAC nucleotides
  std::string out; out.reserve(s_in.size()); had_bad = false;
  for (unsigned char c0 : s_in) {
    if (c0 == '\0') { had_bad = true; continue; }
    unsigned char c = static_cast<unsigned char>(std::toupper(c0));
    if (c >= 'A' && c <= 'Z') {
      if (allowed.find(c) != std::string::npos) out.push_back(char(c));
      else { had_bad = true; /* map unknown to N */ out.push_back('N'); }
    } else {
      had_bad = true;
    }
  }
  return out;
}

std::pair<std::shared_ptr<ncbi::blast::SSeqLoc>, ncbi::CSeq_entry_Handle>  QuickBLAST::Impl::CreateSSeqLocFromType(const FastaSequenceData& fasta_data, ncbi::CRef<ncbi::CScope> parent_scope)
{
  // AVOID STRING COPIES: Use the fasta_data fields directly.
  const TSeqPos seqlen = fasta_data.seq.length();
  if(seqlen >= std::numeric_limits<TSeqPos>::max()){
    Rcpp::stop("[CreateSSeqLocFromType()] seqlen >= std::numeric_limits<TSeqPos>::max().");
  }
  
  // DIRECT ID CONSTRUCTION: Pass the header directly to CSeq_id
  CRef<CSeq_id> id(new CSeq_id(
      fasta_data.header, 
      CSeq_id::fParse_RawText | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal
  ));
  
  CRef<CSeq_loc> cseq_loc_obj(new CSeq_loc());
  
  // STREAMLINED LOCATION SETUP: Use SetWhole instead of Select(e_Whole) + SetId(*id)
  cseq_loc_obj->SetWhole(*id);
  
  if (seq_type == ESeqType::eProtein)
  {
    cseq_loc_obj->SetStrand(eNa_strand_unknown);
  }
  else
  {
    switch (strand)
    {
    case EStrand::ePlus:     cseq_loc_obj->SetStrand(eNa_strand_plus); break;
    case EStrand::eMinus:    cseq_loc_obj->SetStrand(eNa_strand_minus); break;
    case EStrand::eUnknown:  cseq_loc_obj->SetStrand(eNa_strand_unknown); break;
    case EStrand::eBoth:     cseq_loc_obj->SetStrand(eNa_strand_both); break;
    case EStrand::eBoth_rev: cseq_loc_obj->SetStrand(eNa_strand_both_rev); break;
    case EStrand::eOther:    cseq_loc_obj->SetStrand(eNa_strand_other); break;
    }
  }
  
  CRef<CSeq_data> seq_data(new CSeq_data());
  
  // ELIMINATE REDUNDANT SELECT: The CIUPAC constructors combined with SetIupac* handle the choice automatically.
  if (seq_type == ESeqType::eProtein) {
    seq_data->SetIupacaa(CIUPACaa(fasta_data.seq));
  } else {
    seq_data->SetIupacna(CIUPACna(fasta_data.seq));
  }
  
  CRef<CSeq_inst> seq_inst(new CSeq_inst());
  seq_inst->SetSeq_data(*seq_data);
  seq_inst->SetLength(seqlen);
  seq_inst->SetRepr(CSeq_inst::eRepr_raw);
  seq_inst->SetTopology(CSeq_inst::eTopology_linear);
  
  // CONDENSED INST SETUP
  seq_inst->SetMol(seq_type == ESeqType::eProtein ? CSeq_inst_Base::eMol_aa : CSeq_inst_Base::eMol_na);
  
  CRef<CBioseq> bioseq(new CBioseq());
  bioseq->SetId().push_back(id);      
  bioseq->SetInst(*seq_inst);
  
  CRef<CSeq_entry> ret_entry(new CSeq_entry());
  ret_entry->SetSeq(*bioseq);
  
  CSeq_entry_Handle tse_handle = parent_scope->AddTopLevelSeqEntry(*ret_entry);
  
  auto sseq = std::make_shared<ncbi::blast::SSeqLoc>(cseq_loc_obj, parent_scope);
  
  // MOVE SEMANTICS: Avoid ref-count bumps by moving handles and smart pointers
  return std::make_pair(std::move(sseq), std::move(tse_handle));
}

std::string QuickBLAST::Impl::GetSSeqLocSequence(const SSeqLoc &seq_loc)
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

SEXP QuickBLAST::Impl::Hits2RList_internal(std::shared_ptr<arrow::Array> array, std::shared_ptr<arrow::DataType> type, const std::string &field_name)
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
    }
    
    return strings;
  }
  else if (type->id() == arrow::Type::INT32)
  {
    
    auto int_array = std::static_pointer_cast<arrow::Int32Array>(array);
    
    Rcpp::IntegerVector ints(int_array->length());
    
    for (int i = 0; i < int_array->length(); ++i)
    {
      
      if (int_array->IsValid(i))
      {
        ints[i] = int_array->Value(i);
      }
    }
    
    return ints;
  }
  else if (type->id() == arrow::Type::INT64)
  {
    
    auto int_array = std::static_pointer_cast<arrow::Int64Array>(array);
    
    Rcpp::IntegerVector ints(int_array->length());
    
    for (int i = 0; i < int_array->length(); ++i)
    {
      
      if (int_array->IsValid(i))
      {
        ints[i] = int_array->Value(i);
      }
    }
    
    return ints;
  }
  else if (type->id() == arrow::Type::DOUBLE)
  {
    auto double_array = std::static_pointer_cast<arrow::DoubleArray>(array);
    Rcpp::NumericVector doubles(double_array->length());
    
    for (int i = 0; i < double_array->length(); ++i)
    {
      
      if (double_array->IsValid(i))
      {
        doubles[field_name] = double_array->Value(i);
      }
    }
    
    return doubles;
  }
  else
  {
    // For other data types that don't have a direct conversion, return R_NilValue (NA)
    return R_NilValue;
  }
}

SEXP QuickBLAST::Impl::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
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

SEXP QuickBLAST::Impl::Hits2RList(const arrow::RecordBatchVector &rb_vector)
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

// template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return arrow_wrapper->SplitFilesIntoEntries(filename, delim, num_threads, Entry_callback, return_values);
}

// template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::StreamFile(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<FastaSequenceData>)> &Entry_callback, bool return_values)
{
  return pImpl->StreamFile(filename, delim, num_threads, Entry_callback, return_values);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, int batch_size, bool verbose) 
{
  
  auto blast_interrupt_ctx = std::make_shared<InterruptContext>();
  safe_jthread interrupt_check_thread(std::thread([this, num_threads, blast_interrupt_ctx](){
    try {
      const unsigned int thread_wait = num_threads > 1 ? 50 : 0;
      while (!blast_interrupt_ctx->stop.load() && quickblast_running.load()) {
        if (RcppThread::isInterrupted()) {  // safe, non-throwing
          blast_interrupt_ctx->stop.store(true);
          break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(thread_wait));
      }
      if(!quickblast_running.load()){
        blast_interrupt_ctx->stop.store(true);
      }
    } catch(...) {
      blast_interrupt_ctx->stop.store(true);
    }
  }));
  // interrupt_check_thread.detach();
  
  try{
    quickblast_running.store(true); 
    unsigned int n_threads = num_threads; 
    SetThreadCount(n_threads);
    
    unsigned int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);
    
    unsigned int s_seq_count = arrow_wrapper->CountCharacter(subjectFile, '>', n_threads);
    
    if((s_seq_count > 10000 || q_seq_count > 10000) && outFile.empty())
      Rcpp::Rcerr << "Warning: Queries/Subjects > 10000, large inputs can crash during return" << std::endl << std::flush; 
    
    arrow_wrapper->AddFASTAMetadata("BLAST locality", "CBlastOptions::eLocal");
    arrow_wrapper->AddFASTAMetadata("Input source", "files");
    arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);
    if (!outfile_sts.ok())
    {
      /* std::cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
                   << outfile_sts.message() << std::endl; */
      // cerr << "ERROR : Could not create output file stream : " << outfile_sts.detail() << std::endl
      //      << outfile_sts.message() << std::endl;
      // REprintf("ERROR : Could not create output file stream : %s \n %s \n", outfile_sts.detail()->ToString().c_str(), outfile_sts.message().c_str());
      Rcpp::Rcerr << std::string("[BLAST_files()] ERROR : Could not create output file stream :") << outfile_sts.detail()->ToString() << std::endl << outfile_sts.message() << std::endl << std::flush;
      quickblast_running.store(false); 
      blast_interrupt_ctx->stop.store(false);
      // interrupt_check_thread.join();
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    const unsigned int totalIterations = q_seq_count * s_seq_count;
    
    if(q_seq_count > n_threads){
      n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
      n_threads = n_threads >= 1 ? n_threads : 1;
    }
    
    if(n_threads > q_seq_count + s_seq_count)
      n_threads = q_seq_count + s_seq_count;
    SetThreadCount(n_threads);
    if(verbose){
      Rcpp::Rcout << "Num Threads: " << n_threads << std::endl << std::flush; //DEBUG
      // std::cout << "BLAST Sequence Limit: " << blast_sequence_limit << std::endl << std::flush; //DEBUG
      Rcpp::Rcout << "Total Records (Q + S): " << q_seq_count + s_seq_count << " (" << q_seq_count << " + " << s_seq_count << ")"<< std::endl << std::flush; //DEBUG
    }
    if(totalIterations <= 0){
      Rcpp::Rcerr << "[BLAST_files()] Improperly formatted FASTA file. No records detected." << std::endl << std::flush;
      quickblast_running.store(false); 
      blast_interrupt_ctx->stop.store(false);
      // interrupt_check_thread.join();
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    if(n_threads > 1 && batch_size <= 0)
      batch_size = n_threads + 1;
    else if(batch_size <= 1)
      batch_size = 2;
    
    arrow_wrapper->SetVerbosity(verbose);
    arrow_wrapper->SetBatchSize(batch_size);
    
    // Progress progress_bar(totalIterations, show_progress);
    // InitProgressBar(totalIterations, show_progress);
    
    RcppThread::ProgressBar progress_bar(totalIterations, verbose);
    
    std::shared_ptr<arrow::RecordBatchVector> final_ret = StreamFile(
      queryFile, ">", n_threads, [this, n_threads, &progress_bar, &blast_interrupt_ctx, subjectFile, batch_size, return_values, verbose](std::shared_ptr<FastaSequenceData> data_q) 
      {
        if(blast_interrupt_ctx->stop.load() || !quickblast_running.load()){
          progress_bar++;
          return std::make_shared<arrow::RecordBatchVector>();
        }
        RcppThread::checkUserInterrupt();
        
        bool had_bad_q = false;
        switch (seq_type)
        {
        case ESeqType::eProtein:
        {
          data_q->seq = sanitize_protein(data_q->seq, had_bad_q);
          break;
        }
        case ESeqType::eNucleotide:
        {
          data_q->seq = sanitize_nucleotide(data_q->seq, had_bad_q);
          break;
        }
        }
        if (had_bad_q) {
          if(verbose)
            Rcpp::Rcout << "Warning: removed/normalized invalid bytes from Query sequence " << data_q->header << std::endl;
        }
        if(data_q->seq.empty()){
          return std::make_shared<arrow::RecordBatchVector>();
        }
        
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_lock_t scopeLock;
        omp_lock_t subjects_loc_vecLock;
        omp_lock_t subjects_seqent_vecLock;
        omp_init_lock(&scopeLock);
        omp_init_lock(&subjects_loc_vecLock);
        omp_init_lock(&subjects_seqent_vecLock);
#endif
        
        CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
        auto [ query_loc, query_seq_entry ] = CreateSSeqLocFromType(*data_q, scope);
        if (!query_loc) {
          Rcpp::Rcerr << "BLAST_files: CreateSSeqLocFromType(query) returned NULL." << std::endl << std::flush;
          return std::make_shared<arrow::RecordBatchVector>();
        }
        if(!query_loc->seqloc.NotEmpty()){
          Rcpp::Rcerr << "BLAST_files: CreateSSeqLocFromType(query) is empty." << std::endl << std::flush;
          return std::make_shared<arrow::RecordBatchVector>();
        }
        
        std::shared_ptr<TSeqLocVector> subjects_loc_vec(new TSeqLocVector());
        std::shared_ptr<vector<CSeq_entry_Handle>> subjects_seqent_vec =  std::make_shared<vector<CSeq_entry_Handle>>();
        
        std::shared_ptr<arrow::RecordBatchVector> ret_results = StreamFile(
          subjectFile, ">", n_threads, [this, query_loc,
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                       &scopeLock, &subjects_loc_vecLock, &subjects_seqent_vecLock,
#endif
                                       &scope, &data_q, &subjects_loc_vec, &subjects_seqent_vec, &progress_bar, &blast_interrupt_ctx, batch_size, return_values, verbose](std::shared_ptr<FastaSequenceData> data_s) //
                                       {
                                         if(blast_interrupt_ctx->stop.load() || !quickblast_running.load()){
                                           progress_bar++;
                                           return std::make_shared<arrow::RecordBatchVector>();
                                         }
                                         
                                         bool had_bad_s = false;
                                         switch (seq_type)
                                         {
                                         case ESeqType::eProtein:
                                         {
                                           data_s->seq = sanitize_protein(data_s->seq, had_bad_s);
                                           break;
                                         }
                                         case ESeqType::eNucleotide:
                                         {
                                           data_s->seq = sanitize_nucleotide(data_s->seq, had_bad_s);
                                           break;
                                         }
                                         }
                                         if (had_bad_s) {
                                           if(verbose)
                                             Rcpp::Rcout << "Warning: removed/normalized invalid bytes from Query sequence " << data_s->header << std::endl << std::flush;
                                         }
                                         if(data_s->seq.empty()){
                                           return std::make_shared<arrow::RecordBatchVector>();
                                         }
                                         
                                         auto [ subject_loc, subject_seq_entry] = CreateSSeqLocFromType(*data_s, scope);
                                         if (!subject_loc) {
                                           Rcpp::Rcerr << "BLAST_files: CreateSSeqLocFromType(subject) returned NULL." << std::endl << std::flush;
                                           return std::make_shared<arrow::RecordBatchVector>();
                                         }
                                         if(!subject_loc->seqloc.NotEmpty()){
                                           Rcpp::Rcerr << "BLAST_files: CreateSSeqLocFromType(subject) is empty." << std::endl << std::flush;
                                           return std::make_shared<arrow::RecordBatchVector>();
                                         }
                                         
                                         progress_bar++;
                                         
                                         if (strcmp(subject_loc->seqloc->GetId()->GetSeqIdString(true).c_str(), query_loc->seqloc->GetId()->GetSeqIdString(true).c_str()) != 0)
                                         {
                                           arrow_wrapper->AddProcRecordCount();
                                           
                                           std::unique_ptr<CBl2Seq> blaster = nullptr;
                                           
                                           try
                                           {
                                             
                                             switch (batch_size)
                                             {
                                             case 0:
                                           {
                                             blaster = std::make_unique<CBl2Seq>(*query_loc, *subject_loc, this->GetQuickBLASTOptions());
                                             blaster->SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
                                             
                                             TSeqAlignVector alignments;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp critical
#endif
{
  try{
    alignments = blaster->Run();
  }catch (const ncbi::CException& e) {
    Rcpp::Rcerr << std::string("[BLAST_files()] 1. BLAST Execution error: Check Sequence type mismatch (proteins != nucleotides) : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
    blast_interrupt_ctx->stop.store(true);      
    quickblast_running.store(false); 
  }
}
// AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);

arrow::RecordBatchVector tmp_rbv = { ExtractHits(alignments, *query_loc, *subject_loc, scope, return_values) }; 

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
omp_set_lock(&scopeLock);
#endif                                
scope->ResetHistory(CScope::EActionIfLocked::eKeepIfLocked);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
omp_unset_lock(&scopeLock);
#endif
subject_loc.reset();

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
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                             omp_set_lock(&subjects_loc_vecLock);
                                             omp_set_lock(&subjects_seqent_vecLock);
#endif
                                             subjects_loc_vec->emplace_back(std::move(*subject_loc));
                                             subjects_seqent_vec->emplace_back(subject_seq_entry);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                             omp_unset_lock(&subjects_loc_vecLock);
                                             omp_unset_lock(&subjects_seqent_vecLock);
#endif
                                             
                                             if (subjects_loc_vec->size() >= batch_size || (arrow_wrapper->GetPendingRecordCount() <= batch_size && subjects_loc_vec->size() > 0))
                                             {
                                               TSeqLocVector subjects_buffer_vec;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                               omp_set_lock(&subjects_loc_vecLock);
                                               
#endif
                                               subjects_buffer_vec.reserve(subjects_loc_vec->size());
                                               subjects_buffer_vec.swap(*subjects_loc_vec);
                                               subjects_loc_vec->clear();
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
                                               omp_unset_lock(&subjects_loc_vecLock);
#endif
                                               blaster = std::make_unique<CBl2Seq>(*query_loc, subjects_buffer_vec, this->GetQuickBLASTOptions(), /*db_scan*/ false);
                                               blaster->SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
                                               
                                               AddHitCount(subjects_buffer_vec.size());
                                               TSeqAlignVector alignments;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp critical
#endif
{
  try{
    alignments = blaster->Run();
  }catch (const ncbi::CException& e) {
    Rcpp::Rcerr << std::string("[BLAST_files()] 2. BLAST Execution error: Check Sequence type mismatch (proteins != nucleotides) : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
    blast_interrupt_ctx->stop.store(true);      
    quickblast_running.store(false); 
  }
}
// AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);

std::shared_ptr<arrow::RecordBatchVector> tmp_rbv = ExtractHits(alignments, *query_loc, subjects_buffer_vec, scope, return_values); 

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
omp_set_lock(&scopeLock);
#endif 
scope->ResetHistory(CScope::EActionIfLocked::eKeepIfLocked);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
omp_unset_lock(&scopeLock);
#endif 
subjects_buffer_vec.clear();
subjects_buffer_vec.shrink_to_fit();

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
                                           }catch (const CException& e){
                                             Rcpp::Rcerr << std::string("[BLAST_files()]: 1. NCBI Exception: Stopping run :")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << std::flush;
                                             blast_interrupt_ctx->stop.store(true);      
                                             quickblast_running.store(false); 
                                           }catch(const Rcpp::exception &e){
                                             Rcpp::Rcerr << std::string("[BLAST_files()] - 1. Rcpp Exception : ") + e.what() << std::endl << std::flush;
                                             blast_interrupt_ctx->stop.store(true);      
                                             quickblast_running.store(false); 
                                           }catch (const std::runtime_error& e){
                                             Rcpp::Rcerr << std::string("[BLAST_files()]: 1. C++ Runtime Error: ") + e.what() << std::endl << std::flush;
                                             blast_interrupt_ctx->stop.store(true);      
                                             quickblast_running.store(false); 
                                           }catch (const std::exception& e){
                                             Rcpp::Rcerr << std::string("[BLAST_files()]: 1. C++ Exception: ") + e.what() << std::endl << std::flush;
                                             blast_interrupt_ctx->stop.store(true);      
                                             quickblast_running.store(false); 
                                           }catch(...){
                                             Rcpp::Rcerr << "[BLAST_files()] - 1. Unknown Exception" << std::endl << std::flush;
                                             blast_interrupt_ctx->stop.store(true);      
                                             quickblast_running.store(false); 
                                           }
                                           
                                           
                                         }
                                         
                                         return std::make_shared<arrow::RecordBatchVector>(); // EMPTY ERROR Return
                                       },
                                       return_values);
        
        if (subjects_loc_vec->size() > 0)
        {
          std::unique_ptr<CBl2Seq> blaster = std::make_unique<CBl2Seq>(*query_loc, *subjects_loc_vec, this->GetQuickBLASTOptions(), /*db_scan*/ false);
          blaster->SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
          AddHitCount(subjects_loc_vec->size());
          TSeqAlignVector alignments;
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp critical
#endif
{
  try{
    alignments = blaster->Run();
  }catch (const ncbi::CException& e) {
    Rcpp::Rcerr << std::string("[BLAST_files()] 3. BLAST Execution error: Check Sequence type mismatch (proteins != nucleotides): ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << std::flush;
    quickblast_running.store(false); 
  }
}
// AddAllAvailableScoresToSeqAlignVector(alignments, scope, 0);

std::shared_ptr<arrow::RecordBatchVector> ret_vec = ExtractHits(alignments, *query_loc, *subjects_loc_vec, scope, return_values); 

subjects_loc_vec->clear();
subjects_loc_vec->shrink_to_fit();

if (return_values)
{
  ret_results->insert(ret_results->end(), ret_vec->begin(), ret_vec->end());
}
else
{
  ret_vec->clear();
}
        }
        
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_set_lock(&scopeLock);
#endif 
        scope->ResetHistory(CScope::EActionIfLocked::eRemoveIfLocked);
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_unset_lock(&scopeLock);
#endif 
        
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
        omp_destroy_lock(&scopeLock);
        omp_destroy_lock(&subjects_loc_vecLock);
        omp_destroy_lock(&subjects_seqent_vecLock);
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
    
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp barrier
#endif
    
    // std::cout << "Final Batch Size: " << arrow_wrapper->GetBatchSize() << std::endl << std::flush;  //DEBUG
    
    static_cast<void>(arrow_wrapper->FinishOutputStream());
    
    if(verbose)
      Rcpp::Rcout << "Total Records Processed: " << arrow_wrapper->GetProcRecordCount() << std::endl << std::flush;
    
    arrow_wrapper->ResetProcRecordCount();
    quickblast_running.store(false); 
    blast_interrupt_ctx->stop.store(false);
    // interrupt_check_thread.join();
    
    if (return_values)
    {
      return final_ret;
    }
    else
    {
      final_ret->clear();
      final_ret->shrink_to_fit();
      return std::make_shared<arrow::RecordBatchVector>();
    }
  }catch (const ncbi::CException &e) {
    // NCBI toolkit exceptions
    Rcpp::Rcerr << std::string("[BLAST_files()]: 2. NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << std::flush;
  }catch (const std::runtime_error &e) {
    std::string msg = "[BLAST_files()]: 2. C++ Runtime Error : ";
    Rcpp::Rcerr << msg + e.what() << std::endl << std::flush; 
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[BLAST_files()] - 2. Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch (const std::exception &e) {
    Rcpp::Rcerr << std::string("[BLAST_files()] - 2. C++ exception : ") + e.what() << std::endl << std::flush;
  }catch (...) {
    Rcpp::Rcerr << "[BLAST_files()]: 2. Unknown exception" << std::endl << std::flush;
  }
  static_cast<void>(arrow_wrapper->FinishOutputStream());
  arrow_wrapper->ResetProcRecordCount();
  quickblast_running.store(false); 
  // interrupt_check_thread.join();
  return std::make_shared<arrow::RecordBatchVector>();
}

static inline int GetSubjectOID(const CSeqDB& db,
                                const CSeq_id& id)
{
  vector<int> oids;
  db.SeqidToOids(id, oids);
  if (oids.empty())
    NCBI_THROW(CException, eUnknown, "Seq-id not in DB");
  return oids.front();
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_dbs(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, unsigned int batch_size, const bool& enable_chunking, unsigned int chunk_size, unsigned int overlap, bool verbose){ //const bool show_progress
  return pImpl->BLAST_dbs(queryFile, subjectFile, outFile, outFormat, num_threads, return_values, batch_size, enable_chunking, chunk_size, overlap, verbose); 
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_f2db(const std::string &queryFile, const std::string &subjectDB, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, unsigned int batch_size, bool verbose) {
  return pImpl->BLAST_f2db(queryFile, subjectDB, outFile, outFormat, num_threads, return_values, batch_size, verbose); 
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_f2db(const std::string &queryFile, const std::string &subjectDB, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, unsigned int batch_size, bool verbose) {
  auto blast_interrupt_ctx = std::make_shared<InterruptContext>();
  safe_jthread interrupt_check_thread(std::thread([this, num_threads, blast_interrupt_ctx](){
    try {
      const unsigned int thread_wait = num_threads > 1 ? 50 : 0;
      while (!blast_interrupt_ctx->stop.load() && quickblast_running.load()) {
        if (RcppThread::isInterrupted()) {  // safe, non-throwing
          blast_interrupt_ctx->stop.store(true);
          break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(thread_wait));
      }
      if(!quickblast_running.load()){
        blast_interrupt_ctx->stop.store(true);
      }
    } catch(...) {
      blast_interrupt_ctx->stop.store(true);
    }
  }));
  try{
    quickblast_running.store(true); 
    
    unsigned int n_threads = num_threads; 
    SetThreadCount(n_threads);
    
    CSeqDB::ESeqType seqdbType;
    CSearchDatabase::EMoleculeType seqType;
    Rcpp::Rcout << seq_type << std::endl << std::flush; //DEBUG
    switch(seq_type){
    case QuickBLAST::ESeqType::eNucleotide: {
      seqType = CSearchDatabase::EMoleculeType::eBlastDbIsNucleotide;
      seqdbType = CSeqDB::eNucleotide;
      break;
    }
    case QuickBLAST::ESeqType::eProtein: {
      seqType = CSearchDatabase::EMoleculeType::eBlastDbIsProtein;
      seqdbType = CSeqDB::eProtein;
      break;
    }
    }
    
    CRef<CSeqDB> s_seqdb_(new CSeqDB(subjectDB, seqdbType));
    
    
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    s_seqdb_->SetNumberOfThreads(n_threads, /*force_mt*/ true);
#else
    s_seqdb_->SetNumberOfThreads(1, /*force_mt*/ false);
#endif
    
    CRef<CSearchDatabase> s_serdb_(new CSearchDatabase(subjectDB, seqType));
    
    unsigned int q_seq_count = arrow_wrapper->CountCharacter(queryFile, '>', n_threads);
    
    unsigned int s_seq_count = s_seqdb_->GetNumSeqs(); //arrow_wrapper->CountCharacter(subjectDB, '>', n_threads);
    
    if((q_seq_count > 10000) && outFile.empty() && verbose)
      Rcpp::Rcerr << "Warning: Queries > 10000, large inputs can crash during return" << std::endl << std::flush; 
    
    arrow_wrapper->AddFASTAMetadata("BLAST locality", "CBlastOptions::eLocal");
    arrow_wrapper->AddFASTAMetadata("Input source", "file2DB");
    arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);
    if (!outfile_sts.ok())
    {
      quickblast_running.store(false);
      // interrupt_check_thread.join();
      Rcpp::Rcerr << std::string("[BLAST_f2db()] ERROR : Could not create output file stream ") + outfile_sts.detail()->ToString() + std::string("\n") + outfile_sts.message() << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    const unsigned int totalIterations = q_seq_count * s_seq_count;
    arrow_wrapper->SetBLASTSeqLimit(batch_size);
    
    if(q_seq_count > n_threads){
      n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
      n_threads = n_threads >= 1 ? n_threads : 1;
    }
    
    if(n_threads > q_seq_count + s_seq_count)
      n_threads = q_seq_count + s_seq_count;
    SetThreadCount(n_threads);
    
    if(verbose){
      Rcpp::Rcout << "Num Threads: " << n_threads << std::endl << std::flush; //DEBUG
      // std::cout << "BLAST Sequence Limit: " << blast_sequence_limit << std::endl << std::flush; //DEBUG
      Rcpp::Rcout << "Total Records (Q + S): " << q_seq_count + s_seq_count << " (" << q_seq_count << " + " << s_seq_count << ")"<< std::endl << std::flush; //DEBUG
    }
    
    if(totalIterations <= 0){
      quickblast_running.store(false);
      // interrupt_check_thread.join();
      Rcpp::Rcerr << "[BLAST_f2db()] Improperly formatted FASTA file. No records detected." << std::endl << std::flush;
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    if(n_threads > 1 && batch_size <= 0)
      batch_size = n_threads + 1;
    else if(batch_size <= 1)
      batch_size = 2;
    
    arrow_wrapper->SetVerbosity(verbose);
    arrow_wrapper->SetBatchSize(batch_size);
    
    // Progress progress_bar(totalIterations, show_progress);
    // InitProgressBar(totalIterations, show_progress);
    
    RcppThread::ProgressBar progress_bar(totalIterations, verbose);
    
    std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> file_ptrs = arrow_wrapper->MMapFile(queryFile, ">");
    
    char *start_of_file = std::get<1>(*file_ptrs).get();
    char *end_of_file = std::get<3>(*file_ptrs);
    char* delim = ">";
    unsigned int rec_no = 1;
    
    std::shared_ptr<arrow::RecordBatchVector> final_ret = std::make_shared<arrow::RecordBatchVector>();
    
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp parallel shared(final_ret, file_ptrs, q_seq_count, s_seq_count, progress_bar, delim, blast_interrupt_ctx, subjectDB, seqType)
#endif
{ 
  try{
    // Thread-Local Builders (Hoisted for memory reuse)
    arrow::Int64Builder hsp_offset_builder, length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder, qlen_builder, slen_builder, num_alignments_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    arrow::StringBuilder frames_builder, strand_builder, qseqid_builder, sseqid_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder, qhsp_builder, shsp_builder;
    
    std::shared_ptr<arrow::RecordBatchVector> local_ret = std::make_shared<arrow::RecordBatchVector>();
    
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp for schedule(dynamic)
#endif
    for(unsigned int rec_no = 1; rec_no < q_seq_count; rec_no += batch_size)
    {
      if(!quickblast_running.load()) continue;
      if(rec_no > q_seq_count) continue;
      
      std::shared_ptr<std::list<FastaSequenceData>> fasta_batch = arrow_wrapper->FetchRecordByBatch(file_ptrs, batch_size, rec_no, delim);
      CRef<ncbi::CScope> scope(new ncbi::CScope(*CObjectManager::GetInstance()));
      scope->AddDefaults();
      CRef<CBlastQueryVector> blastquery_batch(new CBlastQueryVector());
      CRef<CSearchDatabase> lcl_s_serdb(new CSearchDatabase(subjectDB, seqType));
      CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
      CRef<ncbi::blast::CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, /*verbose*/ false);
      lcl_blast_opts->SetDbLength(s_seq_count);
      
      for(FastaSequenceData data_q : *fasta_batch) {
        RcppThread::checkUserInterrupt();
        
        bool had_bad_q = false;
        if (seq_type == ESeqType::eProtein) data_q.seq = sanitize_protein(data_q.seq, had_bad_q);
        else data_q.seq = sanitize_nucleotide(data_q.seq, had_bad_q);
        
        if(data_q.seq.empty()) continue;
        
        CRef<CSeq_id> id(new CSeq_id(std::string(data_q.header.data()), CSeq_id::fParse_RawText | CSeq_id::fParse_PartialOK | CSeq_id::fParse_AnyLocal));
        
        CRef<CSeq_loc> cseq_loc_obj(new CSeq_loc());
        cseq_loc_obj->SetWhole(*id);
        
        CRef<CSeq_data> seq_data(new CSeq_data());
        if (seq_type == ESeqType::eProtein) seq_data->SetIupacaa(CIUPACaa(std::string(data_q.seq.data())));
        else seq_data->SetIupacna(CIUPACna(std::string(data_q.seq.data())));
        
        CRef<CSeq_inst> seq_inst(new CSeq_inst());
        seq_inst->SetSeq_data(*seq_data);
        seq_inst->SetLength(data_q.seq.length());
        seq_inst->SetRepr(CSeq_inst::eRepr_raw);
        seq_inst->SetTopology(CSeq_inst::eTopology_linear);
        seq_inst->SetMol(seq_type == ESeqType::eProtein ? CSeq_inst_Base::eMol_aa : CSeq_inst_Base::eMol_na);
        
        CRef<CBioseq> q_bioseq(new CBioseq());
        q_bioseq->SetId().push_back(id);
        q_bioseq->SetInst(*seq_inst);
        
        CRef<CSeq_entry> q_seqentry(new CSeq_entry());
        q_seqentry->SetSeq(*q_bioseq);
        scope->AddTopLevelSeqEntry(*q_seqentry);
        
        CRef<CBlastSearchQuery> q(new CBlastSearchQuery(*cseq_loc_obj, *scope));
        blastquery_batch->AddQuery(q);
      }
      
      if (blastquery_batch->Empty()) continue;
      CRef<CSearchResultSet> results;
      try {
        CRef<ncbi::blast::CObjMgr_QueryFactory> query_factory(new ncbi::blast::CObjMgr_QueryFactory(*blastquery_batch));
        CLocalBlast lcl_blaster(query_factory, lcl_blast_opts, s_db_adapter);  
        lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
        results = lcl_blaster.Run();
      } catch (const ncbi::CException& e) {
        Rcpp::Rcerr << "[BLAST_f2db()] Execution error: " << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
        quickblast_running.store(false);
        continue;
      }
      
      for (int pbi = 0; pbi < batch_size; pbi++) {  
        progress_bar++;
      }
      
      if(results->empty()) continue;
      
      // --- PRE-COMPUTE TOTAL ALIGNMENTS ---
      int64_t total_hsp_count = 0;
      for (const auto& res_ref : *results) {
        if (res_ref && res_ref->HasAlignments()) {
          auto seq_align_set = res_ref->GetSeqAlign();
          if (seq_align_set && !seq_align_set->IsEmpty()) {
            total_hsp_count += seq_align_set->Get().size();
          }
        }
      }
      
      if (total_hsp_count == 0) continue;
      
      // --- RESERVE CAPACITY ---
      static_cast<void>(hsp_offset_builder.Reserve(total_hsp_count));
      static_cast<void>(length_builder.Reserve(total_hsp_count));
      static_cast<void>(mismatch_builder.Reserve(total_hsp_count));
      static_cast<void>(gapopen_builder.Reserve(total_hsp_count));
      static_cast<void>(qstart_builder.Reserve(total_hsp_count));
      static_cast<void>(qend_builder.Reserve(total_hsp_count));
      static_cast<void>(sstart_builder.Reserve(total_hsp_count));
      static_cast<void>(send_builder.Reserve(total_hsp_count));
      static_cast<void>(gaps_builder.Reserve(total_hsp_count));
      static_cast<void>(nident_builder.Reserve(total_hsp_count));
      static_cast<void>(positive_builder.Reserve(total_hsp_count));
      static_cast<void>(n_splices_builder.Reserve(total_hsp_count));
      static_cast<void>(hsp_cnt_builder.Reserve(total_hsp_count));
      static_cast<void>(negative_count_builder.Reserve(total_hsp_count));
      static_cast<void>(pident_builder.Reserve(total_hsp_count));
      static_cast<void>(pident_gap_builder.Reserve(total_hsp_count));
      static_cast<void>(evalue_builder.Reserve(total_hsp_count));
      static_cast<void>(bitscore_builder.Reserve(total_hsp_count));
      static_cast<void>(score_builder.Reserve(total_hsp_count));
      static_cast<void>(qcovhsp_builder.Reserve(total_hsp_count));
      static_cast<void>(blast_score_builder.Reserve(total_hsp_count));
      static_cast<void>(aln_len01_builder.Reserve(total_hsp_count));
      static_cast<void>(sum_evalue_builder.Reserve(total_hsp_count));
      static_cast<void>(product_coverage_builder.Reserve(total_hsp_count));
      static_cast<void>(overall_identity_builder.Reserve(total_hsp_count));
      static_cast<void>(matches_builder.Reserve(total_hsp_count));
      static_cast<void>(high_quality_percent_coverage_builder.Reserve(total_hsp_count));
      static_cast<void>(exon_identity_builder.Reserve(total_hsp_count));
      static_cast<void>(consensus_splices_builder.Reserve(total_hsp_count));
      static_cast<void>(comp_adj_method_builder.Reserve(total_hsp_count));
      static_cast<void>(qlen_builder.Reserve(total_hsp_count));
      static_cast<void>(slen_builder.Reserve(total_hsp_count));
      static_cast<void>(num_alignments_builder.Reserve(total_hsp_count));
      
      static_cast<void>(frames_builder.Reserve(total_hsp_count));
      static_cast<void>(strand_builder.Reserve(total_hsp_count));
      static_cast<void>(qseqid_builder.Reserve(total_hsp_count));
      static_cast<void>(sseqid_builder.Reserve(total_hsp_count));
      static_cast<void>(qseq_builder.Reserve(total_hsp_count));
      static_cast<void>(sseq_builder.Reserve(total_hsp_count));
      static_cast<void>(qhsp_builder.Reserve(total_hsp_count));
      static_cast<void>(shsp_builder.Reserve(total_hsp_count));
      
      // --- PRE-ALLOCATE REUSABLE STRINGS ---
      std::string q_full, s_full, strand_str, frames, q_aligned, s_aligned;
      q_full.reserve(4096); s_full.reserve(4096);
      q_aligned.reserve(1024); s_aligned.reserve(1024);
      
      int num_rows = 0;
      
      for (const auto &res_ref : *results) {
        if (!res_ref || !res_ref->HasAlignments()) continue;
        
        const auto seq_align_set = res_ref->GetSeqAlign();
        if (!seq_align_set || seq_align_set->IsEmpty()) continue;
        
        std::string qseq_id = res_ref->GetSeqId()->GetSeqIdString(true);
        auto seq_aligns = seq_align_set->Get();
        int64_t parent_list_size = seq_aligns.size();
        ncbi::sequence::CDeflineGenerator defline_gen;
        q_full.clear(); s_full.clear();
        
        if (save_sequences) {
          ncbi::objects::CBioseq_Handle bh = scope->GetBioseqHandle(*res_ref->GetSeqId());
          if (bh) {
            ncbi::objects::CSeqVector v = bh.GetSeqVector();
            v.SetIupacCoding();
            v.GetSeqData(0, v.size(), q_full);
          }
        }
        
        for (const auto &it : seq_aligns) {
          if (!it || it.IsNull() || !it.NotEmpty() || !it->CanGetScore()) continue;
          
          q_aligned.clear(); s_aligned.clear();
          
          std::string sseq_id;
          CBioseq_Handle s_bh = scope->GetBioseqHandle(it->GetSeq_id(1));
          if (s_bh) {
            // 1. Pull the original FASTA line from the Bioseq's title/description
            std::string full_header = defline_gen.GenerateDefline(s_bh);
            // 2. FASTA headers often have a structure like ">ID Description..."
            // If you ONLY want the "ID" part (the first word), split it at the first space.
            size_t space_pos = full_header.find_first_of(" \t");
            if (space_pos != std::string::npos) {
              sseq_id = full_header.substr(0, space_pos);
            } else {
              sseq_id = full_header;
            }
            // 3. Fallback just in case the title is truly empty
            if (sseq_id.empty()) {
              sseq_id = it->GetSeq_id(1).GetSeqIdString(true);
            }
          } else {
            sseq_id = it->GetSeq_id(1).GetSeqIdString(true);
          }
          
          // Use Append for strings
          static_cast<void>(qseqid_builder.Append(qseq_id));
          static_cast<void>(sseqid_builder.Append(sseq_id));
          
          TSeqPos qlen = scope->GetBioseqHandle(it->GetSeq_id(0)).GetBioseqLength();
          TSeqPos slen = scope->GetBioseqHandle(it->GetSeq_id(1)).GetBioseqLength();
          static_cast<void>(qlen_builder.UnsafeAppend(qlen));
          static_cast<void>(slen_builder.UnsafeAppend(slen));
          
          ncbi::objects::ENa_strand q_strand = it->GetSeqStrand(0);
          ncbi::objects::ENa_strand s_strand = it->GetSeqStrand(1);
          char q_strand_char = (q_strand == ncbi::objects::eNa_strand_plus) ? '+' : (q_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
          char s_strand_char = (s_strand == ncbi::objects::eNa_strand_plus) ? '+' : (s_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
          strand_str = std::string(1, q_strand_char) + "/" + s_strand_char;
          
          if (it->GetSegs().IsDenseg()) {
            const auto& dseg = it->GetSegs().GetDenseg();
            if (save_sequences && dseg.CanGetIds()){
              if (dseg.GetIds().size() > 1)
                GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[1]), scope, s_full);
            }
            if (save_hsp_sequences) {
              // GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned);
              GetHSPSequencesGracefully(dseg, scope, q_aligned, s_aligned);
            }
          }
          
          double score=0, n_splices=0, num_ident=0, aln_len=0, gaps=0, mismatches=0, positive=0, negative_count=0;
          double bits=0, evalue=0, blast_score=0, pident=0, aln_len01=0, pident_gap=0, qcovhsp=0;
          double sum_evalue=0, product_coverage=0, overall_identity=0, matches=0, high_quality_percent_coverage=0, exon_identity=0, consensus_splices=0, comp_adj_method=0;
          
          if(!it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len)) aln_len = it->GetAlignLength(true);
          it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); 
          it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
          bool hasid = it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
          
          if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident) && hasid) pident = 100.0 * num_ident / it->GetAlignLength(false); 
          if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap) && hasid) pident_gap = 100.0 * num_ident / it->GetAlignLength(true); 
          if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps)){ 
            // gaps = it->GetTotalGapCount(-1);
            // Proactively check the segment type to prevent NCBI_THROW from logging errors
            const auto& segs = it->GetSegs();
            
            if (segs.IsDenseg() || segs.IsPacked()) {
              // Standard alignments (blastn, blastp, blastx) support this natively
              gaps = it->GetTotalGapCount(-1);
            } else {
              // tblastx uses Std-seg or Dendiag which represent ungapped blocks.
              // Since there are no gaps in these specific alignment formats, it is exactly 0.
              gaps = 0;
            }
          }
          it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); 
          if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches)) mismatches = it->GetAlignLength(true) - num_ident - gaps;
          
          qcovhsp = q_full.length() > 0 ? (static_cast<double>(it->GetAlignLength(false)) / static_cast<double>(q_full.length())) : 0.0;
          
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
          int qstart = it->GetSeqStart(0);
          int qend = it->GetSeqStop(0);
          int sstart = it->GetSeqStart(1);
          int send = it->GetSeqStop(1);
          
          frames = std::to_string(GetFrame(qstart, aln_len, q_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, s_strand));
          
          // Append strings
          static_cast<void>(qhsp_builder.Append(q_aligned));
          static_cast<void>(shsp_builder.Append(s_aligned));
          static_cast<void>(frames_builder.Append(frames));
          static_cast<void>(strand_builder.Append(strand_str));
          static_cast<void>(qseq_builder.Append(save_sequences ? q_full : ""));
          static_cast<void>(sseq_builder.Append(save_sequences ? s_full : ""));
          
          // UnsafeAppend primitives
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
          static_cast<void>(num_alignments_builder.UnsafeAppend(parent_list_size));
          static_cast<void>(hsp_offset_builder.UnsafeAppend(1));
          
          num_rows++;
        }
      }
      
      if (num_rows > 0) {
        std::shared_ptr<arrow::Array> qhsp_array, shsp_array, frames_array, pident_array, pident_gap_array, evalue_array, length_array, qstart_array, qend_array, sstart_array, send_array, aln_len01_array, bitscore_array, score_array, qcovhsp_array, blast_score_array, gaps_array, nident_array, mismatch_array, positive_array, n_splices_array, hsp_cnt_array, sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array, matches_array, high_quality_percent_coverage_array, exon_identity_array, consensus_splices_array, comp_adj_method_array;
        
        static_cast<void>(qhsp_builder.Finish(&qhsp_array));
        static_cast<void>(shsp_builder.Finish(&shsp_array));
        static_cast<void>(frames_builder.Finish(&frames_array));
        static_cast<void>(pident_builder.Finish(&pident_array));
        static_cast<void>(pident_gap_builder.Finish(&pident_gap_array));
        static_cast<void>(evalue_builder.Finish(&evalue_array));
        static_cast<void>(length_builder.Finish(&length_array));
        static_cast<void>(qstart_builder.Finish(&qstart_array));
        static_cast<void>(qend_builder.Finish(&qend_array));
        static_cast<void>(sstart_builder.Finish(&sstart_array));
        static_cast<void>(send_builder.Finish(&send_array));
        static_cast<void>(aln_len01_builder.Finish(&aln_len01_array));
        static_cast<void>(bitscore_builder.Finish(&bitscore_array));
        static_cast<void>(score_builder.Finish(&score_array));
        static_cast<void>(qcovhsp_builder.Finish(&qcovhsp_array));
        static_cast<void>(blast_score_builder.Finish(&blast_score_array));
        static_cast<void>(gaps_builder.Finish(&gaps_array));
        static_cast<void>(nident_builder.Finish(&nident_array));
        static_cast<void>(mismatch_builder.Finish(&mismatch_array));
        static_cast<void>(positive_builder.Finish(&positive_array));
        static_cast<void>(n_splices_builder.Finish(&n_splices_array));
        static_cast<void>(hsp_cnt_builder.Finish(&hsp_cnt_array));
        static_cast<void>(sum_evalue_builder.Finish(&sum_evalue_array));
        static_cast<void>(product_coverage_builder.Finish(&product_coverage_array));
        static_cast<void>(overall_identity_builder.Finish(&overall_identity_array));
        static_cast<void>(negative_count_builder.Finish(&negative_count_array));
        static_cast<void>(matches_builder.Finish(&matches_array));
        static_cast<void>(high_quality_percent_coverage_builder.Finish(&high_quality_percent_coverage_array));
        static_cast<void>(exon_identity_builder.Finish(&exon_identity_array));
        static_cast<void>(consensus_splices_builder.Finish(&consensus_splices_array));
        static_cast<void>(comp_adj_method_builder.Finish(&comp_adj_method_array));
        
        arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({
          qhsp_array, shsp_array, pident_array, pident_gap_array, frames_array, evalue_array, length_array, aln_len01_array, qstart_array, qend_array, sstart_array, send_array, bitscore_array, score_array, qcovhsp_array, blast_score_array, gaps_array, nident_array, mismatch_array, positive_array, n_splices_array, hsp_cnt_array, sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array, matches_array, high_quality_percent_coverage_array, exon_identity_array, consensus_splices_array, comp_adj_method_array},
          {"qhsp", "shsp", "pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});
        
        if (!aln_struct_array.ok()) {
          Rcpp::Rcerr << "[BLAST_f2db()] 1. Failed to build StructArray: " << aln_struct_array.status().ToString() << std::endl << std::flush;
          quickblast_running.store(false); 
          continue;
        }
        std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
        
        std::shared_ptr<arrow::Array> qseqid_array, sseqid_array, qseq_array, sseq_array, qlen_array, slen_array, strand_array, num_alignment_array;
        static_cast<void>(qseqid_builder.Finish(&qseqid_array));
        static_cast<void>(sseqid_builder.Finish(&sseqid_array));
        static_cast<void>(qseq_builder.Finish(&qseq_array));
        static_cast<void>(sseq_builder.Finish(&sseq_array));
        static_cast<void>(qlen_builder.Finish(&qlen_array));
        static_cast<void>(slen_builder.Finish(&slen_array));
        static_cast<void>(strand_builder.Finish(&strand_array));
        static_cast<void>(num_alignments_builder.Finish(&num_alignment_array));
        
        std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
        std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
        std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
        
        arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make(
        {num_alignment_array, seqids_struct_array, seqs_struct_array, strand_array, lengths_struct_array},
        {"num_alignments", "seqids", "seqs", "strands", "lengths"});
        
        if(!seq_info_array.ok()){
          Rcpp::Rcerr << "[BLAST_f2db()] 2. Failed to build StructArray: " << seq_info_array.status().ToString() << std::endl << std::flush;
          quickblast_running.store(false); 
          continue;
        }
        
        std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
        
        std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(), num_rows, {seq_info_array_, aln_struct_array_});
        
        if(alignment_rb && alignment_rb->ValidateFull().ok()){
          local_ret->emplace_back(alignment_rb);
        }
        
        // Reset Builders
        hsp_offset_builder.Reset(); length_builder.Reset(); mismatch_builder.Reset(); gapopen_builder.Reset(); qstart_builder.Reset(); qend_builder.Reset(); sstart_builder.Reset(); send_builder.Reset(); gaps_builder.Reset(); nident_builder.Reset(); positive_builder.Reset(); n_splices_builder.Reset(); hsp_cnt_builder.Reset(); negative_count_builder.Reset(); pident_builder.Reset(); pident_gap_builder.Reset(); evalue_builder.Reset(); bitscore_builder.Reset(); score_builder.Reset(); qcovhsp_builder.Reset(); blast_score_builder.Reset(); aln_len01_builder.Reset(); sum_evalue_builder.Reset(); product_coverage_builder.Reset(); overall_identity_builder.Reset(); matches_builder.Reset(); high_quality_percent_coverage_builder.Reset(); exon_identity_builder.Reset(); consensus_splices_builder.Reset(); comp_adj_method_builder.Reset(); frames_builder.Reset(); strand_builder.Reset(); qseqid_builder.Reset(); sseqid_builder.Reset(); qseq_builder.Reset(); sseq_builder.Reset(); qhsp_builder.Reset(); shsp_builder.Reset(); qlen_builder.Reset(); slen_builder.Reset(); num_alignments_builder.Reset();
      }
      
      scope->ResetHistory(CScope::EActionIfLocked::eRemoveIfLocked);
    }
    
    // Aggregate local batches safely
    if (!local_ret->empty()) {
#pragma omp critical(final_ret_insert)
{
  final_ret->insert(final_ret->end(), std::make_move_iterator(local_ret->begin()), std::make_move_iterator(local_ret->end()));
}
    }
    
  }catch (const ncbi::CException& e) {
    Rcpp::Rcerr << "[BLAST_f2db()] 1. NCBI CException: " << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch (const std::runtime_error& e) {
    Rcpp::Rcerr << "[BLAST_f2db()] 1. C++ Runtime error: " << e.what() << std::endl << std::flush;
  }catch (const std::exception& e) {
    Rcpp::Rcerr << "[BLAST_f2db()] 1. C++ Exception: " << e.what() << std::endl << std::flush;
  }catch (...) {
    Rcpp::Rcerr << "[BLAST_f2db()] 1. Unknown error" << std::endl << std::flush;
  }
  
}

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp barrier
#endif


static_cast<void>(arrow_wrapper->FinishOutputStream());

if(verbose)
  Rcpp::Rcout << "Total Records Processed: " << arrow_wrapper->GetProcRecordCount() << std::endl << std::flush;  

arrow_wrapper->ResetProcRecordCount();
quickblast_running.store(false); 
// interrupt_check_thread.join();
if (return_values)
{
  return final_ret;
}
else
{
  final_ret->clear();
  final_ret->shrink_to_fit();
  return std::make_shared<arrow::RecordBatchVector>();
}
  }catch (const ncbi::CException &e) {
    // NCBI toolkit exceptions
    Rcpp::Rcerr << std::string("[BLAST_f2db()]: 2. NCBI CException :")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch (const std::runtime_error &e) {
    std::string msg = "[BLAST_f2db()]: 2. C++ Runtime Error : ";
    Rcpp::Rcerr << msg + e.what() << std::endl << std::flush; 
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[BLAST_f2db()] - 2. Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch (const std::exception &e) {
    Rcpp::Rcerr << std::string("[BLAST_f2db()] - 2. C++ exception : ") + e.what() << std::endl << std::flush;
  }catch (...) {
    Rcpp::Rcerr <<"[BLAST_f2db()]: 2. Unknown exception" << std::endl << std::flush;
  }
  static_cast<void>(arrow_wrapper->FinishOutputStream());
  arrow_wrapper->ResetProcRecordCount();
  quickblast_running.store(false);
  // interrupt_check_thread.join();
  return std::make_shared<arrow::RecordBatchVector>();
}

unsigned int QuickBLAST::SizeOfDB(const std::string &dbName){
  return pImpl->SizeOfDB(dbName);
}

unsigned int QuickBLAST::Impl::SizeOfDB(const std::string &dbName){
  CSeqDB::ESeqType seqdbType;
  CSearchDatabase::EMoleculeType seqType;
  switch(seq_type){
  case QuickBLAST::ESeqType::eNucleotide: {
    seqType = CSearchDatabase::EMoleculeType::eBlastDbIsNucleotide;
    seqdbType = CSeqDB::eNucleotide;
    break;
  }
  case QuickBLAST::ESeqType::eProtein: {
    seqType = CSearchDatabase::EMoleculeType::eBlastDbIsProtein;
    seqdbType = CSeqDB::eProtein;
    break;
  }
  }
  
  CRef<CSeqDB> seqdb(new CSeqDB(dbName, seqdbType));
  return seqdb->GetNumSeqs();
}

void QuickBLAST::Impl::GetHSPSequencesGracefully(const CDense_seg& dseg, 
                                                 CRef<CScope> scope,
                                                 string& q_aligned, 
                                                 string& s_aligned) 
{
  // CAlnVec handles the mapping between coordinates and gaps for you
  CAlnVec alnVec(dseg, *scope);
  alnVec.SetGapChar('-');
  
  // row 0 is query, row 1 is subject
  alnVec.GetWholeAlnSeqString(0, q_aligned);
  alnVec.GetWholeAlnSeqString(1, s_aligned);
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractHitsDB(CRef<CSearchResultSet> results, BlastArrowBuilders &arrow_builders, CRef<CScope> scope){
  // --- PRE-COMPUTE TOTAL ALIGNMENTS ---
  int64_t total_hsp_count = 0;
  for (const auto& res : *results) {
    if (res && res->HasAlignments()) {
      const auto& align_set = res->GetSeqAlign();
      if (align_set && !align_set->IsEmpty()) {
        total_hsp_count += align_set->Get().size();
      }
    }
  }
  
  if (total_hsp_count == 0) //continue; // Skip to next batch, nothing to build
    return empty_rb;
  
  arrow_builders.Reserve(total_hsp_count);
  
  // --- PRE-ALLOCATE REUSABLE STRINGS ---
  std::string q_full, s_full, strand_str, frames, q_aligned, s_aligned;
  q_full.reserve(4096); s_full.reserve(4096);
  q_aligned.reserve(1024); s_aligned.reserve(1024);
  
  int num_rows = 0;
  ncbi::sequence::CDeflineGenerator defline_gen;
  for (const auto& res : *results) {
    if (!res->HasAlignments()) continue;
    RcppThread::checkUserInterrupt();
    
    const CSeq_align_set& align_set = *res->GetSeqAlign();
    std::string qseq_id;
    int64_t parent_list_size = align_set.Size();
    
    // Get Query Label
    CBioseq_Handle q_bh = scope->GetBioseqHandle(*res->GetSeqId());
    if (q_bh) {
      // 1. Pull the original FASTA line from the Bioseq's title/description
      std::string full_header = defline_gen.GenerateDefline(q_bh);
      // 2. FASTA headers often have a structure like ">ID Description..."
      // If you ONLY want the "ID" part (the first word), split it at the first space.
      size_t space_pos = full_header.find_first_of(" \t");
      if (space_pos != std::string::npos) {
        qseq_id = full_header.substr(0, space_pos);
      } else {
        qseq_id = full_header;
      }
      // 3. Fallback just in case the title is truly empty
      if (qseq_id.empty()) {
        qseq_id = res->GetSeqId()->GetSeqIdString(true);
      }
    } else {
      qseq_id = res->GetSeqId()->GetSeqIdString(true);
    }
    
    q_full.clear(); s_full.clear(); // Reset for this result
  
    if (save_sequences) {
      ncbi::objects::CBioseq_Handle bh = scope->GetBioseqHandle(*res->GetSeqId());
      if (bh) {
        ncbi::objects::CSeqVector v = q_bh.GetSeqVector();
        v.SetIupacCoding();
        v.GetSeqData(0, v.size(), q_full);
      }
    }
    
    for (const auto& it : align_set.Get()) {
      RcppThread::checkUserInterrupt();
      
      q_aligned.clear(); s_aligned.clear(); // Clean buffers
      
      std::string sseq_id;
      CBioseq_Handle s_bh = scope->GetBioseqHandle(it->GetSeq_id(1));
      if (s_bh) {
        // 1. Pull the original FASTA line from the Bioseq's title/description
        std::string full_header = defline_gen.GenerateDefline(s_bh);
        // 2. FASTA headers often have a structure like ">ID Description..."
        // If you ONLY want the "ID" part (the first word), split it at the first space.
        size_t space_pos = full_header.find_first_of(" \t");
        if (space_pos != std::string::npos) {
          sseq_id = full_header.substr(0, space_pos);
        } else {
          sseq_id = full_header;
        }
        // 3. Fallback just in case the title is truly empty
        if (sseq_id.empty()) {
          sseq_id = it->GetSeq_id(1).GetSeqIdString(true);
        }
      } else {
        sseq_id = it->GetSeq_id(1).GetSeqIdString(true);
      }
      
      // STRINGS MUST USE Append() - Primitives use UnsafeAppend()
      static_cast<void>(arrow_builders.qseqid.Append(qseq_id));
      static_cast<void>(arrow_builders.sseqid.Append(sseq_id));
      
      TSeqPos qlen = scope->GetBioseqHandle(it->GetSeq_id(0)).GetBioseqLength();
      TSeqPos slen = scope->GetBioseqHandle(it->GetSeq_id(1)).GetBioseqLength();
      static_cast<void>(arrow_builders.qlen.UnsafeAppend(qlen));
      static_cast<void>(arrow_builders.slen.UnsafeAppend(slen));
      
      ncbi::objects::ENa_strand q_strand = it->GetSeqStrand(0);
      ncbi::objects::ENa_strand s_strand = it->GetSeqStrand(1);
      
      char q_strand_char = (q_strand == ncbi::objects::eNa_strand_plus) ? '+' : (q_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
      char s_strand_char = (s_strand == ncbi::objects::eNa_strand_plus) ? '+' : (s_strand == ncbi::objects::eNa_strand_minus) ? '-' : '*';
      strand_str = std::string(1, q_strand_char) + "/" + s_strand_char;
      
      if (it->GetSegs().IsDenseg()) {
        const auto& dseg = it->GetSegs().GetDenseg();
        if(save_sequences && dseg.CanGetIds() ){
          if (dseg.GetIds().size() > 0) 
            GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[0]), scope, q_full);
          if (dseg.GetIds().size() > 1) 
            GetFullSequenceString(const_cast<CRef<CSeq_id>&>(dseg.GetIds()[1]), scope, s_full);
        }
        if (save_hsp_sequences) {
          // GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned);
          GetHSPSequencesGracefully(dseg, scope, q_aligned, s_aligned);
        }
      }
      
      double score=0, n_splices=0, num_ident=0, aln_len=0, gaps=0, mismatches=0, positive=0, negative_count=0;
      double bits=0, evalue=0, blast_score=0, pident=0, aln_len01=0, pident_gap=0, qcovhsp=0;
      double sum_evalue=0, product_coverage=0, overall_identity=0, matches=0, high_quality_percent_coverage=0, exon_identity=0, consensus_splices=0, comp_adj_method=0;
      
      if(!it->GetNamedScore(CSeq_align::EScoreType::eScore_AlignLength, aln_len)) aln_len = it->GetAlignLength(true);
      it->GetNamedScore(CSeq_align::EScoreType::eScore_BitScore, bits); 
      it->GetNamedScore(CSeq_align::EScoreType::eScore_Blast, blast_score);
      bool hasid = it->GetNamedScore(CSeq_align::EScoreType::eScore_IdentityCount, num_ident);
      
      if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped, pident) && hasid) pident = 100.0 * num_ident / it->GetAlignLength(false); 
      if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_PercentIdentity, pident_gap) && hasid) pident_gap = 100.0 * num_ident / it->GetAlignLength(true); 
      if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_GapCount, gaps)){
        // gaps = it->GetTotalGapCount(-1);
        // Proactively check the segment type to prevent NCBI_THROW from logging errors
        const auto& segs = it->GetSegs();
        
        if (segs.IsDenseg() || segs.IsPacked()) {
          // Standard alignments (blastn, blastp, blastx) support this natively
          gaps = it->GetTotalGapCount(-1);
        } else {
          // tblastx uses Std-seg or Dendiag which represent ungapped blocks.
          // Since there are no gaps in these specific alignment formats, it is exactly 0.
          gaps = 0;
        }
      }
      
      it->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue); 
      if (!it->GetNamedScore(CSeq_align::EScoreType::eScore_MismatchCount, mismatches)) mismatches = it->GetAlignLength(true) - num_ident - gaps;
      qcovhsp = q_full.length() > 0 ? (static_cast<double>(it->GetAlignLength(false)) / static_cast<double>(q_full.length())) : 0.0;
      
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
      int qstart = it->GetSeqStart(0);
      int qend = it->GetSeqStop(0);
      int sstart = it->GetSeqStart(1);
      int send = it->GetSeqStop(1);
      
      frames = std::to_string(GetFrame(qstart, aln_len, q_strand)) + "/" + std::to_string(GetFrame(sstart, aln_len, s_strand));
      
      // MIXED APPENDS: Append for Strings, UnsafeAppend for primitives
      static_cast<void>(arrow_builders.qhsp.Append(q_aligned));
      static_cast<void>(arrow_builders.shsp.Append(s_aligned));
      static_cast<void>(arrow_builders.frames.Append(frames));
      static_cast<void>(arrow_builders.strand.Append(strand_str));
      static_cast<void>(arrow_builders.qseq.Append(save_sequences ? q_full : ""));
      static_cast<void>(arrow_builders.sseq.Append(save_sequences ? s_full : ""));
      
      static_cast<void>(arrow_builders.qstart.UnsafeAppend(qstart));
      static_cast<void>(arrow_builders.qend.UnsafeAppend(qend));
      static_cast<void>(arrow_builders.sstart.UnsafeAppend(sstart));
      static_cast<void>(arrow_builders.send.UnsafeAppend(send));
      static_cast<void>(arrow_builders.pident.UnsafeAppend(pident));
      static_cast<void>(arrow_builders.evalue.UnsafeAppend(evalue));
      static_cast<void>(arrow_builders.length.UnsafeAppend(aln_len));
      static_cast<void>(arrow_builders.aln_len01.UnsafeAppend(aln_len01));
      static_cast<void>(arrow_builders.bitscore.UnsafeAppend(bits));
      static_cast<void>(arrow_builders.score.UnsafeAppend(score));
      static_cast<void>(arrow_builders.qcovhsp.UnsafeAppend(qcovhsp));
      static_cast<void>(arrow_builders.blast_score.UnsafeAppend(blast_score));
      static_cast<void>(arrow_builders.pident_gap.UnsafeAppend(pident_gap));
      static_cast<void>(arrow_builders.gaps.UnsafeAppend(gaps));
      static_cast<void>(arrow_builders.nident.UnsafeAppend(num_ident));
      static_cast<void>(arrow_builders.mismatch.UnsafeAppend(mismatches));
      static_cast<void>(arrow_builders.positive.UnsafeAppend(positive));
      static_cast<void>(arrow_builders.n_splices.UnsafeAppend(n_splices));
      static_cast<void>(arrow_builders.hsp_cnt.UnsafeAppend(num_rows + 1));
      static_cast<void>(arrow_builders.sum_evalue.UnsafeAppend(sum_evalue));
      static_cast<void>(arrow_builders.product_coverage.UnsafeAppend(product_coverage));
      static_cast<void>(arrow_builders.overall_identity.UnsafeAppend(overall_identity));
      static_cast<void>(arrow_builders.negative_count.UnsafeAppend(negative_count));
      static_cast<void>(arrow_builders.matches.UnsafeAppend(matches));
      static_cast<void>(arrow_builders.high_quality_percent_coverage.UnsafeAppend(high_quality_percent_coverage));
      static_cast<void>(arrow_builders.exon_identity.UnsafeAppend(exon_identity));
      static_cast<void>(arrow_builders.consensus_splices.UnsafeAppend(consensus_splices));
      static_cast<void>(arrow_builders.comp_adj_method.UnsafeAppend(comp_adj_method));
      static_cast<void>(arrow_builders.num_alignments.UnsafeAppend(parent_list_size));
      static_cast<void>(arrow_builders.hsp_offset.UnsafeAppend(1));
      
      num_rows++;
    }
  }
  
  // --- GENERATE RECORD BATCH ---
  if (num_rows > 0) {
    std::shared_ptr<arrow::Array> qhsp_array, shsp_array, frames_array, pident_array, pident_gap_array, evalue_array, length_array, qstart_array, qend_array, sstart_array, send_array, aln_len01_array, bitscore_array, score_array, qcovhsp_array, blast_score_array, gaps_array, nident_array, mismatch_array, positive_array, n_splices_array, hsp_cnt_array, sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array, matches_array, high_quality_percent_coverage_array, exon_identity_array, consensus_splices_array, comp_adj_method_array;
    
    static_cast<void>(arrow_builders.qhsp.Finish(&qhsp_array));
    static_cast<void>(arrow_builders.shsp.Finish(&shsp_array));
    static_cast<void>(arrow_builders.frames.Finish(&frames_array));
    static_cast<void>(arrow_builders.pident.Finish(&pident_array));
    static_cast<void>(arrow_builders.pident_gap.Finish(&pident_gap_array));
    static_cast<void>(arrow_builders.evalue.Finish(&evalue_array));
    static_cast<void>(arrow_builders.length.Finish(&length_array));
    static_cast<void>(arrow_builders.qstart.Finish(&qstart_array));
    static_cast<void>(arrow_builders.qend.Finish(&qend_array));
    static_cast<void>(arrow_builders.sstart.Finish(&sstart_array));
    static_cast<void>(arrow_builders.send.Finish(&send_array));
    static_cast<void>(arrow_builders.aln_len01.Finish(&aln_len01_array));
    static_cast<void>(arrow_builders.bitscore.Finish(&bitscore_array));
    static_cast<void>(arrow_builders.score.Finish(&score_array));
    static_cast<void>(arrow_builders.qcovhsp.Finish(&qcovhsp_array));
    static_cast<void>(arrow_builders.blast_score.Finish(&blast_score_array));
    static_cast<void>(arrow_builders.gaps.Finish(&gaps_array));
    static_cast<void>(arrow_builders.nident.Finish(&nident_array));
    static_cast<void>(arrow_builders.mismatch.Finish(&mismatch_array));
    static_cast<void>(arrow_builders.positive.Finish(&positive_array));
    static_cast<void>(arrow_builders.n_splices.Finish(&n_splices_array));
    static_cast<void>(arrow_builders.hsp_cnt.Finish(&hsp_cnt_array));
    static_cast<void>(arrow_builders.sum_evalue.Finish(&sum_evalue_array));
    static_cast<void>(arrow_builders.product_coverage.Finish(&product_coverage_array));
    static_cast<void>(arrow_builders.overall_identity.Finish(&overall_identity_array));
    static_cast<void>(arrow_builders.negative_count.Finish(&negative_count_array));
    static_cast<void>(arrow_builders.matches.Finish(&matches_array));
    static_cast<void>(arrow_builders.high_quality_percent_coverage.Finish(&high_quality_percent_coverage_array));
    static_cast<void>(arrow_builders.exon_identity.Finish(&exon_identity_array));
    static_cast<void>(arrow_builders.consensus_splices.Finish(&consensus_splices_array));
    static_cast<void>(arrow_builders.comp_adj_method.Finish(&comp_adj_method_array));
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> aln_struct_array = arrow::StructArray::Make({
      qhsp_array, shsp_array, pident_array, pident_gap_array, frames_array, evalue_array, length_array, aln_len01_array, qstart_array, qend_array, sstart_array, send_array, bitscore_array, score_array, qcovhsp_array, blast_score_array, gaps_array, nident_array, mismatch_array, positive_array, n_splices_array, hsp_cnt_array, sum_evalue_array, product_coverage_array, overall_identity_array, negative_count_array, matches_array, high_quality_percent_coverage_array, exon_identity_array, consensus_splices_array, comp_adj_method_array},
      {"qhsp", "shsp", "pident", "pident_gap", "frames", "evalue", "length", "length01", "qstart", "qend", "sstart", "send", "bitscore", "score", "qcovhsp", "blast_score", "gaps", "nident", "mismatch", "positive", "n_splices", "hsp_num", "sum_evalue", "product_coverage", "overall_identity", "negative_count", "matches", "high_quality_percent_coverage", "exon_identity", "consensus_splices", "comp_adj_method"});
    
    if (!aln_struct_array.ok()) {
      Rcpp::Rcerr << "[BLAST_dbs()] 1. Failed to build StructArray: " << aln_struct_array.status().ToString() << std::endl << std::flush;
      quickblast_running.store(false); 
      // continue;
      return empty_rb;
    }
    
    std::shared_ptr<arrow::StructArray> aln_struct_array_ = aln_struct_array.ValueOrDie();
    
    std::shared_ptr<arrow::Array> qseqid_array, sseqid_array, qseq_array, sseq_array, qlen_array, slen_array, strand_array, num_alignment_array;
    static_cast<void>(arrow_builders.qseqid.Finish(&qseqid_array));
    static_cast<void>(arrow_builders.sseqid.Finish(&sseqid_array));
    static_cast<void>(arrow_builders.qseq.Finish(&qseq_array));
    static_cast<void>(arrow_builders.sseq.Finish(&sseq_array));
    static_cast<void>(arrow_builders.qlen.Finish(&qlen_array));
    static_cast<void>(arrow_builders.slen.Finish(&slen_array));
    static_cast<void>(arrow_builders.strand.Finish(&strand_array));
    static_cast<void>(arrow_builders.num_alignments.Finish(&num_alignment_array));
    
    std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make(
    {num_alignment_array, seqids_struct_array, seqs_struct_array, strand_array, lengths_struct_array},
    {"num_alignments", "seqids", "seqs", "strands", "lengths"});
    
    if(!seq_info_array.ok()){
      Rcpp::Rcerr << "[BLAST_dbs()] 2. Failed to build StructArray: " << seq_info_array.status().ToString() << std::endl << std::flush;
      quickblast_running.store(false); 
      // continue;
      return empty_rb;
    }
    
    std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
    
    std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(), num_rows, {seq_info_array_, aln_struct_array_});
    
    if(alignment_rb && alignment_rb->ValidateFull().ok()){
      // local_ret->emplace_back(alignment_rb);
      return alignment_rb;
    }else{
      return empty_rb;
    }
    
    // Reset Builders for the next batch (retains memory capacity for speed)
    // arrow_builders.Reset();
    // hsp_offset_builder.Reset(); length_builder.Reset(); mismatch_builder.Reset(); gapopen_builder.Reset(); qstart_builder.Reset(); qend_builder.Reset(); sstart_builder.Reset(); send_builder.Reset(); gaps_builder.Reset(); nident_builder.Reset(); positive_builder.Reset(); n_splices_builder.Reset(); hsp_cnt_builder.Reset(); negative_count_builder.Reset(); pident_builder.Reset(); pident_gap_builder.Reset(); evalue_builder.Reset(); bitscore_builder.Reset(); score_builder.Reset(); qcovhsp_builder.Reset(); blast_score_builder.Reset(); aln_len01_builder.Reset(); sum_evalue_builder.Reset(); product_coverage_builder.Reset(); overall_identity_builder.Reset(); matches_builder.Reset(); high_quality_percent_coverage_builder.Reset(); exon_identity_builder.Reset(); consensus_splices_builder.Reset(); comp_adj_method_builder.Reset(); frames_builder.Reset(); strand_builder.Reset(); qseqid_builder.Reset(); sseqid_builder.Reset(); qseq_builder.Reset(); sseq_builder.Reset(); qhsp_builder.Reset(); shsp_builder.Reset(); qlen_builder.Reset(); slen_builder.Reset(); num_alignments_builder.Reset();
  }
  
  return empty_rb;
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::BLAST_dbs(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, unsigned int batch_size, const bool& enable_chunking, unsigned int chunk_size, unsigned int overlap, bool verbose) //const bool show_progress
{
  auto blast_interrupt_ctx = std::make_shared<InterruptContext>();
  safe_jthread interrupt_check_thread(std::thread([this, num_threads, blast_interrupt_ctx](){
    try {
      const unsigned int thread_wait = num_threads > 1 ? 50 : 0;
      while (!blast_interrupt_ctx->stop.load() && quickblast_running.load()) {
        if (RcppThread::isInterrupted()) {  // safe, non-throwing
          blast_interrupt_ctx->stop.store(true);
          break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(thread_wait));
      }
      if(!quickblast_running.load()){
        blast_interrupt_ctx->stop.store(true);
      }
    } catch(...) {
      blast_interrupt_ctx->stop.store(true);
    }
  }));
  // interrupt_check_thread.detach();
  try{
    
    quickblast_running.store(true);      
    
    std::atomic<bool> s_batches_done{false};
    std::atomic<bool> q_batches_done{false};
    
    // Create per-thread CSeqDB
    CSeqDB::ESeqType seqdbType;
    CSearchDatabase::EMoleculeType seqType;
    CBlastDbDataLoader::EDbType edbType;
    switch(seq_type){
    case QuickBLAST::ESeqType::eNucleotide: {
      seqType = CSearchDatabase::EMoleculeType::eBlastDbIsNucleotide;
      seqdbType = CSeqDB::eNucleotide;
      edbType =  CBlastDbDataLoader::EDbType::eNucleotide;
      break;
    }
    case QuickBLAST::ESeqType::eProtein: {
      seqType = CSearchDatabase::EMoleculeType::eBlastDbIsProtein;
      seqdbType = CSeqDB::eProtein;
      edbType = CBlastDbDataLoader::EDbType::eProtein;
      break;
    }
    }
    
    CRef<CObjectManager> objMgr(CObjectManager::GetInstance());
    if (!objMgr) {
      Rcpp::Rcerr << "BLAST_dbs: CObjectManager::GetInstance() returned NULL." << std::endl << std::flush;
      quickblast_running.store(false); 
      // interrupt_check_thread.join();
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    std::string loader_name = CBlastDbDataLoader::RegisterInObjectManager(
      *objMgr, queryFile, edbType, true, 
      CObjectManager::eDefault, CObjectManager::kPriority_NotSet
    ).GetLoader()->GetName();
    
    CRef<CSeqDB> q_seqdb_(new CSeqDB(queryFile, seqdbType));
    CRef<CSeqDB> s_seqdb_(new CSeqDB(queryFile, seqdbType));
    // CRef<CSearchDatabase> s_serdb_(new CSearchDatabase(subjectFile, seqType));
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    q_seqdb_->SetNumberOfThreads(num_threads, /*force_mt*/ true);
    s_seqdb_->SetNumberOfThreads(num_threads, /*force_mt*/ true);
#else
    q_seqdb_->SetNumberOfThreads(1, /*force_mt*/ false);
    s_seqdb_->SetNumberOfThreads(1, /*force_mt*/ false);
#endif
    
    int num_queries = q_seqdb_->GetNumOIDs();
    
    arrow_wrapper->AddFASTAMetadata("BLAST locality", "CBlastOptions::eLocal");
    arrow_wrapper->AddFASTAMetadata("Input source", "DBs");     
    arrow::Status outfile_sts = arrow_wrapper->CreateOutputStream(outFile, outFormat);
    if (!outfile_sts.ok())
    {
      Rcpp::Rcerr << std::string("[BLAST_dbs()] ERROR : Could not create output file stream : ") << outfile_sts.detail()->ToString() << std::endl << outfile_sts.message() << std::endl << std::flush;
      quickblast_running.store(false); 
      // interrupt_check_thread.join();
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    unsigned int n_threads = num_threads; // 1;
    SetThreadCount(n_threads);
    
    const unsigned int q_seq_count = q_seqdb_->GetNumSeqs();
    const unsigned int s_seq_count = s_seqdb_->GetNumSeqs();
    const unsigned int totalIterations = q_seq_count * s_seq_count;
    
    arrow_wrapper->SetBLASTSeqLimit(batch_size);
    
    if(q_seq_count > n_threads){
      n_threads = int(ceil(n_threads / 2) - 2) <= 0 ? 1 : int(ceil(n_threads / 2) - 2);
    }
    if(n_threads > q_seq_count + s_seq_count)
      n_threads = q_seq_count + s_seq_count;
    SetThreadCount(n_threads);
    
    if(verbose){    
      Rcpp::Rcout << "Num Threads: " << n_threads << std::endl << std::flush; //DEBUG
      // std::cout << "BLAST Sequence Limit: " << blast_sequence_limit << std::endl << std::flush; //DEBUG
      Rcpp::Rcout << "Total Records (Q + S): " << q_seq_count + s_seq_count << " (" << q_seq_count << " + " << s_seq_count << ")"<< std::endl << std::flush; //DEBUG
      if(enable_chunking){
        Rcpp::Rcout << "Chunking Enabled..."<< std::endl << std::flush; //DEBUG
      }
    }
    if(totalIterations <= 0){
      Rcpp::Rcerr << "[BLAST_dbs()] Improperly formatted DB file. No records detected." << std::endl << std::flush;
      quickblast_running.store(false); 
      // interrupt_check_thread.join();
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    if (enable_chunking) {
      if (chunk_size <= 0) {
        Rcpp::Rcerr << "QuickBLAST Error: chunk_size must be greater than 0." << std::endl << std::flush;
        quickblast_running.store(false); 
        return std::make_shared<arrow::RecordBatchVector>();
      }
      if (overlap >= chunk_size) {
        Rcpp::Rcerr << "QuickBLAST Error: overlap must be strictly less than chunk_size." << std::endl << std::flush;
        quickblast_running.store(false); 
        return std::make_shared<arrow::RecordBatchVector>();
      }
    }
    
    // Calculate safe step size globally
    const TSeqPos step_size = enable_chunking ? (chunk_size - overlap) : 0;
    
    if(n_threads > 1 && batch_size <= 0)
      batch_size = n_threads + 1;
    else if(batch_size <= 1)
      batch_size = 2;
    
    arrow_wrapper->SetVerbosity(verbose);
    arrow_wrapper->SetBatchSize(batch_size);
    
    std::shared_ptr<arrow::RecordBatchVector> final_ret = std::make_shared<arrow::RecordBatchVector>();
    
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
    
    // // Check if program and sequence type match
    // // Get the expected query type for the BLAST program
    // ncbi::blast::EProgram program_enum = ncbi::blast::ProgramNameToEnum(program);
    // bool expects_protein = ncbi::blast::Blast_QueryIsProtein(program_enum);
    Progress progress_bar(q_seq_count, verbose);
    int chunk_counter = 0;
    int batch_counter = 0;
    
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
    #pragma omp parallel
    #endif
    {
      try {
        BlastArrowBuilders arrow_builders;
        
        // Thread-local storage for hits (Safe to append without locks!)
        std::shared_ptr<arrow::RecordBatchVector> local_ret = std::make_shared<arrow::RecordBatchVector>();
        
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    #pragma omp for schedule(dynamic)
    #endif
        for (int i = 0; i < num_queries; i += batch_size) {
          if(!quickblast_running.load()) continue;
          
          // 1. FRESH SCOPE PER BATCH: Prevents memory caching leaks
          CRef<CScope> batch_scope(new CScope(*objMgr));
          batch_scope->AddDataLoader(loader_name);
          
          // 2. FRESH QUERY DB PER BATCH
          CRef<CSeqDB> lcl_q_seqdb(new CSeqDB(queryFile, seqdbType));
          lcl_q_seqdb->SetNumberOfThreads(1, false);
          
          // 3. FRESH SUBJECT DB PER BATCH: Stops the Double-Free SIGSEGV!
          // Make sure to replace `subjectFile` and `seqdbType` with your actual variables
          CRef<CSearchDatabase> lcl_s_serdb(new CSearchDatabase(subjectFile, seqType)); 
          
          int current_batch_end = std::min<int>(i + batch_size, num_queries);
          
          try {
            for (int j = i; j < current_batch_end; ++j) {
              RcppThread::checkUserInterrupt();
              CRef<CSeq_id> id = lcl_q_seqdb->GetSeqIDs(j).front();
              TSeqPos seq_len = lcl_q_seqdb->GetSeqLength(j);
              
              if (enable_chunking && seq_len > chunk_size) {
                
                for (TSeqPos start = 0; start < seq_len; start += step_size) { //(chunk_size - overlap)
                  TSeqPos end = std::min(start + chunk_size, seq_len) - 1;
                  
                  CRef<CSeq_interval> seq_int(new CSeq_interval());
                  seq_int->SetId().Assign(*id);
                  seq_int->SetFrom(start);
                  seq_int->SetTo(end);
                  
                  CRef<CSeq_loc> chunk_loc(new CSeq_loc());
                  chunk_loc->SetInt(*seq_int);
                  
                  // Bind query to our fresh batch scope
                  CRef<CBlastQueryVector> single_chunk_query(new CBlastQueryVector());
                  single_chunk_query->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*chunk_loc, *batch_scope)));
                  
                  CRef<IQueryFactory> chunk_factory(new CObjMgr_QueryFactory(*single_chunk_query));
                  
                  // Use the perfectly safe, thread-local Subject DB
                  CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
                  CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, false);
                  
                  CLocalBlast lcl_blaster(chunk_factory, lcl_blast_opts, s_db_adapter);
                  lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
                  CRef<CSearchResultSet> results = lcl_blaster.Run();
                  // Pass batch_scope to ExtractHitsDB
                  std::shared_ptr<arrow::RecordBatch> arb = ExtractHitsDB(results, arrow_builders, batch_scope);
                  
                  if(arb && arb->ValidateFull().ok()) {
                    // NO CRITICAL SECTION NEEDED HERE! local_ret is already thread-local.
                    local_ret->emplace_back(arb);
                    
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
    #pragma omp atomic update
    #endif
                    chunk_counter++;
                  }
                }
                arrow_builders.Reset();
              } else { //without chunking
                CRef<CSeq_loc> loc(new CSeq_loc());
                loc->SetWhole(*id);
                CRef<CBlastQueryVector> queries(new CBlastQueryVector());
                queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*loc, *batch_scope)));
                CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
                CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(*queries));
                CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, /*verbose*/ false);
                // 5. RUN ENGINE
                CLocalBlast lcl_blaster(query_factory, lcl_blast_opts, s_db_adapter);
                lcl_blaster.SetBatchNumber(i);
                lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
                CRef<CSearchResultSet> results;
                results = lcl_blaster.Run();
                if(verbose)
                  Rcpp::Rcout << "Batch Hits: " << ((results.NotEmpty() && results->GetNumResults() > 0 && (*results)[0].HasAlignments()) ? (*results)[0].GetSeqAlign()->Get().size() : 0) << std::endl;
                std::shared_ptr<arrow::RecordBatch> arb = ExtractHitsDB(results, arrow_builders, batch_scope);
                if(arb && arb->ValidateFull().ok()){
                  local_ret->emplace_back(arb);
                }
                
#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp atomic update
#endif
                batch_counter++;
              }
              arrow_wrapper->AddProcRecordCount();
            }       
          } catch (const ncbi::CException& e) {
            quickblast_running.store(false);
            Rcpp::Rcerr << "[BLAST_dbs()] 1. Execution error: " << e.GetMsg() << std::endl << std::flush;
            continue;
          }catch (const std::exception& e) {
            quickblast_running.store(false);
            Rcpp::Rcerr << "[BLAST_dbs()] 1. C++ Exception: " << e.what() << std::endl << std::flush;
            continue;
          }
          
          for (int pbi = 0; pbi < current_batch_end; pbi++) {  
    #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
    #pragma omp critical (pgbr_incr)
    #endif
    {
      progress_bar.increment();
    }
          }
          
          // batch_scope gracefully dies here, freeing all cached sequences!
          batch_scope->ResetHistory(CScope::EActionIfLocked::eRemoveIfLocked);
        } // end of batch loop
        
        arrow_builders.Reset();
        // Aggregate local batches safely into the global return vector
        if (!local_ret->empty()) {
    #pragma omp critical(final_ret_insert)
    {
      final_ret->insert(final_ret->end(), std::make_move_iterator(local_ret->begin()), std::make_move_iterator(local_ret->end()));
    }
        }
      }catch (const ncbi::CException& e) {
        Rcpp::Rcerr << "[BLAST_dbs()] 2. NCBI CException: " << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
      }catch (const std::runtime_error& e) {
        Rcpp::Rcerr << "[BLAST_dbs()] 2. C++ Runtime error: " << e.what() << std::endl << std::flush;
      }catch (const std::exception& e) {
        Rcpp::Rcerr << "[BLAST_dbs()] 2. C++ Exception: " << e.what() << std::endl << std::flush;
      }catch (...) {
        Rcpp::Rcerr << "[BLAST_dbs()] 2. Unknown error" << std::endl << std::flush;
      }
    } //end of omp parallel

// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
// #pragma omp parallel
// #endif
// {
//   try{
//     
//     CRef<CScope> thread_scope(new CScope(*objMgr));
//     thread_scope->AddDataLoader(loader_name);
//     BlastArrowBuilders arrow_builders;
//     
//     std::shared_ptr<arrow::RecordBatchVector> local_ret = std::make_shared<arrow::RecordBatchVector>();
//     
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
// #pragma omp for schedule(dynamic)
// #endif
//     for (int i = 0; i < num_queries; i += batch_size) {
//       if(!quickblast_running.load()) continue;
//       // CRef<CScope> scope(new CScope(*objMgr));
//       // scope->AddDataLoader(loader_name);
//       // // scope->AddDefaults();
//       
//       CRef<CSeqDB> lcl_q_seqdb(new CSeqDB(queryFile, seqdbType));
//       lcl_q_seqdb->SetNumberOfThreads(1,false);
//       int current_batch_end = std::min<int>(i + batch_size, num_queries);
//       
//       try {
//         for (int j = i; j < current_batch_end; ++j) {
//           RcppThread::checkUserInterrupt();
//           CRef<CSeq_id> id = lcl_q_seqdb->GetSeqIDs(j).front();
//           
//           if (enable_chunking) {
//             Rcpp::Rcout << "HERE1" << std::endl << std::flush; //DEBUG
//             TSeqPos seq_len = lcl_q_seqdb->GetSeqLength(j);
//             for (TSeqPos start = 0; start < seq_len; start += (chunk_size - overlap)) {
//               TSeqPos end = std::min(start + chunk_size, seq_len) - 1;
//               
//               CRef<CSeq_interval> seq_int(new CSeq_interval());
//               seq_int->SetId().Assign(*id);
//               seq_int->SetFrom(start);
//               seq_int->SetTo(end);
//               
//               CRef<CSeq_loc> chunk_loc(new CSeq_loc());
//               chunk_loc->SetInt(*seq_int);
//               
//               // 1. FRESH QUERY VECTOR (Holds ONLY this 50kb chunk)
//               CRef<CBlastQueryVector> single_chunk_query(new CBlastQueryVector());
//               single_chunk_query->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*chunk_loc, *thread_scope)));
//               
//               // 2. FRESH ENGINE PIPELINE for this chunk
//               CRef<IQueryFactory> chunk_factory(new CObjMgr_QueryFactory(*single_chunk_query));
//               
//               CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*s_serdb_));
//               CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, false);
//               CLocalBlast lcl_blaster(chunk_factory, lcl_blast_opts, s_db_adapter);
//               lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
//               CRef<CSearchResultSet> results;
//               Rcpp::Rcout << "HERE1.2" << std::endl << std::flush; //DEBUG  
//               // RUN THE SEARCH (Now it only processes 50kb, perfectly safe!)
//               results = lcl_blaster.Run();
//               if(verbose)
//                 Rcpp::Rcout << "Chunk Hits: " << ((results.NotEmpty() && results->GetNumResults() > 0 && (*results)[0].HasAlignments()) ? (*results)[0].GetSeqAlign()->Get().size() : 0) << std::endl;
//               std::shared_ptr<arrow::RecordBatch> arb = ExtractHitsDB(results, arrow_builders, thread_scope);
// if(arb && arb->ValidateFull().ok()){
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
// #pragma omp critical (ret_append)
// #endif
// {
//   local_ret->emplace_back(arb);
// }
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
// #pragma omp atomic update
//   chunk_counter++;
// #else
//   chunk_counter++;
// #endif
//   Rcpp::Rcout << "HERE3" << chunk_counter << std::endl << std::flush; //DEBUG
// }
//             }
//             // Reset Builders for the next batch (retains memory capacity for speed)
//             arrow_builders.Reset();
//           }else{
//             CRef<CSeq_loc> loc(new CSeq_loc());
//             loc->SetWhole(*id);
//             CRef<CBlastQueryVector> queries(new CBlastQueryVector());
//             queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*loc, *thread_scope)));
//             CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*s_serdb_));
//             CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(*queries));
//             CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, /*verbose*/ false);
//             // 5. RUN ENGINE
//             CLocalBlast lcl_blaster(query_factory, lcl_blast_opts, s_db_adapter);
//             lcl_blaster.SetBatchNumber(i);
//             lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
//             CRef<CSearchResultSet> results;
//             Rcpp::Rcout << "HERE3" << std::endl << std::flush; //DEBUG
//             results = lcl_blaster.Run();
//             if(verbose)
//               Rcpp::Rcout << "Batch Hits: " << ((results.NotEmpty() && results->GetNumResults() > 0 && (*results)[0].HasAlignments()) ? (*results)[0].GetSeqAlign()->Get().size() : 0) << std::endl;
// Rcpp::Rcout << "HERE4" << std::endl << std::flush; //DEBUG
// std::shared_ptr<arrow::RecordBatch> arb = ExtractHitsDB(results, arrow_builders, thread_scope);
// if(arb && arb->ValidateFull().ok()){
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
// #pragma omp critical (ret_append)
// #endif
// {
//   local_ret->emplace_back(arb);
// }
// }
//           }
//           arrow_wrapper->AddProcRecordCount();
//         }       
//         
//         // CRef<CSeq_loc> loc(new CSeq_loc());
//         // loc->SetWhole(*id);
//         // CRef<CBlastQueryVector> queries(new CBlastQueryVector());
//         // queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*loc, *scope)));
//         // // Per-thread OPTIONS
//         // CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, false);
//         // // Per-thread SUBJECT DATABASE
//         // CRef<CSearchDatabase> lcl_s_serdb(new CSearchDatabase(subjectFile, seqType));
//         // CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
//         // 
//         // CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(*queries));
//         
//         //     if(!enable_chunking){          
//         //         CRef<CSeq_loc> loc(new CSeq_loc());
//         //         loc->SetWhole(*id);
//         //         queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*loc, *scope)));
//         //         arrow_wrapper->AddProcRecordCount();
//         //         
//         //         CRef<CSearchDatabase> lcl_s_serdb(new CSearchDatabase(subjectFile, seqType));
//         //         CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
//         //         CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(*queries));
//         //         CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, /*verbose*/ false);
//         //         // 5. RUN ENGINE
//         //         CLocalBlast lcl_blaster(query_factory, lcl_blast_opts, s_db_adapter);
//         //         lcl_blaster.SetBatchNumber(i);
//         //         lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
//         //         Rcpp::Rcout << "HERE3" << std::endl << std::flush; //DEBUG
//         //         CRef<CSearchResultSet> results = lcl_blaster.Run();
//         //         Rcpp::Rcout << "HERE4" << std::endl << std::flush; //DEBUG
//         //         std::shared_ptr<arrow::RecordBatch> arb = ExtractHitsDB(results, arrow_builders, scope);
//         //         if(arb && arb->ValidateFull().ok()){
//         // #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
//         // #pragma omp critical (ret_append)
//         // #endif
//         // {
//         //   local_ret->emplace_back(arb);
//         // }
//         //         }
//         //     } //if(!enable_chunking)
//         
//       }catch (const ncbi::CException& e) {
//         quickblast_running.store(false);
//         Rcpp::Rcerr << "[BLAST_dbs()] Execution error: " << e.GetMsg() << std::endl << std::flush;
//         continue;
//       }catch (const std::exception& e) {
//         quickblast_running.store(false);
//         Rcpp::Rcerr << "[BLAST_dbs()] C++ Exception: " << e.what() << std::endl << std::flush;
//         continue;
//       }
//       
//       Rcpp::Rcout << "HERE2" << std::endl << std::flush; //DEBUG
//       for (int pbi = 0; pbi < current_batch_end; pbi++) {  
// #if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)    
// #pragma omp critical (pgbr_incr)
// #endif
// {
//   progress_bar.increment();
// }
//       }
//       
//       // CRef<CSearchDatabase> lcl_q_serdb(new CSearchDatabase(queryFile, seqType));
//       // CRef<CBlastQueryVector> queries(new CBlastQueryVector());
//       // int current_batch_end = std::min<int>(i + batch_size, num_queries);
//       // 
//       // for (int j = i; j < current_batch_end; ++j) {
//       //   RcppThread::checkUserInterrupt();
//       //   CRef<CSeq_id> id = lcl_q_serdb->GetSeqIDs(j).front();
//       //   
//       //   CRef<CSeq_loc> loc(new CSeq_loc());
//       //   loc->SetWhole(*id);
//       //   queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*loc, *scope)));
//       //   arrow_wrapper->AddProcRecordCount();
//       // }
//       // CRef<CSearchDatabase> lcl_s_serdb(new CSearchDatabase(subjectFile, seqType));
//       // CRef<CLocalDbAdapter> s_db_adapter(new CLocalDbAdapter(*lcl_s_serdb));
//       // CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(*queries));
//       // CRef<CSearchResultSet> results;
//       // try {
//       //   CRef<CBlastOptionsHandle> lcl_blast_opts = MakeQuickBLASTOptions(program, blast_options, CBlastOptions::eLocal, /*verbose*/ false);
//       //   CLocalBlast lcl_blaster(query_factory, lcl_blast_opts, s_db_adapter);
//       //   lcl_blaster.SetBatchNumber(i);
//       //   lcl_blaster.SetInterruptCallback(&BlastInterruptFn, static_cast<void*>(blast_interrupt_ctx.get()));
//       //   results = lcl_blaster.Run();
//       // } catch (const ncbi::CException& e) {
//       //   Rcpp::Rcerr << "[BLAST_dbs()] Execution error: " << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
//       //   quickblast_running.store(false);
//       //   continue;
//       // }
//       
//       // scope->ResetHistory(CScope::EActionIfLocked::eRemoveIfLocked);
//     } // end of batch loop
//     // Reset Builders for the next batch (retains memory capacity for speed)
//     arrow_builders.Reset();
//     // Aggregate local batches safely
//     if (!local_ret->empty()) {
// #pragma omp critical(final_ret_insert)
// {
//   final_ret->insert(final_ret->end(), std::make_move_iterator(local_ret->begin()), std::make_move_iterator(local_ret->end()));
// }
//     }
//   }catch (const ncbi::CException& e) {
//     Rcpp::Rcerr << "[BLAST_dbs()] 1. NCBI CException: " << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
//   }catch (const std::runtime_error& e) {
//     Rcpp::Rcerr << "[BLAST_dbs()] 1. C++ Runtime error: " << e.what() << std::endl << std::flush;
//   }catch (const std::exception& e) {
//     Rcpp::Rcerr << "[BLAST_dbs()] 1. C++ Exception: " << e.what() << std::endl << std::flush;
//   }catch (...) {
//     Rcpp::Rcerr << "[BLAST_dbs()] 1. Unknown error" << std::endl << std::flush;
//   }
// } //end of omp parallel

#if defined(_OPENMP) && !defined(WIN32) && !defined(MINGW32)
#pragma omp barrier
#endif


// std::cout << "Final Batch Size: " << arrow_wrapper->GetBatchSize() << std::endl << std::flush;  //DEBUG

static_cast<void>(arrow_wrapper->FinishOutputStream());

if(verbose && enable_chunking)
  Rcpp::Rcout << "Processed Chunks:" << chunk_counter << std::endl << std::flush; //DEBUG
if(verbose && !enable_chunking)
  Rcpp::Rcout << "Processed Batches:" << batch_counter << std::endl << std::flush; //DEBUG



if(verbose)
  Rcpp::Rcout << "Total Records Processed: " << arrow_wrapper->GetProcRecordCount() << std::endl << std::flush;  

arrow_wrapper->ResetProcRecordCount();
quickblast_running.store(false); 
// interrupt_check_thread.join();
if (return_values)
{
  return final_ret;
}
else
{
  final_ret->clear();
  final_ret->shrink_to_fit();
  return std::make_shared<arrow::RecordBatchVector>();
}

  }catch (const ncbi::CException &e) {
    // NCBI toolkit exceptions
    Rcpp::Rcerr << std::string("[BLAST_dbs()] 2. NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch (const std::runtime_error &e) {
    std::string msg = "[BLAST_dbs()]: 2. C++ Runtime Error : ";
    Rcpp::Rcerr << msg + e.what() << std::endl << std::flush; 
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[BLAST_dbs()] - 2. Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch (const std::exception &e) {
    Rcpp::Rcerr << std::string("[BLAST_dbs()] - 2. C++ exception : ") + e.what() << std::endl << std::flush;
  }catch (...) {
    Rcpp::Rcerr << "[BLAST_dbs()]: 2. Unknown exception" << std::endl << std::flush;
  }
  static_cast<void>(arrow_wrapper->FinishOutputStream());
  arrow_wrapper->ResetProcRecordCount();
  quickblast_running.store(false);
  // interrupt_check_thread.join();
  return std::make_shared<arrow::RecordBatchVector>();
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::BLAST_seqs(const std::string &query, const std::string &subject, bool verbose)
{
  try {
    // quick user interrupt check if this runs under R
    RcppThread::checkUserInterrupt();
    // ensure arrow_wrapper exists
    if (!this->arrow_wrapper) {
      Rcpp::Rcerr << "BLAST_seqs: arrow_wrapper is null (not initialized)." << std::endl << std::flush;
      return empty_rb;
    }
    arrow_wrapper->SetVerbosity(verbose);
    quickblast_running.store(true); 
    // convert inputs via arrow wrapper and validate
    auto q_type = this->arrow_wrapper->CastToType(query);
    arrow_wrapper->AddRecordCount();
    auto s_type = this->arrow_wrapper->CastToType(subject);
    arrow_wrapper->AddRecordCount();
    if (q_type.header.empty() || q_type.seq.empty()) {
      Rcpp::Rcerr << "BLAST_seqs: query header/sequence is empty." << std::endl << std::flush;
      return empty_rb;
    }
    if (s_type.header.empty() || s_type.seq.empty()) {
      Rcpp::Rcerr << "BLAST_seqs: subject header/sequence is empty." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    
    // get CObjectManager instance and create scope
    CRef<CObjectManager> objMgr(CObjectManager::GetInstance());
    if (!objMgr) {
      Rcpp::Rcerr << "BLAST_seqs: CObjectManager::GetInstance() returned NULL." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    
    CRef<ncbi::CScope> scope(new ncbi::CScope(*objMgr));
    scope->AddDefaults();
    
    auto [ query_seqloc, query_seq_entry ] = CreateSSeqLocFromType(q_type, scope);
    if (!query_seqloc) {
      Rcpp::Rcerr << "BLAST_seqs: CreateSSeqLocFromType(query) returned NULL." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    if(!query_seqloc->seqloc.NotEmpty()){
      Rcpp::Rcerr << "BLAST_seqs: CreateSSeqLocFromType(query) is empty." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    
    auto [ subject_seqloc, subject_seq_entry ] = CreateSSeqLocFromType(s_type, scope);
    if (!subject_seqloc) {
      Rcpp::Rcerr << "BLAST_seqs: CreateSSeqLocFromType(subject) returned NULL." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    if(!subject_seqloc->seqloc.NotEmpty()){
      Rcpp::Rcerr << "BLAST_seqs: CreateSSeqLocFromType(subject) is empty." << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    // create the blaster and run
    CBl2Seq blaster(*query_seqloc, *subject_seqloc, this->GetQuickBLASTOptions());
    TSeqAlignVector alignments;
    
    try{  
      alignments = blaster.Run();
    }catch (const ncbi::CException& e) {
      Rcpp::Rcerr << std::string("[BLAST_seqs()] 1. BLAST Execution error: Check Sequence type mismatch (proteins != nucleotides) : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    
    RcppThread::ProgressBar progress_bar(1, verbose);
    quickblast_running.store(false); 
    return this->ExtractHits(alignments, *query_seqloc, *subject_seqloc, scope, true); 
  }catch (const ncbi::CException &e) {
    Rcpp::Rcerr << std::string("[BLAST_seqs()]: 1. NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch (const std::runtime_error &e) {
    std::string msg = "[BLAST_seqs()]: 1. C++ Runtime Error : ";
    Rcpp::Rcerr << msg + e.what() << std::endl << std::flush; 
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[BLAST_seqs()] - 1. Rcpp Exception : ") + e.what() << std::endl << std::flush; 
  }catch (const std::exception &e) {
    Rcpp::Rcerr << std::string("[BLAST_seqs()] - 1. C++ exception : ") + e.what() << std::endl << std::flush; 
  }catch (...) {
    Rcpp::Rcerr << "[BLAST_seqs()]: 1. Unknown exception" << std::endl << std::flush; 
  }
  quickblast_running.store(false); 
  return empty_rb;
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
 // ' @param blast_sequence_limit (int) Batch Size to BLAST at a time  - Not used for Sequence BLAST
 //' @return Nested List of BLAST Hits
 auto QuickBLAST::Impl::BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, bool verbose)
 {
   
   if(!std::filesystem::exists(query)){
     Rcpp::stop("[BLAST()] query file/folder does not exist.");
   }
   if(!std::filesystem::exists(subject)){
     Rcpp::stop("[BLAST()] subject file/folder does not exist.");
   }
   
   switch (input_type)
   {
   case QuickBLAST::EInputType::eFile:
   {
     unsigned int n_threads = num_threads; 
     int batch_size = 96 * num_threads;
     batch_size = 32 * n_threads;
     return BLAST_files(query, subject, outputFile, outFormat, n_threads, true, batch_size, verbose);
   }
     break;
   case QuickBLAST::EInputType::eSequenceString:
   {
     arrow::RecordBatchVector ret_val;
     ret_val.emplace_back(BLAST_seqs(query, subject, verbose));
     return std::make_shared<arrow::RecordBatchVector>(ret_val);
   }
     break;
   default:
   {
     Rcpp::Rcerr << "input_type must be QuickBLAST::EInputType::eFile (0) OR QuickBLAST::EInputType::eSequenceString (1) !" << std::endl << std::flush;
     return std::make_shared<arrow::RecordBatchVector>();
   }
   break;
   }
 }

CRef<ncbi::blast::CBlastOptionsHandle> QuickBLAST::SetQuickBLASTOptions(const std::string &program_name, const std::string &options, CBlastOptions::EAPILocality locality = CBlastOptions::eLocal, const bool& verbose = true)
{
  return pImpl->SetQuickBLASTOptions(program_name, options, locality, verbose);
}

unsigned int QuickBLAST::GetObjectID()
{
  return pImpl->GetObjectID();
}

void QuickBLAST::SetObjectID(unsigned int id)
{
  pImpl->SetObjectID(id);
}

std::shared_ptr<arrow::Schema> QuickBLAST::GetSchema() { return pImpl->GetSchema(); };

std::string QuickBLAST::GetProgram(){ return pImpl->GetProgram(); };

std::string QuickBLAST::Impl::GetProgram(){
  return program;
}

void QuickBLAST::SetThreadCount(unsigned int num_threads)
{
  pImpl->SetThreadCount(num_threads);
}
unsigned int QuickBLAST::GetThreadCount()
{
  return pImpl->GetThreadCount();
}
int QuickBLAST::GetHitCount()
{
  return pImpl->GetHitCount();
}
void QuickBLAST::AddHitCount(int val)
{
  pImpl->AddHitCount(val);
}
ncbi::blast::CBlastOptionsHandle &QuickBLAST::GetQuickBLASTOptions()
{
  return pImpl->GetQuickBLASTOptions();
}
void QuickBLAST::ResetHitCount() { pImpl->ResetHitCount(); }

SEXP QuickBLAST::Hits2RList(const std::shared_ptr<arrow::RecordBatch> &rb)
{
  return pImpl->Hits2RList(rb);
}

SEXP QuickBLAST::Hits2RList(const arrow::RecordBatchVector &rb_vector)
{
  return pImpl->Hits2RList(rb_vector);
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::Impl::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const TSeqLocVector &sloc, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values) 
{
  try{
    RcppThread::checkUserInterrupt();
    if(alignments.empty()){
      return std::make_shared<arrow::RecordBatchVector>();
    }
    
    std::shared_ptr<arrow::RecordBatchVector> ret_rbv;
    
    if (!qloc.seqloc) {
      Rcpp::Rcerr << "ERROR: ExtractHits: qloc.seqloc is NULL" << std::endl << std::flush;
      ret_rbv->emplace_back(empty_rb);
      return ret_rbv;
    }
    
    std::shared_ptr<arrow::RecordBatchVector> recBth_vec = std::make_shared<arrow::RecordBatchVector>();
    
    for (const auto &s_it : sloc)
    {
      RcppThread::checkUserInterrupt();
      
      std::shared_ptr<arrow::RecordBatch> rb = ExtractHits(alignments, qloc, s_it, scope, return_values); 
      if(return_values)
        if (rb)
        {
          recBth_vec->emplace_back(std::move(rb));
        }
    }
    
    return recBth_vec;
    
  }catch(const ncbi::CException &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[ExtractHits()] - Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }catch(...){
    Rcpp::Rcerr << "[ExtractHits()]: Unknown Exception" << std::endl << std::flush;
  }
  return std::make_shared<arrow::RecordBatchVector>();
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::Impl::ExtractHits(const TSeqAlignVector &alignments, const SSeqLoc &qloc, const SSeqLoc &sloc, ncbi::CRef<ncbi::objects::CScope> scope, const bool &return_values)
{
  try{
    RcppThread::checkUserInterrupt();
    if (alignments.empty()) {
      // return an empty but typed record batch
      return empty_rb;
    }
    
    if (!qloc.seqloc) {
      Rcpp::Rcerr << "ERROR: ExtractHits: qloc.seqloc is NULL" << std::endl << std::flush;
      return empty_rb; 
    }
    if (!sloc.seqloc) {
      Rcpp::Rcerr << "ERROR: ExtractHits: sloc.seqloc is NULL" << std::endl << std::flush;
      return empty_rb; 
    }
    
    std::string strand;
    
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
    
    // Arrow Builders
    arrow::Int64Builder hsp_offset_builder, length_builder, mismatch_builder, gapopen_builder, qstart_builder, qend_builder, sstart_builder, send_builder, gaps_builder, nident_builder, positive_builder, n_splices_builder, hsp_cnt_builder, negative_count_builder;
    arrow::DoubleBuilder pident_builder, pident_gap_builder, evalue_builder, bitscore_builder, score_builder, qcovhsp_builder, blast_score_builder, aln_len01_builder, sum_evalue_builder, product_coverage_builder, overall_identity_builder, matches_builder, high_quality_percent_coverage_builder, exon_identity_builder, consensus_splices_builder, comp_adj_method_builder;
    arrow::StringBuilder frames_builder, strand_builder, qseqid_builder, sseqid_builder;
    arrow::LargeStringBuilder qseq_builder, sseq_builder, qhsp_builder, shsp_builder;
    arrow::Int64Builder qlen_builder, slen_builder, num_alignments_builder;
    
    std::string qseq_id = qloc.seqloc->GetId()->GetSeqIdString(true);
    std::string sseq_id = sloc.seqloc->GetId()->GetSeqIdString(true);
    
    CScoreBuilder scorer;
    // if (effective_search_space > 0.0) scorer.SetEffectiveSearchSpace(effective_search_space);
    
    // Compute batch scores (AddScore has an overload for list)
    // We'll ask for a set of scores in a loop to leverage internal batching
    std::vector<CSeq_align::EScoreType> score_types = {
      CSeq_align::EScoreType::eScore_AlignLength,
      CSeq_align::EScoreType::eScore_BitScore,
      CSeq_align::EScoreType::eScore_Blast,
      CSeq_align::EScoreType::eScore_PercentIdentity_Ungapped,
      CSeq_align::EScoreType::eScore_PercentIdentity_Gapped,
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
    
    // --- Pre-calculate total rows to pre-allocate Arrow builders ---
    int64_t estimated_rows = 0;
    for (const auto &align_set_ref : alignments) {
      if (align_set_ref && align_set_ref->IsSet() && align_set_ref->CanGet()) {
        estimated_rows += align_set_ref->Get().size();
      }
    }
    
    if (estimated_rows == 0) {
      return empty_rb;  //// CORRECT RETURN, NO ALIGNMENTS
    }
    
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
    
    std::string q_full, s_full, q_aligned, s_aligned;
    std::string strand_str, frames;
    
    // Pre-reserve common string capacities to avoid growth during appending
    q_full.reserve(4096); s_full.reserve(4096);
    q_aligned.reserve(1024); s_aligned.reserve(1024);
    
    int num_rows = 0;
    
    for (const auto &align_set_ref : alignments) {
      if (!align_set_ref || align_set_ref.IsNull() || !align_set_ref->IsSet()) continue;
      RcppThread::checkUserInterrupt();
      
      auto &seq_align_list = align_set_ref->Get(); 
      int64_t parent_list_size = seq_align_list.size();
      
      for (const auto &seq_align : seq_align_list) {
        if (!seq_align || seq_align.IsNull() || !seq_align.NotEmpty()) continue;
        seq_align->Validate(true);
        RcppThread::checkUserInterrupt();
        
        // Clear hoisted strings for buffer reuse
        q_full.clear(); s_full.clear();
        q_aligned.clear(); s_aligned.clear();
        qseq_id.clear(); sseq_id.clear();
        
        try {
          qseq_id = seq_align->GetSeq_id(0).GetSeqIdString(true);
        } catch (...) { qseq_id = "(unknown)"; }
        try {
          sseq_id = seq_align->GetSeq_id(1).GetSeqIdString(true);
        } catch (...) { sseq_id = "(unknown)"; }
        
        // --- Branchless/Direct Strand Construction ---
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
            // GetHSPSequencesFromDenseg(dseg, scope, q_hsp, s_hsp, &q_aligned, &s_aligned
            GetHSPSequencesGracefully(dseg, scope, q_aligned, s_aligned);
          }
        }
        
        // Variable grouping (Primitive types are cheap to initialize locally)
        double score = 0, n_splices = 0, num_ident = 0, aln_len = 0, gaps = 0, mismatches = 0, positive = 0, negative_count = 0;
        double bits = 0, evalue = 0, blast_score = 0, pident = 0, aln_len01 = 0, pident_gap = 0, qcovhsp = 0;
        double sum_evalue = 0, product_coverage = 0, overall_identity = 0, matches = 0, high_quality_percent_coverage = 0, exon_identity = 0, consensus_splices = 0, comp_adj_method = 0;
        
        // NCBI Score Extractions (Unchanged logic, but now writing to local primitives)
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
        static_cast<void>(qhsp_builder.Append(q_aligned));
        static_cast<void>(shsp_builder.Append(s_aligned));
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
        
        // TSeqPos qlen = scope->GetBioseqHandle(seq_align->GetSeq_id(0)).GetBioseqLength();
        // TSeqPos slen = scope->GetBioseqHandle(seq_align->GetSeq_id(1)).GetBioseqLength();
        int qlen = scope->GetSequenceLength(ncbi::sequence::GetId(*qloc.seqloc, scope));
        int slen = scope->GetSequenceLength(ncbi::sequence::GetId(*sloc.seqloc, scope));
        static_cast<void>(qlen_builder.UnsafeAppend(qlen));
        static_cast<void>(slen_builder.UnsafeAppend(slen));
        static_cast<void>(num_alignments_builder.UnsafeAppend(parent_list_size));
        static_cast<void>(hsp_offset_builder.UnsafeAppend(1));
        
        num_rows++;
        
      }
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
    
    if (!aln_struct_array.ok()) {
      Rcpp::Rcerr << "[ExtractHits()] 1. Failed to build StructArray: " << aln_struct_array.status().ToString() << std::endl << std::flush;
      quickblast_running.store(false); 
      return empty_rb;
    }
    
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
    // std::shared_ptr<arrow::StructArray> seqtitle_struct_array = *arrow::StructArray::Make({qseq_title_array, sseq_title_array}, {arrow::field("qseq_title", arrow::utf8()), arrow::field("sseq_title", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqids_struct_array = *arrow::StructArray::Make({qseqid_array, sseqid_array}, {arrow::field("qseqid", arrow::utf8()), arrow::field("sseqid", arrow::utf8())});
    std::shared_ptr<arrow::StructArray> seqs_struct_array = *arrow::StructArray::Make({qseq_array, sseq_array}, {arrow::field("qseq", arrow::large_utf8()), arrow::field("sseq", arrow::large_utf8())});
    std::shared_ptr<arrow::StructArray> lengths_struct_array = *arrow::StructArray::Make({qlen_array, slen_array}, {arrow::field("qlen", arrow::int64()), arrow::field("slen", arrow::int64())});
    
    arrow::Result<std::shared_ptr<arrow::StructArray>> seq_info_array = arrow::StructArray::Make({num_alignment_array,
                                                                                                 seqids_struct_array,
                                                                                                 seqs_struct_array,
                                                                                                 strand_array,
                                                                                                 lengths_struct_array},
                                                                                                 {"num_alignments", "seqids", "seqs", "strands", "lengths"});
    
    if(!seq_info_array.ok()){
      Rcpp::Rcerr << "[ExtractHits()] 2. Failed to build StructArray: " << seq_info_array.status().ToString() << std::endl << std::flush;
      return empty_rb;
    }
    
    std::shared_ptr<arrow::StructArray> seq_info_array_ = seq_info_array.ValueOrDie();
    
    std::shared_ptr<arrow::RecordBatch> alignment_rb = arrow::RecordBatch::Make(arrow_wrapper->GetBLASTSchema(),
                                                                                num_rows,
                                                                                {seq_info_array_, aln_struct_array_});
    
    if(alignment_rb->num_rows() <= 0){
      Rcpp::Rcerr << "ExtractHits() - arrow::RecordBatch() - No alignments could be computed." << std::endl << std::flush;
      return empty_rb;
    }
    
    arrow::Status align_sts = alignment_rb->ValidateFull();
    if(!align_sts.ok()){
      // std::cout << align_sts.message()  << std::endl << align_sts.ToString() << std::endl << "rows:" << alignment_rb->num_rows() << "\ncols:" << alignment_rb->num_columns()  << std::endl << std::flush; //DEBUG
      Rcpp::Rcerr << std::string("ExtractHits() - arrow::RecordBatch() - Alignments failed validation.") + align_sts.detail()->ToString() + "\n" + align_sts.message() << std::endl << std::flush;
      return empty_rb;
    }                                                             
    
    if (alignment_rb)
    {
      const auto &wrt_sts = arrow_wrapper->AddRB2Batch(alignment_rb);
      if (!wrt_sts.ok())
      {
        Rcpp::Rcerr << std::string("ExtractHits() - Error adding RecordBatch to write buffer...") + wrt_sts.detail()->ToString() + "\n" + wrt_sts.message() << std::endl << std::flush; 
        return empty_rb;
      }
      
      if(return_values){
        return alignment_rb;
      }else{
        Rcpp::Rcerr << "ExtractHits() - Empty alignment_rb..." << std::endl << std::flush; //DEBUG
        return empty_rb;
      }
    }
    
    return empty_rb;
  }catch(const ncbi::CException &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: NCBI CException : ")  << e.GetFunction() << std::endl << e.GetErrCodeString() << std::endl << e.GetErrCode() << std::endl << e.GetModule() << std::endl << e.GetPredecessor() << std::endl << e.GetFile() << std::endl << e.GetLine() << std::endl << e.GetMsg() << std::endl << e.GetStackTrace() << std::endl << e.GetStackTraceLevel() << std::endl << e.GetClass() << std::endl << e.what() << std::endl << std::flush;
  }catch(const std::runtime_error &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: C++ Runtime Error : ") + e.what() << std::endl << std::flush;
  }catch(const Rcpp::exception &e){
    Rcpp::Rcerr << std::string("[ExtractHits()] - Rcpp Exception : ") + e.what() << std::endl << std::flush;
  }catch(const std::exception &e){
    Rcpp::Rcerr << std::string("[ExtractHits()]: C++ Exception : ") + e.what() << std::endl << std::flush;
  }catch(...){
    Rcpp::Rcerr << "[ExtractHits()]: Unknown Exception" << std::endl << std::flush;
  }
  return empty_rb;
}

// Return true on success, writes full sequence text to out_seq
bool QuickBLAST::Impl::GetFullSequenceString(CRef<CSeq_id> id, ncbi::CRef<ncbi::objects::CScope> scope, std::string &out_seq)
{
  if (!id) return false;
  CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(*id);
  if (!idh) return false;
  
  CBioseq_Handle bh = scope->GetBioseqHandle(idh);
  if (!bh) return false;
  
  // If the CBioseq has explicit Seq-data in IUPAC, read it directly.
  if (bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_inst::TSeq_data &sd = bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    // If it's IUPACaa or IUPACna, extract directly:
    if (sd.IsIupacna()) {
      // sd.GetIupacna() returns a CSeq_data::TIupacna (string-like)
      out_seq.assign(sd.GetIupacna().Get());
      return true;
    } else if (sd.IsIupacaa()) {
      out_seq.assign(sd.GetIupacaa().Get());
      return true;
    }
    // otherwise fall back to SeqVector (it will do conversions)
  }
  
  // Fallback: use CSeqVector and force IUPAC coding so we get readable letters.
  CSeqVector sv = bh.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
  out_seq.clear();
  TSeqPos len = bh.GetInst().GetLength();
  if (len == 0) return false;
  
  // Fetch in chunks (use pos, stop overload — stop is inclusive)
  const TSeqPos block = 65536;
  for (TSeqPos pos = 0; pos < len; pos += block) {
    TSeqPos fetch_len = std::min(block, len - pos);
    TSeqPos stop = pos + fetch_len - 1;
    std::string chunk;
    sv.GetSeqData(pos, stop, chunk); // returns ASCII IUPAC letters now
    out_seq.append(chunk);
  }
  return true;
}

bool QuickBLAST::Impl::GetHSPSequencesFromDenseg(const CDense_seg& dseg, ncbi::CRef<ncbi::objects::CScope> scope,
                                                 std::string &q_hsp_ungapped, std::string &s_hsp_ungapped,
                                                 std::string *q_aligned_with_gaps, std::string *s_aligned_with_gaps)
{
  q_hsp_ungapped.clear();
  s_hsp_ungapped.clear();
  if (q_aligned_with_gaps) q_aligned_with_gaps->clear();
  if (s_aligned_with_gaps) s_aligned_with_gaps->clear();
  
  const size_t num_rows = dseg.CheckNumRows();
  const size_t num_segs = dseg.GetNumseg();
  if (num_rows < 2 || num_segs == 0) { return false; }
  
  const size_t q_row = 0, s_row = 1;
  if (!dseg.CanGetIds()) { return false; }
  const auto &ids = dseg.GetIds();
  if (ids.size() <= q_row || ids.size() <= s_row) { return false; }
  
  CSeq_id_Handle q_idh = CSeq_id_Handle::GetHandle(*ids[q_row]);
  CSeq_id_Handle s_idh = CSeq_id_Handle::GetHandle(*ids[s_row]);
  if (!q_idh || !s_idh) { return false; }
  
  CBioseq_Handle q_bh = scope->GetBioseqHandle(q_idh);
  CBioseq_Handle s_bh = scope->GetBioseqHandle(s_idh);
  if (!q_bh || !s_bh) {
    // std::cout << "Fallback: direct handle failed for ids: "
    // << q_idh.AsString() << " / " << s_idh.AsString() << std::endl;
    return false;
  }
  
  // Extract the stored IUPAC strings (if present) — prefer these
  std::string stored_q, stored_s;
  if (q_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_data &sd = q_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    if (sd.IsIupacaa()) stored_q = sd.GetIupacaa().Get();
    else if (sd.IsIupacna()) stored_q = sd.GetIupacna().Get();
  }
  if (s_bh.GetCompleteBioseq()->GetInst().CanGetSeq_data()) {
    const CSeq_data &sd = s_bh.GetCompleteBioseq()->GetInst().GetSeq_data();
    if (sd.IsIupacaa()) stored_s = sd.GetIupacaa().Get();
    else if (sd.IsIupacna()) stored_s = sd.GetIupacna().Get();
  }
  // std::cout << "Stored q len: " << stored_q.size() << " s len: " << stored_s.size() << std::endl << std::flush;
  
  CSeqVector qsv(q_bh);
  CSeqVector ssv(s_bh);
  
  qsv.SetIupacCoding(); 
  ssv.SetIupacCoding(); 
  
  const auto &starts = dseg.GetStarts();
  const auto &lens = dseg.GetLens();
  
  for (size_t seg = 0; seg < num_segs; ++seg) {
    TSeqPos seg_len = lens[seg];
    if (seg_len == 0) continue;
    
    ssize_t q_start = starts[seg * num_rows + q_row];
    ssize_t s_start = starts[seg * num_rows + s_row];
    
    // First try: use stored string for exact ASCII residues (recommended)
    if (q_start >= 0 && stored_q.size() >= (size_t)q_start + seg_len) {
      std::string qchunk = stored_q.substr((size_t)q_start, (size_t)seg_len);
      q_hsp_ungapped.append(qchunk);
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(qchunk);
      // std::cout << "Ungapped qchunk size:" << qchunk.size() << std::endl << std::flush;
    } else if (q_start >= 0) {
      // fallback: pull from CSeqVector (might be encoded; handle with care)
      std::string qchunk;
      qsv.GetSeqData((TSeqPos)q_start, (TSeqPos)seg_len, qchunk); // from,length
      q_hsp_ungapped.append(qchunk);
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(qchunk);
      // std::cout << "Ungapped qchunk (from CSeqVector) size:" << qchunk.size() << std::endl << std::flush;
    } else {
      if (q_aligned_with_gaps) q_aligned_with_gaps->append(seg_len, '-');
    }
    
    if (s_start >= 0 && stored_s.size() >= (size_t)s_start + seg_len) {
      std::string schunk = stored_s.substr((size_t)s_start, (size_t)seg_len);
      s_hsp_ungapped.append(schunk);
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(schunk);
      // std::cout << "Ungapped schunk size:" << schunk.size() << std::endl << std::flush;
    } else if (s_start >= 0) {
      std::string schunk;
      ssv.GetSeqData((TSeqPos)s_start, (TSeqPos)seg_len, schunk);
      s_hsp_ungapped.append(schunk);
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(schunk);
      
    } else {
      if (s_aligned_with_gaps) s_aligned_with_gaps->append(seg_len, '-');
    }
  }
  
  return true;
}

std::shared_ptr<arrow::RecordBatchVector> QuickBLAST::BLAST_files(const std::string &queryFile, const std::string &subjectFile, std::string &outFile, const std::string &outFormat, unsigned int num_threads, const bool return_values, int batch_size, bool verbose)
{
  return pImpl->BLAST_files(queryFile, subjectFile, outFile, outFormat, num_threads, return_values, batch_size, verbose);
}

std::shared_ptr<arrow::RecordBatch> QuickBLAST::BLAST_seqs(const std::string &query, const std::string &subject, bool verbose)
{
  return pImpl->BLAST_seqs(query, subject, verbose);
}

auto QuickBLAST::BLAST(const std::string &query, const std::string &subject, std::string &outputFile, const std::string &outFormat, QuickBLAST::EInputType input_type, bool verbose)
{
  return pImpl->BLAST(query, subject, outputFile, outFormat, input_type, verbose);
}

QuickBLAST::QuickBLAST(QuickBLAST::ESeqType seq_type, QuickBLAST::EStrand strand, std::string program, std::string options, bool save_sequences, bool save_hsp_sequences)
  : pImpl(std::make_unique<Impl>(seq_type, strand, program, options, save_sequences, save_hsp_sequences)) {}

// Destructor: default implementation (pImpl will automatically clean up)
QuickBLAST::~QuickBLAST() = default;



