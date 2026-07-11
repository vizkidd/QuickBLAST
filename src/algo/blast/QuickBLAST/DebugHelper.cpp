#include <algo/blast/QuickBLAST/DebugHelper.hpp>

std::string DebugHelper::SafeSeqIdString(const CSeq_id* idptr) {
  try {
    if (!idptr) return "<null id>";
    return idptr->GetSeqIdString(true);
  } catch (...) {
    return "<seqid-error>";
  }
}

std::string DebugHelper::SafeGetNamedScoreAsString(const CSeq_align& align, CSeq_align::EScoreType scoreType) {
  try {
    double val = 0.0;
    // Many scores are double; use try-get pattern
    align.GetNamedScore(scoreType, val);
    std::ostringstream oss;
    oss << val;
    return oss.str();
  } catch (...) {
    return "<n/a>";
  }
}

long DebugHelper::SafeGetNamedScoreInt(const CSeq_align& align, CSeq_align::EScoreType scoreType, bool &ok) {
  ok = false;
  try {
    int value = 0;
    align.GetNamedScore(scoreType, value);
    ok = true;
    return value;
  } catch (...) {
    ok = false;
    return 0;
  }
}

// Print a TSeqLocVector (queries/subjects as used by BLAST functions). Good to check IDs/strands/lengths.
static void DebugHelper::PrintTSeqLocVector(const TSeqLocVector &vec) {
  cout << "TSeqLocVector: size=" << vec.size() << "\n";
  for (size_t i = 0; i < vec.size(); ++i) {
    const auto &sseqloc = vec[i];
    cout << "  [" << i << "] ";
    try {
      if (!sseqloc.seqloc) {
        cout << "seqloc=null\n";
        continue;
      }
      const CSeq_loc* loc = sseqloc.seqloc.GetPointerOrNull();
      if (!loc) {
        cout << "seqloc pointer null\n";
        continue;
      }
      // Seq id (best effort)
      const CSeq_id* idptr = nullptr;
      try {
        idptr = loc->GetId();
      } catch (...) { idptr = nullptr; }
      cout << "id=" << DebugHelper::SafeSeqIdString(idptr);
      
      // length/strand: try Whole or interval if available
      try {
        // Many CSeq_loc objects represent whole or intervals; this is a best-effort print
        if (loc->Which() == CSeq_loc::e_Whole) {
          cout << " (whole)";
        }
      } catch (...) {}
      
      try {
        auto strand = loc->GetStrand();
        cout << " strand=" << static_cast<int>(strand);
      } catch (...) {
        cout << " strand=<n/a>";
      }
    } catch (const ncbi::CException &e) {
      cout << "  [NCBI exception reading SSeqLoc] " << e.GetMsg();
    } catch (const std::exception &e) {
      cout << "  [std::exception] " << e.what();
    } catch (...) {
      cout << "  [unknown error reading SSeqLoc]";
    }
    cout << "\n";
  }
  cout << flush;
}

// Print a TSeqAlignVector (the common BLAST result container).
// Will print a summary of each CSeq_align_set and a few details for each CSeq_align inside.
static void DebugHelper::PrintTSeqAlignVector(const TSeqAlignVector &alignments, size_t max_align_sets = 10, size_t max_aligns_per_set = 10) {
  cout << "TSeqAlignVector: align_sets=" << alignments.size() << "\n";
  for (size_t ai = 0; ai < alignments.size() && ai < max_align_sets; ++ai) {
    const auto &align_set_ref = alignments[ai];
    cout << "  AlignSet[" << ai << "]: ";
    if (!align_set_ref) {
      cout << "<null align_set>\n";
      continue;
    }
    try {
      // The CSeq_align_set contains a list of CRef<CSeq_align>
      const auto &list = align_set_ref->Get();
      cout << " seq_align_count=" << list.size() << "\n";
      
      size_t idx = 0;
      for (const auto &seq_align_ref : list) {
        if (idx++ >= max_aligns_per_set) {
          cout << "    ... (truncated " << (list.size() - max_aligns_per_set) << " more)\n";
          break;
        }
        if (!seq_align_ref) {
          cout << "    SeqAlign: <null>\n";
          continue;
        }
        const CSeq_align &align = *seq_align_ref;
        
        cout << "    SeqAlign[" << (idx-1) << "]:\n";
        // seq ids for the subjects and queries involved (best-effort)
        try {
          const CSeq_id* id0 = nullptr;
          const CSeq_id* id1 = nullptr;
          try { id0 = align.GetSeq_id(0); } catch(...) { id0 = nullptr; }
          try { id1 = align.GetSeq_id(1); } catch(...) { id1 = nullptr; }
          cout << "      qid=" << SafeSeqIdString(id0) << " sid=" << SafeSeqIdString(id1) << "\n";
        } catch (...) {
          cout << "      seq ids: <error>\n";
        }
        
        // starts/stops and strand - wrap in try/catch
        try {
          int qstart = align.GetSeqStart(0);
          int qstop  = align.GetSeqStop(0);
          int sstart = align.GetSeqStart(1);
          int sstop  = align.GetSeqStop(1);
          cout << "      q:start=" << qstart << " q:stop=" << qstop
               << " s:start=" << sstart << " s:stop=" << sstop << "\n";
        } catch (...) {
          cout << "      start/stop: <n/a>\n";
        }
        
        // Named scores (try common ones)
        try {
          bool ok;
          long nident = SafeGetNamedScoreInt(align, CSeq_align::EScoreType::eScore_IdentityCount, ok);
          cout << "      identity_count=" << (ok ? std::to_string(nident) : "<n/a>") ;
        } catch (...) { cout << "      identity_count=<n/a>"; }
        cout << "  ";
        
        try {
          string evalue = SafeGetNamedScoreAsString(align, CSeq_align::EScoreType::eScore_EValue);
          cout << " evalue=" << evalue;
        } catch (...) { cout << " evalue=<n/a>"; }
        
        try {
          string bits = SafeGetNamedScoreAsString(align, CSeq_align::EScoreType::eScore_BitScore);
          cout << " bitscore=" << bits;
        } catch (...) { cout << " bitscore=<n/a>"; }
        
        cout << "\n";
      } // end seq_align loop
      
    } catch (const ncbi::CException &e) {
      cout << "    [NCBI exception iterating align_set] " << e.GetMsg() << "\n";
    } catch (const std::exception &e) {
      cout << "    [std::exception iterating align_set] " << e.what() << "\n";
    } catch (...) {
      cout << "    [unknown exception iterating align_set]\n";
    }
  } // end align_set loop
  cout << flush;
}

// Pretty-print an Arrow RecordBatch. Uses arrow::PrettyPrint with a row window.
static void DebugHelper::PrintRecordBatch(const std::shared_ptr<arrow::RecordBatch> &rb, int max_rows = 10) {
  if (!rb) {
    cout << "RecordBatch: <null>\n";
    return;
  }
  cout << "RecordBatch: rows=" << rb->num_rows() << " cols=" << rb->num_columns() << "\n";
  try {
    if (rb->schema()) {
      cout << "Schema: " << rb->schema()->ToString() << "\n";
    }
  } catch (...) {
    cout << "Schema: <error printing schema>\n";
  }
  
  // Try PrettyPrint (best convenient human readable format), fallback to column-by-column printing.
  arrow::PrettyPrintOptions opts;
  opts.indent = 2;
  opts.window = max_rows;     // only print window rows for big batches
  Status st = PrettyPrint(*rb, opts, &cout);
  if (!st.ok()) {
    cout << "arrow::PrettyPrint failed: " << st.ToString() << "\n";
    // fallback: print first few cells manually
    int nrows = std::min<int>(rb->num_rows(), max_rows);
    cout << "First " << nrows << " rows (fallback):\n";
    for (int c = 0; c < rb->num_columns(); ++c) {
      cout << " Column[" << c << "] name=" << rb->column_name(c) << " type=" << rb->column(c)->type()->ToString() << "\n";
    }
  }
  cout << flush;
}
