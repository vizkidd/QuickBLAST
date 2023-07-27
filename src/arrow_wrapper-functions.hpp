
ArrowWrapper::ArrowWrapper()
{

  ok_promise.set_value(arrow::Status::OK());

  outputStream = std::make_shared<std::ostringstream>();

  arrow_LFS = arrow::fs::LocalFileSystem();

  parquet::WriterProperties::Builder props_bldr;
  parquet::ArrowWriterProperties::Builder arrow_props_bldr;

  std::string username = getlogin();
  username += "(QuickBLAST)";

  props_bldr.compression(arrow::Compression::LZ4_FRAME);
  props_bldr.created_by(username);
  props_bldr.data_page_version(parquet::ParquetDataPageVersion::V2);
  props_bldr.write_batch_size(1024);
  props_bldr.encoding(parquet::Encoding::RLE);
  props_bldr.version(parquet::ParquetVersion::PARQUET_2_LATEST);
  props_bldr.enable_write_page_index();

  arrow_props_bldr.set_engine_version(parquet::ArrowWriterProperties::EngineVersion::V2);
  arrow_props_bldr.set_use_threads(true);

  parquet_writer_props = props_bldr.build();
  arrow_writer_props = arrow_props_bldr.build();

  ipc_options.allow_64bit = false;
  ipc_options.use_threads = true;
  ipc_options.metadata_version = arrow::ipc::MetadataVersion::V5;
  ipc_options.codec = arrow::util::Codec::Create(arrow::Compression::LZ4_FRAME).ValueOrDie();

  ipc_options.write_legacy_ipc_format = false;

  rbv_batch = std::make_shared<arrow::RecordBatchVector>();
  blast_metadata = std::make_shared<arrow::KeyValueMetadata>();
  AddFASTAMetadata("format", "Arrow Parquet");
  AddFASTAMetadata("Created By", username);
  fasta_schema = arrow::schema({arrow::field("index", arrow::int32()), arrow::field("header", arrow::utf8()), arrow::field("sequence", arrow::utf8())});

  seq_info_type = arrow::struct_({
      arrow::field("num_alignments", arrow::int8()),
      arrow::field("seqids", arrow::struct_({arrow::field("qseqid", arrow::utf8()),
                                             arrow::field("sseqid", arrow::utf8())})),
      arrow::field("seqs", arrow::struct_({arrow::field("qseq", arrow::large_utf8()),
                                           arrow::field("sseq", arrow::large_utf8())})),
      arrow::field("strands", arrow::utf8()),
      arrow::field("lengths", arrow::struct_({arrow::field("qlen", arrow::int8()),
                                              arrow::field("slen", arrow::int8())})),
  });
  hsp_type = arrow::struct_({arrow::field("pident", arrow::float64()),
                             arrow::field("pident_gap", arrow::float64()),
                             arrow::field("frames", arrow::utf8()),
                             arrow::field("evalue", arrow::float64()),
                             arrow::field("length", arrow::int8()),
                             arrow::field("length01", arrow::float64()),
                             arrow::field("qstart", arrow::int8()),
                             arrow::field("qend", arrow::int8()),
                             arrow::field("sstart", arrow::int8()),
                             arrow::field("send", arrow::int8()),
                             arrow::field("bitscore", arrow::float64()),
                             arrow::field("score", arrow::float64()),
                             arrow::field("qcovhsp", arrow::float64()),
                             arrow::field("blast_score", arrow::float64()),
                             arrow::field("gaps", arrow::int8()),
                             arrow::field("nident", arrow::int8()),
                             arrow::field("mismatch", arrow::int8()),
                             arrow::field("positive", arrow::int8()),
                             arrow::field("n_splices", arrow::int8()),
                             arrow::field("hsp_num", arrow::int8())});
  alignment_scores_type = arrow::list({hsp_type});

  blast_schema = arrow::schema({arrow::field("seq_info", seq_info_type),
                                arrow::field("hsps", hsp_type)});
#ifdef _OPENMP
  omp_init_lock(&rec_countLock);
  omp_init_lock(&rbv_batchLock);
#endif
}

ArrowWrapper::~ArrowWrapper()
{
#ifdef _OPENMP
  omp_destroy_lock(&rec_countLock);
#endif
  Rcpp::Rcout << "~ArrowWrapper " << std::endl;
}

std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> ArrowWrapper::MMapFile(const std::string_view &filename, const char *delim)
{
  FILE *file_ptr;
  long fileSize;
  char *end_of_file_ptr;
  std::string delim_str(delim);

  file_ptr = fopen(filename.data(), "r");

  if (!file_ptr)
  {
    Rcpp::Rcerr << "Error: Failed to open file: " << filename.data() << std::endl;
    return nullptr;
  }

  // Get the file size
  fileSize = GetFileSize(file_ptr);

  char *fileData_ptr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file_ptr), 0)); // MAP_SHARED

  if (fileData_ptr == MAP_FAILED)
  {
    Rcpp::Rcerr << "Error: Failed to map file " << filename.data() << std::endl;
    fclose(file_ptr);
    return nullptr;
  }

  end_of_file_ptr = fileData_ptr + fileSize;

  return std::make_shared<std::tuple<FILE *, std::shared_ptr<char>, long, char *>>(std::make_tuple(file_ptr, std::shared_ptr<char>(fileData_ptr, [fileSize, file_ptr](char *ptr)
                                                                                                                                   {
    munmap(ptr, fileSize);
    fclose(file_ptr); }),
                                                                                                   fileSize, end_of_file_ptr));
}

long ArrowWrapper::GetFileSize(FILE *file_ptr)
{
  fseek(file_ptr, 0, SEEK_END);
  long fileSize = ftell(file_ptr);
  rewind(file_ptr);
  return fileSize;
}

void ArrowWrapper::FinishOutputStream()
{
  if (rbv_batch->size() > 0)
  {
    WriteBatch2File();
  }
  rec_writer->Close();
}

arrow::Status ArrowWrapper::WriteBatch2File()
{
  arrow::RecordBatchVector rbv_buffer;
#ifdef _OPENMP
  omp_set_lock(&rbv_batchLock);
#endif
  rbv_buffer.swap(*rbv_batch);
#ifdef _OPENMP
  omp_unset_lock(&rbv_batchLock);
#endif
  for (const auto &rb : rbv_buffer)
  {
    if (rb)
    {
      if (rb->num_rows() > 0)
      {
        arrow::Status rb_sts = rb->ValidateFull();
        if (rb_sts.ok())
        {

          arrow::Status sts = rec_writer->WriteRecordBatch(*rb);
          if (!sts.ok())
          {
            return sts;
          }
        }
        else
        {
          Rcpp::Rcerr << "Warn : Invalid Alignment RB (Not Writing) : " << rb_sts.detail() << std::endl
                      << rb_sts.message() << std::endl
                      << rb->schema()->ToString() << std::endl;
        }
      }
    }
    else
    {
      Rcpp::Rcerr << "ERR : Invalid Alignment RB Ptr (Not Writing)..." << std::endl;
    }
  }

  return arrow::Status::OK();
}

void CountCharacter_thread(const std::string &filename, char character, std::atomic<int> &count, size_t start, size_t end)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    Rcpp::Rcerr << "Failed to open the file." << std::endl;
    return;
  }

  file.seekg(start, std::ios::beg);
  size_t chunkSize = end - start;
  std::vector<char> buffer(chunkSize);

  file.read(buffer.data(), chunkSize);

  for (size_t i = 0; i < chunkSize; ++i)
  {
    if (buffer[i] == character)
    {
      ++count;
    }
  }
}

int ArrowWrapper::CountCharacter(std::string filename, char character, int num_threads)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file)
  {
    Rcpp::Rcerr << "Failed to open the file." << std::endl;
    return 1;
  }

  file.seekg(0, std::ios::end);
  size_t fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  std::atomic<int> totalCount(0);
  std::vector<std::thread> threads;

  size_t chunkSize = fileSize / num_threads;
  size_t start = 0;
  size_t end = chunkSize;

  for (int i = 0; i < num_threads - 1; ++i)
  {
    threads.emplace_back(CountCharacter_thread, filename, character, std::ref(totalCount), start, end);
    start = end;
    end += chunkSize;
  }

  // The last thread might have a slightly larger chunk
  threads.emplace_back(CountCharacter_thread, filename, character, std::ref(totalCount), start, fileSize);

  // Wait for all threads to finish
  for (auto &thread : threads)
  {
    thread.join();
  }

  return totalCount;
}

template <typename T1>
std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::SplitFilesIntoEntries(const std::string_view &filename, const char *delim, const int &num_threads, const std::function<std::shared_ptr<arrow::RecordBatchVector>(std::shared_ptr<T1>)> &Entry_callback)
{

  if constexpr (!std::is_same_v<T1, std::string> || !std::is_same_v<T1, FastaSequenceData>)
  {
    static_assert(std::is_same_v<T1, T1>, "Unsupported type, only std::string & FastaSequenceData are supported");
  }

  std::string delim_str(delim);

  std::shared_ptr<std::tuple<FILE *, std::shared_ptr<char>, long, char *>> file_ptrs = MMapFile(filename, delim);

  char *start_of_file = std::get<1>(*file_ptrs).get();
  char *end_of_file = std::get<3>(*file_ptrs);

  char *p = start_of_file;

  arrow::RecordBatchVector ret_results;

#ifdef _OPENMP
  omp_lock_t pLock;
  omp_lock_t ret_resultsLock;
  omp_init_lock(&pLock);
  omp_init_lock(&ret_resultsLock);
#endif

#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads) shared(end_of_file, start_of_file, delim) // entry_ptr_vec
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait // schedule(dynamic)
#endif
    for (int i = 0; i < num_threads; ++i)
    {

      assert(!Progress::check_abort());
      Rcpp::checkUserInterrupt();
      // Get the thread-specific range to process
      size_t chunk_size = (end_of_file - start_of_file) / num_threads;
      char *thread_start = start_of_file + i * chunk_size;
      char *thread_end = (i == num_threads - 1) ? end_of_file : (thread_start + chunk_size);

      if (thread_start != start_of_file)
      {
        while (strncmp(thread_start, delim, strlen(delim)) != 0)
        {
          --thread_start;
        }
      }
      if (thread_end != end_of_file)
      {
        while (strncmp(thread_end, delim, strlen(delim)) != 0)
        {
          ++thread_end;
        }
        thread_end = thread_end - 1;
      }
      // Process the entries within the thread's range
      char *entryStart = nullptr;
      char *entryEnd = nullptr;

      for (char *p = thread_start; p < thread_end; ++p)
      {
        assert(!Progress::check_abort());

        if (strncmp(p, delim, strlen(delim)) == 0)
        {
          entryStart = strstr(p, delim);
          entryEnd = strstr(p + 1, delim);

          if (entryEnd == nullptr)
          {
            // Handle the case where the delimiter is not found within the thread's range
            entryEnd = end_of_file - 1;
          }

          // Process the entry from entryStart to entryEnd
          std::string full_entry(entryStart, entryEnd - entryStart - 1);
          if (!full_entry.empty())
          {
            AddRecordCount();
            const T1 conv_entry = CastToType<T1>(full_entry);
            if (Entry_callback != nullptr)
            {
              std::shared_ptr<arrow::RecordBatchVector> tmp_result = Entry_callback(std::make_shared<T1>(conv_entry));
#ifdef _OPENMP
              omp_set_lock(&ret_resultsLock);
#endif
              ret_results.insert(ret_results.end(), tmp_result->begin(), tmp_result->end());
#ifdef _OPENMP
              omp_unset_lock(&ret_resultsLock);
#endif
            }
          }
          p = entryEnd - 1; // Move to the next position after the delimiter
        }
      }
    }
#ifdef _OPENMP
#pragma omp barrier
#endif
  }

#ifdef _OPENMP
  omp_destroy_lock(&pLock);
  omp_destroy_lock(&ret_resultsLock);
#endif

  return std::make_shared<arrow::RecordBatchVector>(ret_results);
}

template <typename T>
T ArrowWrapper::CastToType(const std::string &full_entry)
{

  if constexpr (std::is_same_v<T, std::string>)
  {
    return full_entry;
  }
  else if constexpr (std::is_same_v<T, FastaSequenceData>)
  {
    std::regex pattern("^>(\\w+)[\\r\\n]+([\\w\\W]+)$");
    std::regex punct_pattern("[[:punct:][:space:]\\n\\t]");
    std::smatch match;

    std::string full_entry_str(static_cast<std::string>(full_entry));
    FastaSequenceData fasta_data;
    fasta_data.rec_no = GetRecordCount();
    if (std::regex_match(full_entry, match, pattern))
    {
      fasta_data.header = match[1].str();
      fasta_data.seq = std::regex_replace(match[2].str(), punct_pattern, "");
    }
    else
    {

      fasta_data.seq = std::regex_replace(full_entry, punct_pattern, "");
    }
    fasta_data.header = fasta_data.header.empty() ? std::to_string(fasta_data.rec_no) : fasta_data.header;
    AddRecordCount();
    return fasta_data;
  }
  else if constexpr (std::is_same_v<T, int>)
  {
    return std::stoi(full_entry);
  }
  else if constexpr (std::is_same_v<T, long>)
  {
    return std::stol(full_entry);
  }
  else if constexpr (std::is_same_v<T, char>)
  {
    return static_cast<T>(full_entry.front());
  }
  else
  {
    // Handle unsupported types here
    static_assert(std::is_same_v<T, T>, "Unsupported type");
  }
}

int ArrowWrapper::GetColumnCount(const std::string_view &filename, char delim)
{
  std::ifstream file(filename.data());
  if (!file.is_open())
  {
    Rcpp::Rcerr << "Failed to open file: " << filename << std::endl;
    return -1;
  }

  std::string line;
  if (std::getline(file, line))
  {
    std::stringstream ss(line);
    std::string column;
    int count = 0;
    while (std::getline(ss, column, delim))
    {
      count++;
    }
    return count;
  }
  else
  {
    Rcpp::Rcerr << "File is empty: " << filename << std::endl;
    return -1;
  }
}

std::shared_ptr<arrow::RecordBatch> ArrowWrapper::ReadRecordBatchVector(const std::string &file)
{
  std::shared_ptr<arrow::io::ReadableFile> infile = arrow::io::ReadableFile::Open(file).ValueOrDie();
}
