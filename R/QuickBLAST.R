#' Globals

#' Get file path of file inside "libs" folder of package
#'
#' @param file_name file_name inside "libs" folder
#' @return File path inside QuickBLAST package
#' @md
.GetLibsPath <- function(file_name = "") {
  libs_path <- file.path("libs", Sys.getenv("R_ARCH"), file_name)
  return(system.file(libs_path, package = "QuickBLAST", mustWork = T))
}

#' Get file path of file inside "bin" folder of package
#'
#' @param file_name file_name inside "bin" folder
#' @return File path inside QuickBLAST package
#' @md
.GetBinPath <- function(file_name = "") {
  bin_path <- file.path("bin", Sys.getenv("R_ARCH"), file_name)
  return(system.file(bin_path, package = "QuickBLAST", mustWork = T))
}

# # arrow_lfs <- arrow::LocalFileSystem$create()
# ## internal single-file-system getter for the package
# .arrow_fs_singleton <- local({
#   fs <- NULL
#   function() {
#     if (!is.null(fs)) return(fs)
#
#     fs <- tryCatch(
#       {
#         arrow::LocalFileSystem$create()
#       },
#       error = function(e) {
#         # fallback if the "file" scheme is already registered (common when arrow already initialized)
#         if (grepl("Attempted to register factory for scheme 'file'", e$message, fixed = TRUE)) {
#           arrow::FileSystem$from_uri("file:///")
#         } else {
#           stop(e)
#         }
#       }
#     )
#
#     # cache and return
#     force(fs)
#     fs
#   }
# })

# #' Get an Instance of QuickBLAST class and its exposed methods
# #' @note Check BLAST C++ Call in Help for the list of parameters for the exposed BLAST function. Exposed C++ function only takes BLAST options as string.
# #' @seealso [QuickBLAST::CreateQuickBLASTInstance()], [QuickBLAST::GetQuickBLASTEnums()]
# #' @examples
# #' \dontrun{
# #' bw_obj <- QuickBLAST::GetQuickBLASTInstance(list(0, 0, FALSE), "blastn", "-evalue 1e-05")
# #' bw_obj$BLAST("ungrouped.cds", "ungrouped.cds", "out.tmp", 0, 1000, TRUE)
# #' bw_obj$BLAST("AAAAAAAAAAAATTTTTTTTTTTTGGGGGGGGGGGCCCCCCCCC", "TTTTTTTTTTTGGGGGGGGGGGG", "", 1, 1000, FALSE)
# #' }
# #'
# #' @param seq_info Ordered List of 1) (int) Sequence Type, 2) (int) Strand (bool), 3) Save Sequences in BLAST Hits? : TRUE - BLAST Hits have sequences, FALSE - FASTA sequences are not stored in BLAST Hits. Check QuickBLAST::GetQuickBLASTEnums() for available enums
# #' @param program (string) Name of the BLAST program
# #' @param options (string) String of BLAST options - eg, "-evalue 1e-05 -pident 0.75" - check QuickBLAST::GetAvailableBLASTOptions() for available options
# #' @return A new QuickBLAST Object
# #' @md
# #' @export
# GetQuickBLASTInstance <- function(seq_info, program, options) {
#   mod <- Rcpp::Module("blast_module", inline::getDynLib("QuickBLAST"))
#   return(methods::new(mod$QuickBLAST, seq_info[[1]], seq_info[[2]], program, options, seq_info[[3]]))
# }


#' Get a list of Enums used by QuickBLAST
#'
#' @return A List of Enums used by QuickBLAST
#' @examples
#' enums <- QuickBLAST::GetQuickBLASTEnums()
#' print(names(enums))
#' print(enums$ESeqType$eNucleotide)
#' @importFrom Rcpp evalCpp
#' @useDynLib QuickBLAST, .registration=TRUE
#' @export
GetQuickBLASTEnums <- function() {
  warning("Check https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/ident to search/lookup all the enums.")
  return(list("ESeqType" = list(eNucleotide = 0, eProtein = 1), "EStrand" = list(
    ePlus = 0,
    eMinus = 1,
    eBoth = 2,
    eBoth_rev = 3,
    eOther = 4,
    eUnknown = 5
  ), "EInputType" = list(eFile = 0, eSequenceString = 1, eFolder = 2)))
}

#' Get a List of Available BLAST options
#'
#'  Use this function in blast_options to set BLAST defaults for the chosen BLAST program.
#'
#' @note CREATE a NEW LIST with ONLY the OPTIONS THAT YOU NEED
#' @examples
#' opts <- QuickBLAST::GetAvailableBLASTOptions()
#' print(opts)
#' @return A List of Available BLAST options
#' @export
GetAvailableBLASTOptions <- function() {
  warning("These are Defaults, CREATE a NEW LIST with ONLY the OPTIONS THAT YOU NEED (or) Use this function in blast_options to set BLAST defaults for the chosen BLAST program. (Options can also be given as a string '-evalue 1 -hit_list_size 10000')")
  blastOptions <- list(
    "evalue" = double(),
    "pident" = double(),
    "gapped_mode" = logical(),
    "filter_string" = character(),
    "effective_search_space" = integer(),
    "cutoff_score" = integer(),
    "gap_trigger" = double(),
    "gap_x_dropoff" = double(),
    "gap_x_dropoff_final" = double(),
    "hit_list_size" = integer(),
    "low_score_percentage" = double(),
    "max_hsp_per_subject" = integer(),
    "max_hsp_per_sequence" = integer(),
    "qcovhsp_perc" = double(),
    "window_size" = integer()
  )
  return(blastOptions)
}

#' Load BLAST Hits into a data.frame
#'
#' Give the path to a BLAST Hits file (table/ipc/parquet) to load it into a tibble (BLAST HITs Table).
#'
#' @param infile BLAST hits filename (not a connection) (Gzipped files supported)
#' @param sep Delimiter of the BLAST File columns. (Only applies when format == 'table'). Default - '\\t'
#' @param header Does the file have a header? (Only applies when format == 'table'). Default - FALSE
#' @param format Input Format (Required) - 'table' (CSV/TSV)/'ipc' (arrow::ipc)/'parquet' (arrow::parquet) - Default : 'parquet'
#' @examples
#' \donttest{
#' # Assuming 'blast_results.parquet' exists in your working directory
#' # results <- QuickBLAST::LoadBLASTHits("blast_results.parquet", format = "parquet")
#' }
#' @return Data Frame with BLAST Results
#' @export
LoadBLASTHits <- function(infile, sep = "\t", header = F, format = "parquet") {
  if (try(file.exists(infile)) && file.info(infile)$size > 0) {
    # any(grepl(x = class(infile), pattern = "gzfile|connection", ignore.case = T)) #
    # if(gzipped){
    #  infile <- gzfile(description = infile, open = "r")
    # }
    # blast_results <- read.table(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, fill = t,na.strings = NA)
    if (format == "table" || format == "tsv" || format == "csv") {
      blast_results <- iterators::iread.table(file = infile, row.names = NULL, header = header, sep = sep, quote = "", blank.lines.skip = T, fill = T, na.strings = "NA") # data.table::fread(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, nThread = n_threads)
      # return(blast_results)
      # Use read_delim_arrow to properly parse the \t (or user-provided) separator
      blast_results <- arrow::read_delim_arrow(
        file = infile,
        delim = sep,
        as_data_frame = TRUE
      )
      return(blast_results)
    } else if (format == "ipc") {
      # # arrow_lfs <- arrow::LocalFileSystem$create()
      # arrow_lfs <- .arrow_fs_singleton()
      # # arrow_i_stream <- arrow_lfs$OpenInputStream(infile)
      # # batch_reader <- arrow::RecordBatchStreamReader$create(arrow_i_stream)
      # in_stream <- arrow_lfs$OpenInputStream(infile)
      # # compressed_stream <- arrow::CompressedInputStream$create(infile, "gzip")
      # # Use the FileReader (not StreamReader)
      # file_reader <- arrow::RecordBatchFileReader$create(in_stream)
      # # file_reader <- arrow::RecordBatchStreamReader$create(compressed_stream)
      # # # create an iterator over 1..n, mapping to get_batch
      # # n <- file_reader$num_record_batches
      # # batch_iter <- iterators::iter(seq_len(n), function(j) {
      # #   return(file_reader$get_batch(j - 1L)$to_data_frame())
      # # })
      # ret_tib <- dplyr::bind_rows(lapply(file_reader$batches(), FUN = function(x){
      #   return(x$to_data_frame())
      # }))
      #
      # ret_df <- ret_tib %>%
      #   # expand seq_info (it is a tibble/list-col)
      #   tidyr::unnest(cols = c(seq_info)) %>%
      #   # expand seqids and seqs and lengths which are themselves nested tibbles
      #   tidyr::unnest_wider(seqids, names_sep = "_") %>%
      #   tidyr::unnest_wider(seqs,   names_sep = "_") %>%
      #   tidyr::unnest_wider(lengths, names_sep = "_") %>%
      #   # expand hsps which is another nested tibble
      #   tidyr::unnest(cols = c(hsps)) %>%
      #   # if any columns are still list-columns of length-1, unnest them:
      #   dplyr::mutate(dplyr::across(dplyr::where(~ is.list(.) && all(lengths(.) == n())), ~ unlist(.))) %>%
      #   tidyr::as_tibble()
      #
      # # blast_results <- arrow::read_feather(file = infile, mmap = T)
      # return(ret_df) # batch_iter <- iter(function())

      # 1. Safely read the IPC/Feather file directly into a tibble
      # This completely bypasses the LocalFileSystem registry bug!
      ret_tib <- arrow::read_ipc_file(file = infile, as_data_frame = TRUE)

      # 2. Check if the data is nested (i.e. 'seq_info' exists)
      if ("seq_info" %in% names(ret_tib)) {
        ret_df <- ret_tib %>%
          # expand seq_info (it is a tibble/list-col)
          tidyr::unnest(cols = c(seq_info)) %>%
          # expand seqids and seqs and lengths which are themselves nested tibbles
          tidyr::unnest_wider(seqids, names_sep = "_") %>%
          tidyr::unnest_wider(seqs, names_sep = "_") %>%
          tidyr::unnest_wider(lengths, names_sep = "_") %>%
          # expand hsps which is another nested tibble
          tidyr::unnest(cols = c(hsps)) %>%
          # if any columns are still list-columns of length-1, unnest them:
          dplyr::mutate(dplyr::across(dplyr::where(~ is.list(.) && all(lengths(.) == dplyr::n())), ~ unlist(.))) %>%
          tidyr::as_tibble()

        return(ret_df)
      } else {
        # If the backend C++ already flattened the structure (like we did for CSV),
        # return it immediately without the unnesting pipeline
        return(ret_tib)
      }
    } else if (format == "parquet") {
      return(arrow::read_parquet(infile))
    } else {
      stop(paste("Format: Should be one of 'table'/'ipc'/'parquet' "))
    }
  } else {
    stop(paste("File", infile, "does not exist or size 0"))
  }
}

#' Execute one2one QuickBLAST
#'
#' Executes One-to-One QuickBLAST between two lists of organisms/genes/clusters. The BLAST Hits are stored in Arrow::Feather/Parquet format.
#'
#' @examples
#' QuickBLAST::one2one(
#'   first_list = fs::path_package("QuickBLAST", "extdata", "protein_query.fasta"),
#'   second_list = fs::path_package("QuickBLAST", "extdata", "protein_subject.fasta"),
#'   blast_fun = QuickBLAST::BLAST2Files,
#'   seq_type = 0,
#'   strand = 0,
#'   output_dir = "./",
#'   n_threads = 8,
#'   blast_program = "tblastx",
#'   save_sequences = FALSE,
#'   save_hsp_sequences = FALSE,
#'   return_values = FALSE,
#'   min_batch_size = 256,
#'   out_format = "parquet",
#'   blast_options = "",
#'   verbose = TRUE
#' )
#'
#' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
#' @param first_list Vector of FASTA Filenames or Strings
#' @param second_list Vector of FASTA Filenames or Strings
#' @param blast_fun One of QuickBLAST::BLAST2Seqs, QuickBLAST::BLAST2Files, QuickBLAST::BLAST2Folders, QuickBLAST::BLAST2DBs
#' @param seq_type (int) Sequence Type. Check QuickBLAST::GetQuickBLASTEnums()$ESeqType for available enums.
#' @param strand (int) Strand. Check QuickBLAST::GetQuickBLASTEnums()$EStrand for available enums.
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program.
#' @param file_ext File extension of input files. eg- ".cds" or ".fa", Unused if input_type is GetQuickBLASTEnums()$EInputType$eSequencesString
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param output_dir Path to BLAST output
#' @param ... Extended options passed to internal functions, including:
#'   \itemize{
#'     \item \code{blast_options}: BLAST Options to use - QuickBLAST::GetAvailableBLASTOptions()
#'     \item \code{save_sequences}: (bool) Save full sequences to result?
#'     \item \code{save_hsp_sequences}: (bool) Save HSP sequences to result?
#'     \item \code{return_values}: (bool) Return values back to R?
#'     \item \code{min_batch_size}: Minimum batch size. (Default: 256)
#'     \item \code{n_threads}: Number of threads. (Default: 8)
#'     \item \code{out_format}: Output format. ipc/csv/parquet (Default: "parquet")
#'     \item \code{extension}: File extension. (Only for QuickBLAST::BLAST2Folders())
#'     \item \code{reciprocal_hits}: (bool) Reciprocal (Bi-directional) Hits?
#'     \item \code{verbose}: (bool) Print DEBUG Messages?
#'   }
#' @return If \code{return_values = TRUE}, returns a list of data frames corresponding to each alignment query. Otherwise, returns \code{invisible(NULL)} or outputs directly to files.
#' @md
#' @export
one2one <- function(first_list, second_list, blast_fun, seq_type, strand, blast_program, file_ext = ".fa", input_prefix_path = NULL, output_dir = "./", ...) {
  dot_args <- list(...)
  blast_options <- ""
  if ("blast_options" %in% names(dot_args)) {
    blast_options <- dot_args$blast_options
  }
  save_sequences <- F
  if ("save_sequences" %in% names(dot_args)) {
    save_sequences <- dot_args$save_sequences
  }
  save_hsp_sequences <- T
  if ("save_hsp_sequences" %in% names(dot_args)) {
    save_hsp_sequences <- dot_args$save_hsp_sequences
  }
  return_values <- F
  if ("return_values" %in% names(dot_args)) {
    return_values <- dot_args$return_values
  }
  min_batch_size <- 256
  if ("min_batch_size" %in% names(dot_args)) {
    min_batch_size <- dot_args$min_batch_size
  }
  n_threads <- tryCatch(parallel::detectCores(all.tests = T, logical = T), error = function(cond) {
    return(2)
  })
  if ("n_threads" %in% names(dot_args)) {
    n_threads <- dot_args$n_threads
  }
  out_format <- "parquet"
  if ("out_format" %in% names(dot_args)) {
    out_format <- dot_args$out_format
  }
  extension <- ".fa"
  if ("extension" %in% names(dot_args)) {
    extension <- dot_args$extension
  }
  reciprocal_hits <- FALSE
  if ("reciprocal_hits" %in% names(dot_args)) {
    reciprocal_hits <- dot_args$reciprocal_hits
  }
  verbose <- T
  if ("verbose" %in% names(dot_args)) {
    verbose <- dot_args$verbose
  }
  future::plan(future.callr::callr, workers = n_threads)

  if (!is.null(input_prefix_path)) {
    first_list <- paste(input_prefix_path, "/", first_list, file_ext, sep = "")
    second_list <- paste(input_prefix_path, "/", second_list, file_ext, sep = "")
  } else {
    first_list <- paste(first_list, file_ext, sep = "")
    second_list <- paste(second_list, file_ext, sep = "")
  }

  list_combos <- unique(tidyr::crossing(first_list[order(first_list)], second_list[order(second_list)]))

  # return_data <- furrr::future_map2(.x=first_list[order(first_list)], .y=second_list[order(second_list)], .f=function(x,y){
  # parallel::mclapply(seq_along(1:nrow(list_combos)), function(idx) {
  future.apply::future_lapply(seq_len(nrow(list_combos)), function(idx) {
    x <- toString(list_combos[idx, 1])
    y <- toString(list_combos[idx, 2])
    if (input_type == QuickBLAST::GetQuickBLASTEnums()$EInputType$eFile) {
      if (!all(file.exists(x), file.exists(y), file.info(x)$size > 0, file.info(y)$size > 0)) {
        warning(paste(x, "or", y, "missing/empty and input_type is EInputType$eFile, assuming input to be sequences!", sep = " "))
        input_type <- 1
      }
    }
    if (verbose) {
      print(x)
      print(y)
    }

    blast_ptr <- QuickBLAST::CreateQuickBLASTInstance(seq_type = seq_type, strand = strand, program = blast_program, save_sequences = save_sequences, save_hsp_sequences = save_hsp_sequences, options = blast_options)

    # switch(input_type,
    #   { # eFile
    #     return(QuickBLAST::BLAST2Files(ptr = blast_ptr, query = x, subject = y, out_file = paste(output_dir, "/", basename(tools::file_path_sans_ext(x)), ".", basename(tools::file_path_sans_ext(y)), ".hits", sep = ""), out_format=out_format, return_values = return_values, min_batch_size = min_batch_size, num_threads = n_threads, verbose=verbose))
    #   },
    #   { # eSequenceString
    #     return(QuickBLAST::BLAST2Seqs(ptr = blast_ptr, query = x, subject = y, verbose=verbose))
    #   },
    #   { # eFolder
    #   }
    # )
    if (identical(blast_fun, QuickBLAST::BLAST2Files)) {
      return(blast_fun(ptr = blast_ptr, query = x, subject = y, out_file = paste(output_dir, "/", basename(tools::file_path_sans_ext(x)), ".", basename(tools::file_path_sans_ext(y)), ".hits", sep = ""), out_format = out_format, return_values = return_values, min_batch_size = min_batch_size, num_threads = n_threads, verbose = verbose))
    } else if (identical(blast_fun, QuickBLAST::BLAST2Seqs)) {
      return(blast_fun(ptr = blast_ptr, query = x, subject = y, verbose = verbose))
    } else if (identical(blast_fun, QuickBLAST::BLAST2Folders)) {
      return(blast_fun(ptr = blast_ptr, query = x, subject = y, extension = extension, out_folder = output_dir, out_format = out_format, reciprocal_hits = reciprocal_hits, min_batch_size = min_batch_size, num_threads = n_threads, verbose = verbose))
    } else if (identical(blast_fun, QuickBLAST::BLAST2DBs)) {
      return(blast_fun(ptr = blast_ptr, query = x, subject = y, out_file = paste(output_dir, "/", basename(tools::file_path_sans_ext(x)), ".", basename(tools::file_path_sans_ext(y)), ".hits", sep = ""), out_format = out_format, return_values = return_values, min_batch_size = min_batch_size, num_threads = n_threads, verbose = verbose))
    } else {
      stop("[one2one()] blast_fun: Should be one of QuickBLAST::BLAST2Files|QuickBLAST::BLAST2Seqs|QuickBLAST::BLAST2Folders|QuickBLAST::BLAST2DBs")
    }
  })
  # }, mc.cores = n_threads, mc.silent = !verbose)

  # return(return_data)
}


#' Execute all2all QuickBLAST
#'
#' Executes All-to-All QuickBLAST between two lists of organisms/genes/clusters. Output BLAST files are bi-directional and are stored in the filename filename1.filename2.all2all under output_dir. (All-to-All is simply Many-to-Many association)
#'
#' @examples
#' QuickBLAST::all2all(
#'   first_list = fs::path_package("QuickBLAST", "extdata", "protein_query.fasta"),
#'   second_list = fs::path_package("QuickBLAST", "extdata", "protein_subject.fasta"),
#'   blast_fun = QuickBLAST::BLAST2Files,
#'   seq_type = 0,
#'   strand = 0,
#'   output_dir = "./",
#'   n_threads = 8,
#'   blast_program = "tblastx",
#'   save_sequences = FALSE,
#'   save_hsp_sequences = FALSE,
#'   return_values = TRUE,
#'   min_batch_size = 256,
#'   out_format = "parquet",
#'   blast_options = "",
#'   verbose = TRUE
#' )
#'
#' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
#' @param first_list Vector of FASTA Filenames or Strings
#' @param second_list Vector of FASTA Filenames or Strings
#' @param blast_fun One of QuickBLAST::BLAST2Seqs, QuickBLAST::BLAST2Files, QuickBLAST::BLAST2Folders, QuickBLAST::BLAST2DBs
#' @param seq_type (int) Sequence Type. Check QuickBLAST::GetQuickBLASTEnums()$ESeqType for available enums.
#' @param strand (int) Strand. Check QuickBLAST::GetQuickBLASTEnums()$EStrand for available enums.
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program.
#' @param file_ext File extension of input files. eg- ".cds" or ".fa", Unused if input_type is GetQuickBLASTEnums()$EInputType$eSequencesString
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param output_dir Path to BLAST output
#' @param ... Extended options passed to internal functions, including:
#'   \itemize{
#'     \item \code{blast_options}: BLAST Options to use - QuickBLAST::GetAvailableBLASTOptions()
#'     \item \code{save_sequences}: (bool) Save full sequences to result?
#'     \item \code{save_hsp_sequences}: (bool) Save HSP sequences to result?
#'     \item \code{return_values}: (bool) Return values back to R?
#'     \item \code{min_batch_size}: Minimum batch size. (Default: 256)
#'     \item \code{n_threads}: Number of threads. (Default: 8)
#'     \item \code{out_format}: Output format. ipc/csv/parquet (Default: "parquet")
#'     \item \code{extension}: File extension. (Only for QuickBLAST::BLAST2Folders())
#'     \item \code{reciprocal_hits}: (bool) Reciprocal (Bi-directional) Hits?
#'     \item \code{verbose}: (bool) Print DEBUG Messages?
#'   }
#' @return If \code{return_values = TRUE}, returns a list of data frames corresponding to each alignment query. Otherwise, returns \code{invisible(NULL)} or outputs directly to files.
#' @md
#' @export
all2all <- function(first_list, second_list, blast_fun, seq_type, strand, blast_program, file_ext = ".fa", input_prefix_path = NULL, output_dir = "./", ...) {
  dot_args <- list(...)
  blast_options <- ""
  if ("blast_options" %in% names(dot_args)) {
    blast_options <- dot_args$blast_options
  }
  save_sequences <- F
  if ("save_sequences" %in% names(dot_args)) {
    save_sequences <- dot_args$save_sequences
  }
  save_hsp_sequences <- T
  if ("save_hsp_sequences" %in% names(dot_args)) {
    save_hsp_sequences <- dot_args$save_hsp_sequences
  }
  return_values <- F
  if ("return_values" %in% names(dot_args)) {
    return_values <- dot_args$return_values
  }
  min_batch_size <- 256
  if ("min_batch_size" %in% names(dot_args)) {
    min_batch_size <- dot_args$min_batch_size
  }
  n_threads <- tryCatch(parallel::detectCores(all.tests = T, logical = T), error = function(cond) {
    return(2)
  })
  if ("n_threads" %in% names(dot_args)) {
    n_threads <- dot_args$n_threads
  }
  out_format <- "parquet"
  if ("out_format" %in% names(dot_args)) {
    out_format <- dot_args$out_format
  }
  extension <- ".fa"
  if ("extension" %in% names(dot_args)) {
    extension <- dot_args$extension
  }
  reciprocal_hits <- FALSE
  if ("reciprocal_hits" %in% names(dot_args)) {
    reciprocal_hits <- dot_args$reciprocal_hits
  }
  verbose <- T
  if ("verbose" %in% names(dot_args)) {
    verbose <- dot_args$verbose
  }
  if (verbose) {
    cat(paste("All2All QuickBLAST Started...", "\n", sep = ""))
    print(paste(first_list, collapse = ","))
    print(paste(second_list, collapse = ","))
  }
  future::plan(future.callr::callr, workers = n_threads)
  dir.create(path = output_dir, recursive = T, showWarnings = F)

  list_combinations <- unique(tidyr::crossing(first_list, second_list))
  # furrr::future_map2(.x=list_combinations$first_list, .y=list_combinations$second_list, .f=function(first_set,second_set){
  parallel::mclapply(seq_len(nrow(list_combinations)), function(idx) {
    # future.apply::future_lapply(seq_along(1:nrow(list_combinations)), function(idx) {
    first_set <- list_combinations[idx, 1]
    second_set <- list_combinations[idx, 2]

    tryCatch(QuickBLAST::one2one(first_list = first_set, second_list = second_set, blast_fun = blast_fun, seq_type = seq_type, strand = strand, file_ext = file_ext, input_prefix_path = input_prefix_path, reciprocal_hits = reciprocal_hits, n_threads = 1, blast_program = blast_program, output_dir = output_dir, blast_options = blast_options, save_sequences = save_sequences, save_hsp_sequences = save_hsp_sequences, return_values = return_values, min_batch_size = min_batch_size, out_format = out_format, extension = extension, verbose = verbose),
      error = function(cond) {
        if (verbose) {
          message(cond)
        }
      }
    )
    tryCatch(QuickBLAST::one2one(first_list = first_set, second_list = second_set, blast_fun = blast_fun, seq_type = seq_type, strand = strand, file_ext = file_ext, input_prefix_path = input_prefix_path, reciprocal_hits = reciprocal_hits, n_threads = 1, blast_program = blast_program, output_dir = output_dir, blast_options = blast_options, save_sequences = save_sequences, save_hsp_sequences = save_hsp_sequences, return_values = return_values, min_batch_size = min_batch_size, out_format = out_format, extension = extension, verbose = verbose),
      error = function(cond) {
        if (verbose) {
          message(cond)
        }
      }
    )
  }, mc.cores = n_threads, mc.preschedule = T)
}

R_dlls <- c("Riconv", "R", "Rgraphapp", "Rblas", "R", "Rlapack")
R_dll_paths <- paste0(file.path(Sys.getenv("R_HOME"), "bin", Sys.getenv("R_ARCH")), .Platform$file.sep, R_dlls, .Platform$dynlib.ext)

# dlls <- c("libncbi_core", "libncbi_general", "libncbi_pub", "libncbi_seq", "libncbi_trackmgr", "libncbi_eutils", "libncbi_misc", "libsqlitewrapp", "liblmdb", "libefetch", "libncbi_seqext", "libncbi_xreader", "libncbi_xreader_id1", "libncbi_xreader_id2", "libncbi_xreader_cache", "libxxconnect2", "libpsg_client", "libncbi_xloader_genbank", "libncbi_xloader_blastdb", "libncbi_xloader_blastdb_rmt", "libncbi_general", "libncbi_web", "libncbi_align_format", "libutrtprof", "libncbi_algo", tools::file_path_sans_ext(list.files(fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH")),pattern = "*arrow.*dll", full.names = F)), "libncbi_blastinput", "libQuickBLASTcpp")

## Generate the DLL dependencies
# lddtree libQuickBLASTcpp.so | awk '{ print $1 }' | grep -f <(ls -1)
dlls <- c("libxncbi", "libxutil", "libxser", "libgeneral", "libseqcode", "libbiblio", "libmedline", "libpub", "libsequtil", "libseq", "libseqset", "libseqsplit", "libseqedit", "libgenome_collection", "libsubmit", "libxobjmgr", "libscoremat", "libxnetblast", "libxconnect", "libxobjutil", "libxlogging", "libxobjread", "libsqlitewrapp", "liblmdb", "libblastdb", "libseqdb", "libtables", "libconnect", "libcomposition_adjustment", "libutrtprof", "libxconnect", "libxnetblastcli", "libseqmasks_io", "libxalgowinmask", "libxalgodustmask", "libblast", "libxalgoblastdbindex", "libxblast", "libarrow", "libQuickBLASTcpp")
# dlls <- c("libncbi_core", "libncbi_general", "libncbi_pub", "libncbi_seq", "libncbi_trackmgr", "libncbi_eutils", "libncbi_misc", "libsqlitewrapp", "liblmdb", "libefetch", "libncbi_seqext", "libncbi_xreader", "libncbi_xreader_id1", "libncbi_xreader_id2", "libncbi_xreader_cache", "libxxconnect2", "libpsg_client", "libncbi_xloader_genbank", "libncbi_web", "libncbi_align_format", "libutrtprof", "libncbi_algo", "libQuickBLASTcpp")


dll_paths <- paste(fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), .Platform$file.sep, dlls, .Platform$dynlib.ext, sep = "")

dll_obj_list <- list()

.onLoad <- function(libname, pkgname) {
  # Sys.setenv("ASAN_OPTIONS"="detect_leaks=0") #ASAN - remove in prod
  # Sys.setenv("LD_PRELOAD"="/usr/lib/x86_64-linux-gnu/libasan.so.8") #ASAN - remove in prod
  RcppThread::detectCores()
  tzdb::tzdb_initialize()

  if (xfun::is_windows()) {
    Sys.setenv("ARROW_DEFAULT_MEMORY_POOL" = "system")
  } else {
    Sys.setenv("ARROW_DEFAULT_MEMORY_POOL" = "jemalloc")
  }
  Sys.setenv("ARROW_DEBUG_MEMORY_POOL" = "warn")
  Sys.setenv("OMP_NUM_THREADS" = parallel::detectCores(all.tests = T, logical = T))
  Sys.setenv("OMP_DYNAMIC" = "TRUE")
  Sys.setenv("OMP_WAIT_POLICY" = "PASSIVE")
  # Sys.setenv("OMP_DISPLAY_ENV"="VERBOSE")
  Sys.setenv("OMP_DISPLAY_ENV" = "FALSE")
  # packageStartupMessage("QuickBLAST Loaded!")
  # packageStartupMessage("Version: ", utils::packageVersion("QuickBLAST"))
  # packageStartupMessage("Github: https://github.com/vizkidd/QuickBLAST")
  # #R_dll_paths, msys_dll_paths
  # # library.dynam("Rcpp", "Rcpp", fs::path_package("Rcpp","..",".."))
  # if(Sys.info()['sysname'] != "Linux"){
  #   for (dll_path in c(R_dll_paths)) {
  #     dyn.load(dll_path,now = T)
  #   }
  # }
  # for (dll_path in c(dll_paths)) {
  # # for (dll_path in c(dlls)) {
  #   if (!file.exists(dll_path)) {
  #     cat("DLL file not found:", dll_path, "\n")
  #   } else {
  #     # dyn.load(dll_path, local=F, now = T)
  #     # library.dynam(dll_path, "QuickBLAST", fs::path_package("QuickBLAST",".."))
  #     # if(!invisible(is.loaded(dll_path))){
  #       dll_obj_list <- append(dll_obj_list, list(dyn.load(dll_path,now = T)))
  #     # }
  #     cat("Loaded DLL:", dll_path, "\n")
  #   }
  # }
}

.onAttach <- function(libname, pkgname) {
  Sys.setenv("ARROW_DEFAULT_MEMORY_POOL" = "system")
  Sys.setenv("ARROW_DEBUG_MEMORY_POOL" = "warn")
  Sys.setenv("OMP_NUM_THREADS" = parallel::detectCores(all.tests = T, logical = T))
  Sys.setenv("OMP_DYNAMIC" = "TRUE")
  Sys.setenv("OMP_WAIT_POLICY" = "PASSIVE")
  Sys.setenv("OMP_DISPLAY_ENV" = "FALSE")
  # Sys.setenv("OMP_DISPLAY_ENV"="VERBOSE")
  packageStartupMessage("QuickBLAST Loaded!")
  packageStartupMessage("Version: ", utils::packageVersion("QuickBLAST"))
  packageStartupMessage("Github: https://github.com/vizkidd/QuickBLAST")
}

# .onUnload() function
.onUnload <- function(libpath) {
  # # Unload the DLLs when the package is unloaded
  # msys_dll_paths
  packageStartupMessage("Unloading QuickBLAST...")

  # loaded_dlls <- getLoadedDLLs()
  # loaded_dlls <- loaded_dlls[na.omit(match(dlls,names(loaded_dlls)))]
  #
  # for(dll_info in loaded_dlls){
  #   try(dyn.unload(dll_info[["path"]]), silent = T)
  # }
  #
  # # for (dll_path in c(rev(c(dll_paths)))) {
  # #   if(is.loaded(dll_path)){
  # #     if (dyn.unload(dll_path)) {
  # #       packageStartupMessage(cat("Unloaded DLL:", dll_path, "\n"))
  # #     } else {
  # #       packageStartupMessage(cat("Failed to unload DLL:", dll_path, "\n"))
  # #     }
  # #   }
  # # }
  # # detach("package:QuickBLAST", unload = TRUE)
}

.onDetach <- function(libpath) {
  packageStartupMessage("Detaching QuickBLAST...")
}

LoadBLASTHits <- compiler::cmpfun(LoadBLASTHits)
GetAvailableBLASTOptions <- compiler::cmpfun(GetAvailableBLASTOptions)
GetQuickBLASTEnums <- compiler::cmpfun(GetQuickBLASTEnums)
one2one <- compiler::cmpfun(one2one)
all2all <- compiler::cmpfun(all2all)

BLAST2Seqs <- compiler::cmpfun(BLAST2Seqs)
BLAST2Folders <- compiler::cmpfun(BLAST2Folders)
BLAST1Folder <- compiler::cmpfun(BLAST1Folder)
BLAST2Files <- compiler::cmpfun(BLAST2Files)
RemoteBLAST <- compiler::cmpfun(RemoteBLAST)
BLAST2DBs <- compiler::cmpfun(BLAST2DBs)
GetFASTAHeaders <- compiler::cmpfun(GetFASTAHeaders)
BLASTFile2DB <- compiler::cmpfun(BLASTFile2DB)
MakeBLASTDB <- compiler::cmpfun(MakeBLASTDB)
isBLASTDB <- compiler::cmpfun(isBLASTDB)
SetQuickBLASTOptions <- compiler::cmpfun(SetQuickBLASTOptions)
GetQuickBLASTOptions <- compiler::cmpfun(GetQuickBLASTOptions)
DeleteQuickBLASTInstance <- compiler::cmpfun(DeleteQuickBLASTInstance)
GetQuickBLASTInstance <- compiler::cmpfun(GetQuickBLASTInstance)
GetInstanceID <- compiler::cmpfun(GetInstanceID)
GetInstanceCount <- compiler::cmpfun(GetInstanceCount)
CreateQuickBLASTInstance <- compiler::cmpfun(CreateQuickBLASTInstance)
RecordBatchVectorToFlattenedDFList <- compiler::cmpfun(RecordBatchVectorToFlattenedDFList)
isQuickBLASTLoaded <- compiler::cmpfun(isQuickBLASTLoaded)

utils::globalVariables(c("seq_info", "seqids", "seqs", "hsps", "."))

NULL
