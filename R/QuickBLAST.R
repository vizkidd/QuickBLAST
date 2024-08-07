# #' Get an Instance of QuickBLAST class and its exposed methods
# #' @note Check BLAST C++ Call in Help for the list of parameters for the exposed BLAST function. Exposed C++ function only takes BLAST options as string.
# #' @seealso [QuickBLAST::CreateNewBLASTInstance()], [QuickBLAST::GetQuickBLASTEnums()]
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
#' @export
GetQuickBLASTEnums <- function() {
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
#' @return A List of Available BLAST options
#' @export
GetAvailableBLASTOptions <- function() {
  warning("These are Defaults, CREATE a NEW LIST with ONLY the OPTIONS THAT YOU NEED (or) Use this function in blast_options to set BLAST defaults for the chosen BLAST program.")
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
#' Give the path to a BLAST Hits file to load it into a data.frame(BLAST HITs Table). The column names can be provided as col.names. Rows with NAs are automatically removed.
#'
#' Note : If use.feather is enabled, then the functions returns an arrow::RecordBatchStreamReader object. Call wrap it in an iterators::iter() (like batch_iter <- iter(function(){ batch_reader$read_next_batch() })) and call iterators::nextElem(batch_iter) to get each batch of the BLAST Hits.
#'
#'
#' @param infile BLAST hits filename (not a connection) (Gzipped files supported)
#' @param sep Delimiter of the BLAST File columns. Default - '\\t'
#' @param header Does the file have a header? . Default - FALSE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - FALSE
#' @return Data Frame with BLAST Results
#' @export
LoadBLASTHits <- function(infile, sep = "\t", header = F, use.feather = F) {
  if (try(file.exists(infile)) && file.info(infile)$size > 0) { # any(grepl(x = class(infile), pattern = "gzfile|connection", ignore.case = T)) #
    # if(gzipped){
    #  infile <- gzfile(description = infile, open = "r")
    # }
    # blast_results <- read.table(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, fill = t,na.strings = NA)
    if (!use.feather) {
      blast_results <- iterators::iread.table(file = infile, row.names = NULL, header = header, sep = sep, quote = "", blank.lines.skip = T, fill = T, na.strings = "NA") # data.table::fread(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, nThread = n_threads)
      return(blast_results)
    } else {
      arrow_lfs <- arrow::LocalFileSystem$create()
      arrow_i_stream <- arrow_lfs$OpenInputStream(infile)
      batch_reader <- arrow::RecordBatchStreamReader$create(arrow_i_stream)
      # blast_results <- arrow::read_feather(file = infile, mmap = T)
      return(batch_reader) # batch_iter <- iter(function())
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
#' \dontrun{
#' one2one("ungrouped.cds", "ungrouped.cds", seq_info = list(GetQuickBLASTEnums()$ESeqType$eNucleotide, GetQuickBLASTEnums()$EStrand$ePlus, F), file_ext = "fa", blast.sequence.limit = 200, input_type = GetQuickBLASTEnums()$EInputType$eFile, n_threads = 8, blast_program = "tblastx", output_dir = "./", blast_options = QuickBLAST::GetAvailableBLASTOptions(), verbose = T)
#' }
#'
#' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
#' @param first_list Vector of FASTA Filenames or Strings
#' @param second_list Vector of FASTA Filenames or Strings
#' @param seq_info Ordered List of 1) (int) Sequence Type, 2) (int) Strand (bool), 3) Save Sequences in BLAST Hits? : TRUE - BLAST Hits have sequences, FALSE - FASTA sequences are not stored in BLAST Hits. Check QuickBLAST::GetQuickBLASTEnums() for available enums
#' @param file_ext File extension of input files. eg- ".cds" or ".fa", Unused if input_type is GetQuickBLASTEnums()$EInputType$eSequencesString
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param output_dir Path to BLAST output
#' @param blast.sequence.limit Maximum number of sequences to be BLASTed at a time, not used for Seqs
#' @param n_threads Number of threads. Default - 8
#' @param blast_program BLAST Program to use
#' @param blast_options BLAST Options to use - QuickBLAST::GetAvailableBLASTOptions()
#' @param input_type (integer) 0 for Files, 1 for Sequences - QuickBLAST::GetQuickBLASTEnums()
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param input_type (integer) 0 for Files, 1 for Sequences - QuickBLAST::GetQuickBLASTEnums()$EInputType
#' @param verbose Print DEBUG Messages?
#' @md
#' @export
one2one <- function(first_list, second_list, seq_info, file_ext = ".fa", input_prefix_path = NULL, blast.sequence.limit = 1000, input_type, n_threads, blast_program, output_dir = "./", blast_options, verbose = T) {
  if (!is.null(input_prefix_path)) {
    first_list <- paste(input_prefix_path, "/", first_list, file_ext, sep = "")
    second_list <- paste(input_prefix_path, "/", second_list, file_ext, sep = "")
  }else{
    first_list <- paste(first_list, file_ext, sep = "")
    second_list <- paste(second_list, file_ext, sep = "")
  }

  list_combos <- unique(tidyr::crossing(first_list[order(first_list)], second_list[order(second_list)]))

  # return_data <- furrr::future_map2(.x=first_list[order(first_list)], .y=second_list[order(second_list)], .f=function(x,y){
  parallel::mclapply(seq_along(1:nrow(list_combos)), function(idx) {
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

    blast_ptr <- QuickBLAST::CreateNewBLASTInstance(seq_info = seq_info, program = blast_program, options = blast_options)

    switch(input_type,
      { # eFile
        return(QuickBLAST::BLAST2Files(ptr = blast_ptr, query = x, subject = y, out_file = paste(output_dir, "/", basename(tools::file_path_sans_ext(x)), ".", basename(tools::file_path_sans_ext(y)), ".hits", sep = ""), seq_limit = blast.sequence.limit, show_progress = T, return_values = F, min_batch_size = 256, num_threads = n_threads))
      },
      { # eSequenceString
        return(QuickBLAST::BLAST2Seqs(ptr = blast_ptr, query = x, subject = y))
      },
      { # eFolder
      }
    )
  }, mc.cores = n_threads, mc.silent = !verbose)

  # return(return_data)
}


#' Execute all2all QuickBLAST
#'
#' Executes All-to-All QuickBLAST between two lists of organisms/genes/clusters. Output BLAST files are bi-directional and are stored in the filename filename1.filename2.all2all under output_dir. (All-to-All is simply Many-to-Many association)
#'
#' @examples
#' \dontrun{
#' all2all("ungrouped.cds", "ungrouped.cds", seq_info = list(GetQuickBLASTEnums()$ESeqType$eNucleotide, GetQuickBLASTEnums()$EStrand$ePlus, F), file_ext = "fa", blast.sequence.limit = 200, input_type = GetQuickBLASTEnums()$EInputType$eFile, n_threads = 8, blast_program = "tblastx", output_dir = "./", blast_options = QuickBLAST::GetAvailableBLASTOptions(), verbose = T)
#' }
#'
#' @seealso [QuickBLAST::GetAvailableBLASTOptions()], [QuickBLAST::GetQuickBLASTEnums()]
#' @param first_list Vector of FASTA Filenames or Strings
#' @param second_list Vector of FASTA Filenames or Strings
#' @param seq_info Ordered List of 1) (int) Sequence Type, 2) (int) Strand (bool), 3) Save Sequences in BLAST Hits? : TRUE - BLAST Hits have sequences, FALSE - FASTA sequences are not stored in BLAST Hits. Check QuickBLAST::GetQuickBLASTEnums() for available enums
#' @param file_ext (Optional) File extension of input files. eg- ".cds" or ".fa", Unused if input_type is GetQuickBLASTEnums()$EInputType$eSequencesString
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param output_dir Path to BLAST output
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param input_type (integer) 0 for Files, 1 for Sequences - QuickBLAST::GetQuickBLASTEnums()$EInputType
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param blast.sequence.limit Maximum number of sequences allowed in each BLAST file. Default - 0. If the query FASTA sequences > blast.sequence.limit, the sequences are split into multiple files and BLASTed. Use 0 to not split & copy the files into temporary files
#' @param n_threads Number of Threads
#' @param verbose Print Debug Messages?
#' @md
#' @export
all2all <- function(first_list, second_list, input_type, seq_info, blast_program, file_ext = ".fa", blast_options = "", output_dir = "./", input_prefix_path = NULL, blast.sequence.limit = 1000, n_threads = 8, verbose = T) {
  if (verbose) {
    cat(paste("All2All QuickBLAST Started...", "\n", sep = ""))
    print(paste(first_list, collapse = ","))
    print(paste(second_list, collapse = ","))
  }

  dir.create(path = output_dir, recursive = T, showWarnings = F)

  list_combinations <- unique(tidyr::crossing(first_list, second_list))
  # furrr::future_map2(.x=list_combinations$first_list, .y=list_combinations$second_list, .f=function(first_set,second_set){
  parallel::mclapply(seq_along(1:nrow(list_combinations)), function(idx) {
    first_set <- list_combinations[idx, 1]
    second_set <- list_combinations[idx, 2]

    tryCatch(QuickBLAST::one2one(first_list = first_set, second_list = second_set, seq_info = seq_info, file_ext = file_ext, input_prefix_path = input_prefix_path, blast.sequence.limit = blast.sequence.limit, input_type = input_type, n_threads = n_threads, blast_program = blast_program, output_dir = output_dir, blast_options = blast_options, verbose = verbose),
      error = function(cond) {
        if (verbose) {
          message(cond)
        }
      }
    )
    tryCatch(QuickBLAST::one2one(first_set, second_set, seq_info = seq_info, file_ext = file_ext, input_prefix_path = input_prefix_path, blast.sequence.limit = blast.sequence.limit, input_type = input_type, n_threads = n_threads, blast_program = blast_program, output_dir = output_dir, blast_options = blast_options, verbose = verbose),
      error = function(cond) {
        if (verbose) {
          message(cond)
        }
      }
    )
  }, mc.cores = n_threads, mc.preschedule = T)
}

# R_dll_paths <- c(
#   # list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=".dll", full.names = T),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Riconv.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"R.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rgraphapp.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rblas.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"R.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rlapack.dll")
#   # fs::path_package("Rcpp","libs", Sys.getenv("R_ARCH"),"Rcpp.dll")
#   # list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=".dll", full.names = T)
#   # fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rcpp.dll")
#   # fs::path_package("Rcpp","libs", Sys.getenv("R_ARCH"),"Rcpp.dll")
# )

# msys_dlls <- c("msys-bz2-1", "msys-lzo2-2", "msys-z", "msys-pcre-1", "msys-zstd-1", "msys-sqlite3-0", "msys-gomp-1", "msys-stdc++-6", "msys-gcc_s-seh-1", "msys-2.0")

# msys_dll_paths <- c( fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"), paste(msys_dlls, .Platform$dynlib.ext,sep="")) )

# cat(list.files(fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH")),pattern = "*arrow.*dll", full.names = F))
#"Riconv","R","Rgraphapp","Rblas","R","Rlapack"
#"Riconv", "R", "Rgraphapp", "Rblas", "R", "Rlapack", 
R_dlls <- c("Riconv", "R", "Rgraphapp", "Rblas", "R", "Rlapack")
R_dll_paths <- paste0(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),.Platform$file.sep, R_dlls, .Platform$dynlib.ext)

dlls <- c("libncbi_core", "libncbi_general", "libncbi_pub", "libncbi_seq", "libncbi_trackmgr", "libncbi_eutils", "libncbi_misc", "libsqlitewrapp", "liblmdb", "libefetch", "libncbi_seqext", "libncbi_xreader", "libncbi_xreader_id1", "libncbi_xreader_id2", "libncbi_xreader_cache", "libxxconnect2", "libpsg_client", "libncbi_xloader_genbank", "libncbi_xloader_blastdb", "libncbi_xloader_blastdb_rmt", "libncbi_general", "libncbi_web", "libncbi_align_format", "libutrtprof", "libncbi_algo", tools::file_path_sans_ext(list.files(fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH")),pattern = "*arrow.*dll", full.names = F)), "libncbi_blastinput", "libQuickBLASTcpp")
# dlls <- c("libncbi_core", "libncbi_general", "libncbi_pub", "libncbi_seq", "libncbi_trackmgr", "libncbi_eutils", "libncbi_misc", "libsqlitewrapp", "liblmdb", "libefetch", "libncbi_seqext", "libncbi_xreader", "libncbi_xreader_id1", "libncbi_xreader_id2", "libncbi_xreader_cache", "libxxconnect2", "libpsg_client", "libncbi_xloader_genbank", "libncbi_web", "libncbi_align_format", "libutrtprof", "libncbi_algo", "libQuickBLASTcpp")


dll_paths <- paste(fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH")), .Platform$file.sep, dlls, .Platform$dynlib.ext,sep="") 

dll_obj_list <-  list()

.onLoad <- function(libname, pkgname) {
  # # Load the DLLs when the package is loaded
  # require(QuickBLAST)
  #require(QuickBLASTdeps)
  # require(arrow)
  # require(Rcpp)
  # require(RcppProgress)
  #R_dll_paths, msys_dll_paths
  # library.dynam("Rcpp", "Rcpp", fs::path_package("Rcpp","..",".."))
  for (dll_path in c(R_dll_paths)) {
    dyn.load(dll_path,now = T)
  }
  for (dll_path in c(dll_paths)) {
  # for (dll_path in c(dlls)) {
    if (!file.exists(dll_path)) {
      cat("DLL file not found:", dll_path, "\n")
    } else {
      # dyn.load(dll_path, local=F, now = T)
      # library.dynam(dll_path, "QuickBLAST", fs::path_package("QuickBLAST",".."))
      # if(!invisible(is.loaded(dll_path))){
        dll_obj_list <- append(dll_obj_list, list(dyn.load(dll_path,now = T)))
      # }
      cat("Loaded DLL:", dll_path, "\n")
    }
  }
}

.onAttach <- function(libname, pkgname) {
  Sys.setenv("ARROW_DEFAULT_MEMORY_POOL"="system")
  Sys.setenv("ARROW_DEBUG_MEMORY_POOL"="warn")
  Sys.setenv("OMP_NUM_THREADS"=parallel::detectCores(all.tests = T, logical = T))
  Sys.setenv("OMP_DYNAMIC"="TRUE")
  Sys.setenv("OMP_WAIT_POLICY"="PASSIVE")
  Sys.setenv("OMP_DISPLAY_ENV"="VERBOSE")
  packageStartupMessage("QuickBLAST Loaded!")
  packageStartupMessage("Version: ", utils::packageVersion("QuickBLAST"))
}

# .onUnload() function
.onUnload <- function(libpath) {
  # # Unload the DLLs when the package is unloaded
  # detach("package:QuickBLASTdeps", unload = TRUE)
  #msys_dll_paths
  packageStartupMessage("Unloading QuickBLAST...")

  loaded_dlls <- getLoadedDLLs()
  loaded_dlls <- loaded_dlls[na.omit(match(dlls,names(loaded_dlls)))]

  for(dll_info in loaded_dlls){
    dyn.unload(dll_info[["path"]])
  }

  # for (dll_path in c(rev(c(dll_paths)))) {
  #   if(is.loaded(dll_path)){
  #     if (dyn.unload(dll_path)) {
  #       packageStartupMessage(cat("Unloaded DLL:", dll_path, "\n"))
  #     } else {
  #       packageStartupMessage(cat("Failed to unload DLL:", dll_path, "\n"))
  #     }
  #   }
  # }
  # detach("package:QuickBLAST", unload = TRUE)
}

.onDetach <- function(libpath) {
  packageStartupMessage("Detaching QuickBLAST...")
  # rm(list = ls(envir = .GlobalEnv, pattern = "^QuickBLAST"))
  
  # for (dll_path in c(rev(c(dll_paths)))) {
  #   if(is.loaded(dll_path)){
  #     if (dyn.unload(dll_path)) {
  #       packageStartupMessage(cat("Unloaded DLL:", dll_path, "\n"))
  #     } else {
  #       packageStartupMessage(cat("Failed to unload DLL:", dll_path, "\n"))
  #     }
  #   }
  # }
}
