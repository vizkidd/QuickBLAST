% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BLAST2Files}
\alias{BLAST2Files}
\title{BLAST 2 Files}
\usage{
BLAST2Files(
  ptr,
  query,
  subject,
  out_file,
  seq_limit,
  num_threads,
  show_progress,
  return_values,
  min_batch_size = 1024L
)
}
\arguments{
\item{ptr}{(QuickBLAST ptr) Pointer to QuickBLAST object created with CreateNewBLASTInstance().}

\item{query}{(string) Query FASTA File}

\item{subject}{(string) Subject FASTA File}

\item{out_file}{(string) Output Filename (Arrow Feather/IPC Format)}

\item{seq_limit}{(int) Batch Size to BLAST at a time. { -1 = Whole File, 0 - One sequence at a time or > 0 }}

\item{num_threads}{(int) Number of Threads}

\item{show_progress}{(bool) TRUE - Show progress, Set FALSE for multiple instances}

\item{return_values}{(bool) TRUE - Returns values as list, Default - FALSE - Does not return values (Return true on completion)}

\item{min_batch_size}{(int) Minimum Batch Size to start writing to output file. (Default - 1024)}
}
\value{
RecordBatchVector of Asynchronous BLAST Hits (OR) TRUE on successful completion
}
\description{
BLAST 2 Files and return the RecordBatchVector of BLAST Hits Asynchronously
}
\note{
Keep return_values_ = FALSE for large files to avoid choking R.
}
