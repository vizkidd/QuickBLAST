![QuickBLAST Logo](man/figures/logo.png)

# QuickBLAST v1.6.5

R is widely used for data analysis, but running NCBI's standard BLAST tools within R has traditionally been slow. Because the NCBI C++ toolkit is massive and inflexible, existing R packages are forced to run BLAST as an external subprocess, which creates major read/write bottlenecks.

QuickBLAST solves this by building a direct bridge between R and the NCBI C++ toolkit via Rcpp. By bypassing traditional text-based formatting and transporting data directly into memory using Apache Arrow, QuickBLAST performs sequence comparisons exceptionally fast.

## Key Features

-   Zero Subprocesses: Runs entirely natively within your R session. QuickBLAST completely avoids Sys.Call\(\) and does not require pre-installed BLAST executables.
-   True Multi-Threading: Employs a concurrent architecture where file reading (in chunks), sequence alignment, Arrow wrapping, and disk writing all occur simultaneously in separate threads.
-   Memory & I/O Efficiency: Wraps hits natively into Arrow data structures (Arrow::RecordBatches) for large-scale disk writing, or returns an Rcpp::List directly to R for smaller queries.
-   No Length Limits: Removes legacy limits on sequence and header lengths.
-   Versatile: Instantly compare raw sequences, local FASTA files, local databases, or remote NCBI databases.

## Requires

-   GNU GCC \>= 13.3.0 with C++20 support
-   CMake
-   OpenMP support (-fopenmp)
-   R \> 4.4.0
-   Rtools \>= 4.4 (Windows)
-   `sudo apt install libsqlite3-dev libeigen3-dev libboost-dev libfontconfig1-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev cmake` (Linux)

Written in C++ and interfaced with R using Rcpp, the package is wrapped directly around the NCBI-C++ Toolkit's BLAST-specific classes and Apache Arrow, exposing these functions to R with C linkage.

The main difference between this package and legacy wrappers is the data lifecycle. Instead of waiting for a sequence alignment to finish before parsing a TSV, QuickBLAST sets up a sophisticated pipeline:

1)  Producer Threads: Read sequence files in chunks and perform the mathematical sequence comparisons.
2)  Transformer Threads: Immediately wrap the alignments into Arrow data structures in memory.
3)  Consumer Threads: Batch write Arrow::RecordBatches directly to an output file. (Because these operate independently, your CPU is fully utilized without I/O blocking).

## Installation

-   [Install RTools 4.4 or greater](https://cran.r-project.org/bin/windows/Rtools/) - (For Windows. Rtools must be the same Major and Minor version of R)

``` r
devtools::install_github("https://github.com/vizkidd/QuickBLAST", force=T)
```

<a name="blast_options"/>

#### QuickBLAST Options

List of available options can be checked with `QuickBLAST::GetAvailableBLASTOptions()`. Enums used by QuickBLAST in C++ are not exposed in R and only integers are used, check `QuickBLAST::GetQuickBLASTEnums()`.

#### Output Formats

-   [arrow::ipc](https://arrow.apache.org/docs/format/Columnar.html#serialization-and-interprocess-communication-ipc)
-   arrow::csv
-   [arrow::parquet](https://parquet.apache.org/docs/file-format/)

``` r
?QuickBLAST::LoadBLASTHits
```

#### BLAST Scores :

[Currently supported scores](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/include/objects/seqalign/Seq_align.hpp#0128)

## Usage

``` r
??QuickBLAST
```

## Quick Start Guide

#### 1) Initialize QuickBLAST

QuickBLAST uses "instances" to maintain search parameters (like E-values and programs) in the background.

``` r
library(QuickBLAST)

# Create a Nucleotide (blastn) instance
blastn_inst <- QuickBLAST::CreateQuickBLASTInstance(
  seq_type = 0, strand = 0, program = "blastn", options = "-evalue 100000"
)

# Create a Protein (blastp) instance
blastp_inst <- QuickBLAST::CreateQuickBLASTInstance(
  seq_type = 1, strand = 0, program = "blastp", save_sequences = FALSE, save_hsp_sequences = TRUE
)
```

#### 2) Compare Raw Sequences

You can pass raw character strings directly to QuickBLAST without needing to write temporary FASTA files to your disk. Results are returned natively as an Rcpp::List.

``` r
QuickBLAST::BLAST2Seqs(
  blastn_inst, 
  query = "AAAAAAAAAAAATTTTTTTTTTTTGGGGGGGGGGGCCCCCCCCC", 
  subject = "TTTTTTTTTTTGGGGGGGGGGGG"
)
```

#### 3) File and Database Comparisons

QuickBLAST makes large-scale genomics easy with built-in file and database tools.

##### BLASTing Two Files:

``` r
QuickBLAST::BLAST2Files(blastn_inst, query = "query.fasta", subject = "genome.fasta")
```

##### Creating and Searching a Local Database:

``` r
# 1. Compile the database
QuickBLAST::MakeBLASTDB(
  in_seq = "reference_genome.fasta", 
  db_type = "nucl", 
  out_db = "my_custom_db"
)

# 2. Search against it
QuickBLAST::BLAST2DBs(blastn_inst, query = "query.fasta", db = "my_custom_db")
```

#### 4) Remote NCBI Searching

If you don't want to download databases, you can query NCBI's remote servers directly from R:

``` r
QuickBLAST::RemoteBLAST(
  blastp_inst, 
  query_input="MQILLVEDDNTLFQELKKELEQWDFNVAGIEDFG...", 
  database= "pdb", 
  input_type=1, 
  return_values=TRUE
)
```

#### 5) Instance Management

Because QuickBLAST opens direct connections to C++ libraries, it includes utility functions to track and clean up memory.

``` r
# See how many instances are running
QuickBLAST::GetInstanceCount()

# Delete a specific instance by its ID
QuickBLAST::DeleteQuickBLASTInstance(1)
```

Inherits and follows the licenses of Apache Arrow and NCBI-C++-Toolkit.
Parts of the code, optimizations and documentation in the recent versions of QuickBLAST were written with the help of Google Gemini AI.
Developed and maintained by [vizkidd](https://github.com/vizkidd/QuickBLAST).
