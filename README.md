# QuickBLAST v1.2.4

Current BUILD is being tested on linux and is not guaranteed to work on Windows. [Binaries of older version available here](https://github.com/vizkidd/QuickBLAST/releases/tag/binaries)

## Requires

-   GNU GCC >= 13.3.0
-   CMake
-   OpenMP support (-fopenmp)
-   R \> 4.4.0
-   Rtools \>= 4.4 (Windows)
-   `sudo apt install libsqlite3-dev libeigen3-dev libboost-dev libfontconfig1-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev cmake` (Linux)

``` r
??QuickBLAST
```

Written in C++ and interfaced with R using Rcpp, the package is wrapped around ncbi-c++ toolkit's CBl2Seq Class (same with arrow) and exposing the functions to R with C linkage. I use getlogin() to store username in output metadata, this might raise red flags (in ArrowWrapper.cpp). QuickBLAST provides better interoperability with R for NCBI-BLAST. After much poking around, dependent libraries (Apache Arrow and NCBI-C++ Toolkit) are now compiled from scratch (and without Windows APIs on Windows - using MSYS2 and MinGW provided with RTools4.4).

The main difference between this PKG and the rest would be that + Quick blast is multi-threaded with { file reading (as chunks), BLASTing, wrapping hits into Arrow data structures }, and { writing of Arrow::RecordBatches to the output file in batches } is done in seperate threads. Hits are also converted into Rcpp::List if you want values to be returned to R. + QuickBLAST does not use Sys.Calls to invoke BLAST exes. You don't need BLAST programs in you system + BLAST DBs are not explicitly created

~~Cons : + Limited score attributes~~

Let me know if you want more information and please address bugs to me on github.

## Installation (under construction)

-   [Install Apache Arrow](https://arrow.apache.org/install/) - (Optional).
-   [Install RTools 4.4 or greater](https://cran.r-project.org/bin/windows/Rtools/) - (For Windows. Rtools must be the same Major and Minor version of R)

``` r
devtools::install_github("https://github.com/vizkidd/QuickBLAST", force=T)
```

## Usage

``` r

```

<a name="blast_options"/>

#### QuickBLAST Options

List of available options can be checked with `QuickBLAST::GetAvailableBLASTOptions()` (Empty elements from the list are removed and BLAST defaults are set on the c++ side). Inputs and Outputs are provided as parameters and sequence specification(strand, sequence type) can be provided during QuickBLAST object creation with `QuickBLAST::GetQuickBLASTInstance()` (or use the QuickBLAST::BLAST\*() functions in R). Enums used by QuickBLAST in C++ are not exposed in R and only integers are used, check `QuickBLAST::GetQuickBLASTEnums()`.

#### Output Formats

+ [arrow::ipc](https://arrow.apache.org/docs/format/Columnar.html#serialization-and-interprocess-communication-ipc)
+ [~~arrow::csv~~]()
+ [arrow::parquet](https://parquet.apache.org/docs/file-format/)

``` r
?QuickBLAST::LoadBLASTHits
```

#### BLAST Scores :

[Currently supported scores](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/include/objects/seqalign/Seq_align.hpp#0128)

#### Future : (Looking for suggestions)

-   Implement more scores and filtering options
-   Include function for reading the arrow output files
-   Convert from arrow to GRanges (maybe with the use of arrow::Visit() functions)

Disclaimers for disclaimers, legal stuff for legal stuff and respect for respect, wherever it should go.

