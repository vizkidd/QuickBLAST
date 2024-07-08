# QuickBLAST v1.2

Current BUILD is being tested and is not guaranteed to work. [Binaries of older version available here](https://github.com/vizkidd/QuickBLAST/releases/tag/binaries)


```R
??QuickBLAST
```

 Written in C++ and interfaced with R using Rcpp, the package is wrapped around ncbi-c++ toolkit's CBl2Seq Class (same with arrow) and exposing the functions to R with C linkage. I use getlogin() to store username in output metadata, this might raise red flags (in ArrowWrapper.cpp). QuickBLAST provides better interoperability with R for NCBI-BLAST. After much poking around, dependent libraries (Apache Arrow and NCBI-C++ Toolkit) are now compiled from scratch (and without Windows APIs on Windows - using MSYS2 and MinGW provided with RTools4.3).

The main difference between this PKG and the rest would be that
+ Quick blast is multi-threaded with { file reading (as chunks), BLASTing, wrapping hits into Arrow data structures }, and { writing of Arrow::RecordBatches to the output file in batches } is done in seperate threads. Hits are also converted into Rcpp::List if you want values to be returned to R.
+ QuickBLAST does not use Sys.Calls to invoke BLAST. You don't need BLAST programs in you system
+ BLAST DBs are not explicitly created

Cons :
+ Limited score attributes
 
Let me know if you want more information and please address bugs to me on github.

## Installation (under construction)

+ [Install Apache Arrow](https://arrow.apache.org/install/) - (Optional).
+ [Install RTools4.3](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) - (For Windows)

```R
devtools::install_github("https://github.com/vizkidd/QuickBLAST", force=T)
```

## Usage
```R
remotes::install_local("QuickBLAST_1.0_R_x86_64-pc-linux-gnu.tar.gz", build=F)
tblastx_ptr <- QuickBLAST::CreateNewBLASTInstance(seq_info = list(0,0,F), program = "tblastx", options = list("evalue"=1e-05, "pident"=0.75, "qcovhsp_perc"=0.75))
blastn_ptr <- QuickBLAST::CreateNewBLASTInstance(seq_info = list(0,0,F), program = "blastn", options = "")
 QuickBLAST::BLAST2Files(ptr=tblastx_ptr, query="ungrouped.cds", subject="ungrouped.cds", out_file="out.tmp", seq_limit=1000, show_progress=T,return_values=F, num_threads=5)
QuickBLAST::BLAST2Seqs(ptr=tblastx_ptr, query="AAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGG", subject="TTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
QuickBLAST::BLAST2Seqs(ptr=blastn_ptr, query=">i11166\nTGGCACGTCTGGTAGCAGTTTGCAGGGAAGGGGAAGAGGAATACCCGTTTCTCGCCAGACAGATCC", subject=">i11167\nATGGCACGTCTGGTAGCAGTTTGCAGGGAAGGGGAAGAGGAATACCCGTTTCTCGCCAGACAGATCCCCCTCTTCATCGATGACACTCTCACGATGGTGATGGAGTTTTCCGATAGCGTCATGG")
QuickBLAST::BLAST1Folder(ptr = tblastx_ptr, input_folder="test", extension= ".cds", out_folder="test_out", num_threads=7, reciprocal_hits=F)
QuickBLAST::BLAST2Folders(ptr=blastn_ptr, query="query", subject="subject", extension = ".cds", out_folder="test2_out", num_threads=8, reciprocal_hits=F)

```

<a name="blast_options"/>

#### QuickBLAST Options 
    
   Same as BLAST but DB & OUTPUT Format are not available. List of available options can be checked with `QuickBLAST::GetAvailableBLASTOptions()` (Empty elements from the list are removed and BLAST defaults are set on the c++ side). Inputs and Outputs are provided as parameters and sequence specification(strand, sequence type) can be provided during QuickBLAST object creation with `QuickBLAST::GetQuickBLASTInstance()` (or use the QuickBLAST::BLAST*() functions in R). Enums used by QuickBLAST in C++ are not exposed in R and only integers are used, check `QuickBLAST::GetQuickBLASTEnums()`.

#### BLAST Scores :

[Currently supported scores](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/include/objects/seqalign/Seq_align.hpp#0128)

#### Future : (Looking for suggestions)
+ Implement more scores and filtering options 
+ Include function for reading the arrow output files
+ Convert from arrow to GRanges (maybe with the use of arrow::Visit() functions)

Disclaimers for disclaimers, legal stuff for legal stuff and respect for respect, wherever it should go.

[LinkedIN](https://www.linkedin.com/in/vishveshkarthik/)
