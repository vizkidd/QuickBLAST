# QuickBLAST

Requires [NCBI BLAST C++](https://github.com/ncbi/ncbi-cxx-toolkit-public) includes and libs to compile from source. [Binaries available here](https://github.com/vizkidd/QuickBLAST/releases/tag/binaries)


```R
??QuickBLAST
```

 It is written in C++ and interfaced with R using Rcpp. Compiling it from source requires ncbi headers and libs so it's better to use the compiled binary. I am just wrapping around ncbi-c++ toolkits CBl2Seq Class with my own class (same with arrow). Not really doing much in the R side except for exposing functions. I use getlogin() to store username as output metadata, this might raise red flags (in arrow_wrapper-functions.hp). Quick BLAST is a bit faster on soft Benchmarking.  Ncbi-BLAST is in a bit of a fix because you have to compile the whole BLAST suite repo to use it in R and they have their own "ecosystem" so exposing their classes in R would be a bit of a hassle, that is why I just wrapped parts of the BLAST suite. Essentially, the huge size of the [C++ toolkit](https://github.com/ncbi/ncbi-cxx-toolkit-public) makes it infeasible (or just annoying) to include the entire toolkit in this repo, hence the hard coding of ncbi-tools++ headers (which are usually installed in /usr/local/). Moreover CRAN would not accept packages which include other libraries (unless compiled from source) hence I had to provide the binary version. Let me know if you want more information and please address bugs to me on github.

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

