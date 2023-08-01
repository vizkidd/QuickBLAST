# QuickBLAST

Requires NCBI BLAST headers and libs to compile from source. [Binaries available here](https://github.com/vizkidd/QuickBLAST/releases/tag/binaries)


```R
??QuickBLAST
```

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
    Same as BLAST but DB & OUTPUT Format are not available. List of available options can be checked with `QuickBLAST::GetAvailableBLASTOptions()` (Emplty elements from the list are removed and BLAST defaults are set on the c++ side). Inputs and Outputs are provided as parameters and sequence specification(strand, sequence type) can be provided during QuickBLAST object creation with `COMPLETE::GetQuickBLASTInstance()` (or use the BLAST2*() functions in R). Enums used by QuickBLAST in C++ are not exposed in R and only integers are used, check `COMPLETE::GetQuickBLASTEnums()`.
