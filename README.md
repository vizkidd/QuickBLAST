# QuickBLAST

Requires NCBI BLAST headers and libs to compile from source. [Binaries available here](https://github.com/vizkidd/QuickBLAST/releases/tag/binaries)


```R
??QuickBLAST
```

```R
remotes::install_local("QuickBLAST_1.0_R_x86_64-pc-linux-gnu.tar.gz", build=F)
tblastx_ptr <- QuickBLAST::CreateNewBLASTInstance(seq_info = list(0,0,F), program = "tblastx", options = list("evalue"=1e-05, "pident"=0.75, "qcovhsp_perc"=0.75))
blastn_ptr <- QuickBLAST::CreateNewBLASTInstance(seq_info = list(0,0,F), program = "blastn", options = "")
 QuickBLAST::BLAST2Files(ptr=tblastx_ptr, query="ungrouped.cds", subject="ungrouped.cds", outFile="out.tmp", seq_limit=1000, show_progress=T,return_values=F, num_threads=8)
QuickBLAST::BLAST2Seqs(ptr=blastn_ptr, query="AAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC", subject="TTTTTTTTTTTTGGGGGGGGGGGGGGGG")
```

<a name="blast_options"/>

#### QuickBLAST Options 
    Same as BLAST but DB & OUTPUT Format are not available. List of available options can be checked with `QuickBLAST::GetAvailableBLASTOptions()` (Emplty elements from the list are removed and BLAST defaults are set on the c++ side). Inputs and Outputs are provided as parameters and sequence specification(strand, sequence type) can be provided during QuickBLAST object creation with `COMPLETE::GetQuickBLASTInstance()` (or use the BLAST2*() functions in R). Enums used by QuickBLAST in C++ are not exposed in R and only integers are used, check `COMPLETE::GetQuickBLASTEnums()`.
