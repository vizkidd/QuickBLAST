// [[Rcpp::plugins(openmp)]]

#include <stdlib.h>
#include <R_ext/Rdynload.h>
// #include <R_ext/Visibility.h>

#include <RcppCommon.h>
#include <Rcpp.h>

// [[Rcpp::export]]
RcppExport SEXP isQuickBLASTLoaded()
{
    Rcpp::Rcout << "QuickBLAST dependencies Loaded!" << std::endl;
    return Rcpp::wrap(true);
}

static const R_CallMethodDef callMethods[] = {
    {"isQuickBLASTLoaded", (DL_FUNC)&isQuickBLASTLoaded, 0},
    {NULL, NULL, 0}};

/*
// R_CMethodDef
// R_ExternalMethodDef
static const R_CMethodDef cMethods[] = {
    {"test_QBR", (DL_FUNC)&test_QBR, 0},
    {"test_QBR_cpp", (DL_FUNC)&test_QBR_cpp, 0},
    //{"test_QBcpp", (DL_FUNC)&test_QBcpp, 0},
    // {"GetQuickBLASTInstance", (DL_FUNC)&GetQuickBLASTInstance, 1},
    // {"GetInstanceCount", (DL_FUNC)&GetInstanceCount, 0},
    // {"GetInstanceID", (DL_FUNC)&GetInstanceID, 1},
    // {"CreateQuickBLASTInstance", (DL_FUNC)&QB_CreateNewBLASTInstance, 5},
    // {"SetQuickBLASTOptions", (DL_FUNC)&SetQuickBLASTOptions, 3},
    // {"BLAST2Seqs", (DL_FUNC)&QB_BLAST2Seqs, 3},
    // {"BLAST2Files", (DL_FUNC)&QB_BLAST2Files, 9},
    // {"BLAST2Folders", (DL_FUNC)&QB_BLAST2Folders, 8},
    // {"BLAST1Folder", (DL_FUNC)&QB_BLAST1Folder, 7},
    NULL};*/

void R_init_QuickBLASTdeps(DllInfo *dll)
{
    // R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
