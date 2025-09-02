#!/bin/bash

rm -rf BUILD/ inst/

# Rscript -e "Rcpp::compileAttributes(verbose = T)"

Rscript -e "devtools::document(quiet = T)" #roxygen2::roxygenize(clean = T)

# rm -f src/RcppExports.cpp R/RcppExports.R
# 
# Rscript -e "Rcpp::compileAttributes(verbose = T)"

rm -rf BUILD/ inst/

Rscript -e "devtools::build_vignettes()"

Rscript -e "devtools::build(binary = T,vignettes = T)"

rm -rf BUILD/ inst/

Rscript -e "devtools::install(force = T)"

Rscript -e 'library(QuickBLAST); loadeddlls <- getLoadedDLLs(); getDLLRegisteredRoutines(loadeddlls$libQuickBLASTcpp)'



