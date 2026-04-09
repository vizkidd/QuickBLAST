#!/bin/bash

rm -rf BUILD/ inst/libs/

Rscript -e "Rcpp::compileAttributes(verbose = T)"

## Run this in your R console BEFORE devtools::document()
Rscript -e "devtools::install(); Sys.setenv(LD_LIBRARY_PATH = paste('$HOME/R/packages/QuickBLAST/libs', Sys.getenv('LD_LIBRARY_PATH'), sep = ':')); usethis::use_logo('logo.png'), devtools::document(quiet = T); devtools::build_site(quiet = F); pkgdown::build_news(); pkgdown::build_favicons(); " #roxygen2::roxygenize(clean = T)

# rm -f src/RcppExports.cpp R/RcppExports.R
# 
# Rscript -e "Rcpp::compileAttributes(verbose = T)"

rm -rf BUILD/ inst/libs/

Rscript -e "devtools::build_vignettes()"

Rscript -e "devtools::build(binary = T,vignettes = T)"

Rscript -e "devtools::build(vignettes = T)"

rm -rf BUILD/ inst/libs/

Rscript -e "devtools::install(force = T)"

R CMD build .

Rscript -e 'library(QuickBLAST); loadeddlls <- getLoadedDLLs(); getDLLRegisteredRoutines(loadeddlls$libQuickBLASTcpp)'



