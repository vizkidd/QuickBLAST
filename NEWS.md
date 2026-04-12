# News

All notable changes to this project will be documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project follows [Semantic Versioning](https://semver.org/).

---

# Version 1.7.0 - 2026-04-12

## Updated:
- CXX STD to CXX20
- Package directory structure
- Startup message moved to .OnAttach()
- manpages and vignettes

## Added:
- Patch fix to arrow for timezone format with GCC > 14 and C++20
- QuickBLAST logo

## Removed:
- Non-portable flags: -fPIC -O3 -Wall
- pkgdown site files

---

# Version 1.6.5 - 2026-04-09

## Updated:
- README
- vignettes
- R doc site

## Fixed:
- roxygen //' tags
- File existence check for file BLAST
- Race conditions and corruptions in file BLAST when the no. of threads was greater than the input sequence count
- When return_values==F and writing to output file, return of empty record batch is not erroreous.

---

# Version 1.6.4 - 2026-04-08

## Updated:
- Static dependencies + Shared R library as imposed by the R package policy - Shared libs are not built and copied over anymore
- C std=gnu2x, CXX std=gnu++20
- [Windows] Using source built sqlite3 due to missing aux functions
- [Windows] Toolchain injection provides 40+ patches to compile ncbi-cxx-toolkit with Rtools45 toolchain
- Seperated *nix and Windows build flow using preprocessor directives and build vars
- BLAST_dbs() parallelization
- [Windows] OpenMP Support

## Fixed:
- Windows build
- iconv for multi-UTF support & unbounded file mapping and finish of dangling file writer thread in arrow wrapper

## Added:
- [Windows] Dependencies provided as external projects to be compiled by CMake
- multi-os memory mapping, parquet and csv format support to arrow wrapper
- Chunking of large sequences to support genome scale megablast - blast does not work for some reason

## Removed:
- patchelf, R source compilation

---

# Version 1.5.0 - 2026-03-25

## Updated:
- Added optimal Visit() functions that uses the visitor pattern of Arrow to flatten arrow types to R lists

---

# Version 1.4.7 - 2026-03-13
## Updated:
- Considerable optimizations
- Attempt to set $ORIGIN through linker flags failed for ncbi-toolkit
- Tweaking windows build system for r-cran winbuilder
- Removing '-fdiagnostics-color=always'
- Clear verbosity
- shared lib trimming failed but stripping works

## Added:
- Using boost headers for quick vector access in RecordBatch to R List conversion
- auto-joining safe_jthread implementation wherever necessary
- Support for genomes thorugh sequence chunking
- Seperate stuct to simplify arrow container storage 

---

# Version 1.3.2 - 2026-03-08

## Bug
- There is a known bug with RStudio due to libc++ library collisions between QuickBLAST and RStudio

---
