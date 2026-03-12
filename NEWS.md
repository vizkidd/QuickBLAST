# News

All notable changes to this project will be documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project follows [Semantic Versioning](https://semver.org/).

---

## [Version 1.4.7] - 2026-03-13
### Updated:
- Considerable optimizations
- Attempt to set $ORIGIN through linker flags failed for ncbi-toolkit
- Tweaking windows build system for r-cran winbuilder
- Removing '-fdiagnostics-color=always'
- Clear verbosity
- shared lib trimming failed but stripping works
### Added:
- Using boost headers for quick vector access in RecordBatch to R List conversion
- auto-joining safe_jthread implementation wherever necessary
- Support for genomes thorugh sequence chunking
- Seperate stuct to simplify arrow container storage 

---

## [1.3.2] - 2026-03-08

### Bug
- There is a known bug with RStudio due to libc++ library collisions between QuickBLAST and RStudio

---
