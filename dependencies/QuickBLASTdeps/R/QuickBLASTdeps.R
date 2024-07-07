
#' Stub function that always returns true. Only to test the connection of DLLs and C function calls
#' @examples
#' \dontrun{
#' QuickBLASTdeps::isQuickBLASTLoaded()
#' }
#'
#' @return Always TRUE
#' @md
#' @export
isQuickBLASTLoaded <- function() {
  load_result <- .Call("isQuickBLASTLoaded")
  return(load_result)
}

R_dll_paths <- list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=".dll", full.names = T)

dll_paths <- c(             
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_core.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_general.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_pub.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_seq.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_trackmgr.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_eutils.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_misc.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libsqlitewrapp.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"liblmdb.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libefetch.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_seqext.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_xreader.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_xreader_id1.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_xreader_id2.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_xreader_cache.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libxxconnect2.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libpsg_client.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_xloader_genbank.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_web.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_align_format.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libutrtprof.dll.a"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libncbi_algo.dll.a"),
  #"inst/libs", Sys.getenv("R_ARCH"),"libarrow_msys2.dll.a",
  #"inst/libs", Sys.getenv("R_ARCH"),"msys-arrow-1601.dll",
  list.files(fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH")),pattern = "*arrow.*dll", full.names = T),
  R_dll_paths,
  # file.path(Sys.getenv("RTOOLS43_HOME"),"usr","bin","msys-gomp-1.dll"),
  # file.path(Sys.getenv("RTOOLS43_HOME"),"usr","bin","msys-stdc++-6.dll"),
  # file.path(Sys.getenv("RTOOLS43_HOME"),"usr","bin","msys-gcc_s-seh-1.dll"),
  # file.path(Sys.getenv("RTOOLS43_HOME"),"usr","bin","msys-2.0.dll")
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"msys-gomp-1.dll"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"msys-stdc++-6.dll"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"msys-gcc_s-seh-1.dll"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"msys-2.0.dll"),
  fs::path_package("QuickBLASTdeps","libs", Sys.getenv("R_ARCH"),"libQuickBLASTcpp.dll.a"),
  #"inst/libs", Sys.getenv("R_ARCH"),"QuickBLAST.dll"
)

.onLoad <- function(libname, pkgname) {
  # Load the DLLs when the package is loaded
  #library(arrow)
  
  for (dll_path in dll_paths) {
    if (!file.exists(dll_path)) {
      cat("DLL file not found:", dll_path, "\n")
    } else {
      # dyn.load(dll_path, local=F, now = T)
      if(!invisible(is.loaded(dll_path))){
        dyn.load(dll_path,now = T)
      }
      cat("Loaded DLL:", dll_path, "\n")
    }
  }
  #rdyncall::.dynload(dll_paths)
  # dyn.load("inst/libs", Sys.getenv("R_ARCH"),"QuickBLAST.dll", local=F, now = T)
  #dyn.load("inst/libs", Sys.getenv("R_ARCH"),"QuickBLAST.dll",now = T)
}

.onAttach <- function(libname, pkgname) {
  # Load the DLLs when the package is loaded
  #library(arrow)
  
  for (dll_path in dll_paths) {
    if (!file.exists(dll_path)) {
      cat("DLL file not found:", dll_path, "\n")
    } else {
      #dyn.load(dll_path, local=F, now = T)
      if(!invisible(is.loaded(dll_path))){
        dyn.load(dll_path, now = T)
      }
      cat("Loaded DLL:", dll_path, "\n")
    }
  }#
  #rdyncall::.dynload(dll_paths)
  #dyn.load("inst/libs", Sys.getenv("R_ARCH"),"QuickBLASTcpp.dll", local=F, now = T)
}

# .onUnload() function
.onUnload <- function(libpath) {
  # Unload the DLLs when the package is unloaded
  for (dll_path in dll_paths) {
    # if(is.loaded(dll_path)){
      if (dyn.unload(dll_path)) {
        cat("Unloaded DLL:", dll_path, "\n")
      } else {
        cat("Failed to unload DLL:", dll_path, "\n")
      }
    # }
  }
}