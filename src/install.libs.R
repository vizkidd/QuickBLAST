# Determine the correct destination directory (usually inst/libs/x64)
libarch <- if (nzchar(R_ARCH)) paste0("libs", R_ARCH) else "libs"
dest <- file.path(R_PACKAGE_DIR, libarch)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

# if(grepl(Sys.getenv("R_PLATFORM"), pattern ="linux", ignore.case = T)){
#   dest <- file.path(R_PACKAGE_DIR, paste0('libs'))
# }else{
#   dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH)) 
# }
# dir.create(dest, recursive = TRUE, showWarnings = FALSE)

## staged_dest <- file.path(R_PACKAGE_DIR,"..","..","..","QuickBLAST", paste0('libs', R_ARCH))
## dir.create(staged_dest, recursive = TRUE, showWarnings = FALSE)

##C:/Users/vishv/AppData/Local/Temp/RtmpOSDZyI/temp_libpath1bacef43ec2/00LOCK-QuickBLAST_cmake_test/00new/QuickBLAST
##  C:/Users/vishv/AppData/Local/Temp/RtmpO6DZ4p/temp_libpath2fb068301991/QuickBLAST/ 

cat(paste("INSTALLING....", R_PACKAGE_NAME, " FOR R - ", R_ARCH, "\n", sep=""))

# cat(list.files(path = fs::path_package("QuickBLAST"), all.files = T, full.names = T, recursive = T))
# cat(list.files(path = getwd(), all.files = T, full.names = T, recursive = F))
# cat(list.files(path=dest, all.files = T, full.names = T, recursive = F))
# cat(list.files(path = path.package(package = "QuickBLAST"), all.files = T, full.names = T, recursive = F))

# cat(paste0(dest, "\n"))
# cat(paste0(staged_dest, "\n"))
# cat(paste0(R_PACKAGE_DIR, "\n"))
# cat(paste0(R_PACKAGE_SOURCE, "\n"))
# #R_PACKAGE_DIR
# cat(paste0( fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")),"\n"))

# cat(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs",Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT,"\n"))
# cat(paste0(Sys.glob(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs",Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT)),"\n"))

# cat(paste0( fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT,"\n"))
# cat(paste0(Sys.glob(paste0( fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT)),"\n"))

##Sys.glob(paste0("*", SHLIB_EXT)
#if(grepl(Sys.getenv("R_PLATFORM"), pattern ="linux", ignore.case = T)){
#  files <- c(Sys.glob(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs"), .Platform$file.sep,"*", SHLIB_EXT)))
#}else{
#  files <- c(Sys.glob(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs",Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT)))
#}
## if(WINDOWS) files <- c(files, list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=SHLIB_EXT, full.names = T))

## cat(files)
## cat(dest)
#file.copy(files, dest, overwrite = TRUE)
## file.copy(files, staged_dest, overwrite = TRUE)
## file.copy(files, fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), overwrite = FALSE)

# Look for the DLL. Adjust this path if your CMake outputs it elsewhere!
dll_paths <- c(
  "QuickBLAST.so",
  "QuickBLAST.dll",                     # If SHLIB_LINK put it in src/
  "../BUILD/QuickBLAST.so",
  "../BUILD/QuickBLAST.dll",            # If it's still in the build dir
  "../BUILD/libQuickBLAST.so",
  "../BUILD/libQuickBLAST.dll"
)

src_dir <- getwd()
# print(list.files(src_dir))
# print(dest)
dll_found <- FALSE
for (dll in dll_paths) {
  # print(file.path(src_dir,dll))
  # if (fs::file_exists(file.path(src_dir,dll))) {
  if (fs::file_exists(dll)) {
    file.copy(dll, dest, overwrite = TRUE)
    dll_found <- TRUE
    break
  }
}

if (!dll_found) {
  stop("install.libs.R could not find QuickBLAST.dll to copy!")
}

if(file.exists("symbols.rds"))
  file.copy("symbols.rds", dest, overwrite = TRUE)

# # --- Copy makeblastdb executable to the architecture-specific folder ---
# r_arch <- Sys.getenv("R_ARCH") # Evaluates to "/x64" or ""
# 
# # Look for the executable in all the locations CMake might have left it
# exe_paths <- c(
#   "makeblastdb",
#   "makeblastdb.exe",
#   # 1. Look in the exact architecture folders
#   paste0("../inst/libs", r_arch, "/makeblastdb.exe"),
#   paste0("../inst/libs", r_arch, "/makeblastdb"),
#   paste0("../libs", r_arch, "/makeblastdb.exe"),
#   paste0("../libs", r_arch, "/makeblastdb"),
#   # 2. Bulletproof Fallback: Grab it straight from the raw NCBI build output!
#   "../BUILD/ncbi-cxx-toolkit-public/BUILD/bin/makeblastdb.exe",
#   "../BUILD/ncbi-cxx-toolkit-public/BUILD/bin/makeblastdb"
# )
# 
# exe_found <- FALSE
# for (exe in exe_paths) {
#   if (fs::file_exists(exe)) {
#     # Copy it to the same 'dest' (libs/x64) as the DLL
#     file.copy(exe, dest, overwrite = TRUE)
#     exe_found <- TRUE
#     break
#   }
# }
# 
# if (!exe_found) {
#   stop("install.libs.R could not find makeblastdb executable to copy!")
# }

# require(Rcpp)
# require(RcppProgress)

# R_dll_paths <- c(
#   # list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=".dll", full.names = T),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Riconv.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"R.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rgraphapp.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rblas.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"R.dll"),
#   fs::path_package("QuickBLAST","libs", Sys.getenv("R_ARCH"),"Rlapack.dll")
#   # fs::path_package("Rcpp","libs", Sys.getenv("R_ARCH"),"Rcpp.dll")
# )


# for (dll_path in c(R_dll_paths)) {
#     if (!file.exists(dll_path)) {
#       cat("R DLL file not found:", dll_path, "\n")
#     } else {
#       dyn.load(dll_path, local=F, now = T)
#       # if(!invisible(is.loaded(dll_path))){
#       #   dyn.load(dll_path,now = T)
#       # }
#       cat("Loaded R DLL:", dll_path, "\n")
#     }
#   }