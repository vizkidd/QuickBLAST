dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

# staged_dest <- file.path(R_PACKAGE_DIR,"..","..","..","QuickBLAST", paste0('libs', R_ARCH))
# dir.create(staged_dest, recursive = TRUE, showWarnings = FALSE)

#C:/Users/vishv/AppData/Local/Temp/RtmpOSDZyI/temp_libpath1bacef43ec2/00LOCK-QuickBLAST_cmake_test/00new/QuickBLAST
#  C:/Users/vishv/AppData/Local/Temp/RtmpO6DZ4p/temp_libpath2fb068301991/QuickBLAST/ 

cat(paste("INSTALLING....", R_PACKAGE_NAME, "\n", sep=""))
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

#Sys.glob(paste0("*", SHLIB_EXT)
files <- c(Sys.glob(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs",Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT)))

# if(WINDOWS) files <- c(files, list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=SHLIB_EXT, full.names = T))

# cat(files)

file.copy(files, dest, overwrite = TRUE)
# file.copy(files, staged_dest, overwrite = TRUE)
# file.copy(files, fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), overwrite = FALSE)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)

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