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


files <- c(Sys.glob(paste0(file.path(R_PACKAGE_SOURCE,"inst","libs",Sys.getenv("R_ARCH")), .Platform$file.sep,"*", SHLIB_EXT)),Sys.glob(paste0("*", SHLIB_EXT)))

# if(WINDOWS) files <- c(files, list.files(file.path(Sys.getenv("R_HOME"),"bin",Sys.getenv("R_ARCH")),pattern=SHLIB_EXT, full.names = T))

# cat(files)

file.copy(files, dest, overwrite = TRUE)
# file.copy(files, staged_dest, overwrite = TRUE)
# file.copy(files, fs::path_package("QuickBLAST", "libs", Sys.getenv("R_ARCH")), overwrite = FALSE)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)

# require(Rcpp)
# require(RcppProgress)