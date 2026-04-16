#Set NCBI Components which give issues during build
# DO NOT print messages here as this is included with the toolchain and will plague the output, only print for DEBUGGING

if(NOT NCBITK_ROOT_DIR)
    set(NCBITK_ROOT_DIR $ENV{NCBITK_ROOT_DIR})
endif(NOT NCBITK_ROOT_DIR)

#lmdb
if(NOT NCBI_COMPONENT_LMDB_FOUND AND NOT HAVE_LIBLMDB)
    set(NCBI_COMPONENT_LocalLMDB_FOUND LOCAL)
    set(NCBI_COMPONENT_LocalLMDB_INCLUDE ${NCBITK_ROOT_DIR}/include/util/lmdb)
    set(NCBI_COMPONENT_LocalLMDB_NCBILIB lmdb)
    set(NCBI_COMPONENT_LMDB_FOUND ${NCBI_COMPONENT_LocalLMDB_FOUND})
    set(NCBI_COMPONENT_LMDB_INCLUDE ${NCBI_COMPONENT_LocalLMDB_INCLUDE})
    set(NCBI_COMPONENT_LMDB_NCBILIB ${NCBI_COMPONENT_LocalLMDB_NCBILIB})
    set(HAVE_LIBLMDB ${NCBI_COMPONENT_LMDB_FOUND})
endif()

#sqlite
# 1. Step up one directory from ncbi-cxx-toolkit-public to the BUILD directory
get_filename_component(QUICKBLAST_BUILD_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

# 2. Define the paths to your pre-built SQLite library
set(CUSTOM_SQLITE_DIR "${QUICKBLAST_BUILD_DIR}/sqlite")
set(CUSTOM_SQLITE_LIB "${CUSTOM_SQLITE_DIR}/libsqlite3.a")

# Find the amalgamation folder where sqlite3.h lives
file(GLOB CUSTOM_SQLITE_INC "${CUSTOM_SQLITE_DIR}/sqlite-amalgamation-*")

message(STATUS "CUSTOM_SQLITE_LIB: ${CUSTOM_SQLITE_LIB}")

if(EXISTS "${CUSTOM_SQLITE_LIB}")
    message(STATUS "Injecting pre-built SQLite paths into NCBI scanner...")
    
    # INSTEAD of forcing NCBI's cache variables, append your paths to CMake's native search lists.
    # When NCBI runs its normal system scan, it will look here first and organically find it!
    list(APPEND CMAKE_PREFIX_PATH "${CUSTOM_SQLITE_DIR}")
    list(APPEND CMAKE_INCLUDE_PATH "${CUSTOM_SQLITE_INC}")
    list(APPEND CMAKE_LIBRARY_PATH "${CUSTOM_SQLITE_DIR}")
    
    # Also set the Environment variables as a bulletproof backup
    set(ENV{CMAKE_INCLUDE_PATH} "${CUSTOM_SQLITE_INC}:$ENV{CMAKE_INCLUDE_PATH}")
    set(ENV{CMAKE_LIBRARY_PATH} "${CUSTOM_SQLITE_DIR}:$ENV{CMAKE_LIBRARY_PATH}")
    
else()
    message(FATAL_ERROR "Custom SQLite lib NOT FOUND at: ${CUSTOM_SQLITE_LIB}")
endif()