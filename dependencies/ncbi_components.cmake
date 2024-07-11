#Set NCBI Components which give issues during build
# DO NOT print messages here as this is included with the toolchain and will plague the output, only print for DEBUGGING

# set(FOLDER_MODE "-m") #FOR cygpath
# function(cygpath_cmd ret_var _folder_path)
#     block(SCOPE_FOR VARIABLES)
#         # message(STATUS ${_folder_path})
#         execute_process(
#             COMMAND bash -c "cygpath ${FOLDER_MODE} ${_folder_path}" ;
#             OUTPUT_VARIABLE _return_path
#             OUTPUT_STRIP_TRAILING_WHITESPACE
#         )
#         # message(STATUS ${ret_var} ${_return_path})
#         set(${ret_var} ${_return_path})
#         return(PROPAGATE ${ret_var} ${_return_path})
#     endblock()
# endfunction(cygpath_cmd)



if(NOT NCBITK_ROOT_DIR)
    # set(NCBITK_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)
    # set(NCBITK_ROOT_DIR ${NCBITK_INC_ROOT}/include/)
    # set(NCBITK_ROOT_DIR ${NCBITK_INC_ROOT}/../)
    set(NCBITK_ROOT_DIR $ENV{NCBITK_ROOT_DIR})
    # message("NCBI ROOT : " ${NCBITK_ROOT_DIR}) #DBG
endif(NOT NCBITK_ROOT_DIR)

#DLL
#set(BUILD_SHARED_LIBS ON)
# set(NCBI_REQUIRE_DLL_BUILD_FOUND YES)
# set(NCBI_REQUIRE_DLL_FOUND YES)
# set(NCBI_PTBCFG_ALLOW_COMPOSITE ON)
# set(NCBI_DLL_BUILD 1)
# set(NCBI_DLL_SUPPORT 1)


# cygpath_cmd(NCBITK_ROOT_DIR "${NCBITK_ROOT_DIR}")

#lmdb
if(NOT NCBI_COMPONENT_LMDB_FOUND AND NOT HAVE_LIBLMDB)
    set(NCBI_COMPONENT_LocalLMDB_FOUND LOCAL)
    set(NCBI_COMPONENT_LocalLMDB_INCLUDE ${NCBITK_ROOT_DIR}/include/util/lmdb)
    set(NCBI_COMPONENT_LocalLMDB_NCBILIB lmdb)
    set(NCBI_COMPONENT_LMDB_FOUND ${NCBI_COMPONENT_LocalLMDB_FOUND})
    set(NCBI_COMPONENT_LMDB_INCLUDE ${NCBI_COMPONENT_LocalLMDB_INCLUDE})
    set(NCBI_COMPONENT_LMDB_NCBILIB ${NCBI_COMPONENT_LocalLMDB_NCBILIB})
    set(HAVE_LIBLMDB ${NCBI_COMPONENT_LMDB_FOUND})
    # message("LMDB INC : " ${NCBI_COMPONENT_LMDB_INCLUDE}) #DBG
endif()