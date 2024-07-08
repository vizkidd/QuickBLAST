# $Id$

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)
project(QuickBLAST)
# option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
# option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

if(NOT R_SHELL)
    message(FATAL_ERROR "Install through R. Use devtools package. Does not support generic CMAKE compilation")
endif(NOT R_SHELL)

# if(MSVC)
execute_process(
    # COMMAND "ping www.google.com -n 2 > /dev/null"
    COMMAND ${RTOOLS_DIR}/usr/bin/wget --delete-after www.google.com
    ERROR_QUIET
    # COMMAND_ERROR_IS_FATAL ANY
    OUTPUT_QUIET
    RESULT_VARIABLE NO_CONNECTION
)
# else()
#     execute_process(
#         COMMAND ping www.google.com -c 2
#         ERROR_QUIET
#         RESULT_VARIABLE NO_CONNECTION
#     )
# endif()

# message(${NO_CONNECTION})
if(NO_CONNECTION GREATER 0)
    #offline
    message(FATAL_ERROR "Requires an internet connection")
endif()

if(NOT WIN32)    

    execute_process(
    COMMAND nproc;
    OUTPUT_VARIABLE N_PROC
    OUTPUT_STRIP_TRAILING_WHITESPACE
	)

    execute_process(
        COMMAND gcc -dumpversion;
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
    COMMAND git clone https://github.com/apache/arrow
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    message("If build fails, change -DARROW_DEPENDENCY_SOURCE=AUTO TO -DARROW_DEPENDENCY_SOURCE=BUNDLED or install arrow with package manager (Default paths are seached for lib)")

    execute_process(
    COMMAND cmake .. -G "Unix Makefiles" -DARROW_DEPENDENCY_USE_SHARED=OFF -DARROW_BUILD_SHARED=OFF -DARROW_BUILD_STATIC=ON -DARROW_BZ2_USE_SHARED=OFF -DARROW_COMPUTE=ON -DARROW_CSV=ON -DARROW_DEPENDENCY_SOURCE=AUTO -DARROW_ENABLE_THREADING=ON -DARROW_FILESYSTEM=ON -DARROW_HDFS=ON -DARROW_IPC=ON -DARROW_JEMALLOC=ON -DARROW_JSON=ON -DARROW_PARQUET=ON -DARROW_POSITION_INDEPENDENT_CODE=ON -DARROW_TENSORFLOW=ON -DARROW_VERBOSE_THIRDPARTY_BUILD=ON -DARROW_WITH_BACKTRACE=ON -DARROW_WITH_BROTLI=ON -DARROW_WITH_BZ2=ON -DARROW_WITH_LZ4=ON -DARROW_WITH_RE2=ON -DARROW_WITH_SNAPPY=ON -DARROW_WITH_UTF8PROC=ON -DARROW_WITH_ZLIB=OFF -DARROW_WITH_ZSTD=ON -DCMAKE_BUILD_TYPE=release -DBUILD_SHARED_LIBS=OFF 	-DCMAKE_POSITION_INDEPENDENT_CODE=ON #-DARROW_ZSTD_USE_SHARED=OFF -DARROW_JEMALLOC_USE_SHARED=OFF -DARROW_LLVM_USE_SHARED=OFF -DARROW_LZ4_USE_SHARED=OFF -DARROW_OPENSSL_USE_SHARED=OFF -DARROW_THRIFT_USE_SHARED=OFF -DARROW_UTF8PROC_USE_SHARED=OFF -DARROW_BOOST_USE_SHARED=OFF -DARROW_BROTLI_USE_SHARED=OFF -DARROW_PROTOBUF_USE_SHARED=OFF -DARROW_SNAPPY_USE_SHARED=OFF # -DARROW_CUDA=ON -DARROW_DEPENDENCY_SOURCE=BUNDLED  # -DARROW_BUILD_INTEGRATION=ON
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
    COMMAND make -j ${N_PROC}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    #COMMAND_ERROR_IS_FATAL ANY
    )

    # Clone the repository using Git
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND_ERROR_IS_FATAL ANY
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/src " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    # Copy files & folders to the repo clone
    file(COPY ${CMAKE_SOURCE_DIR}/QuickBLAST.prj 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/include " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/include
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    execute_process(
    COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cmake-configure
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    #set(TOOLKIT_DIR ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)
    #set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build)

    execute_process(
    COMMAND ./cmake-configure GCC ${GCC_VERSION} --with-projects=QuickBLAST.prj --without-debug --with-build-root=BUILD # --with-dll
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
    COMMAND make -j ${N_PROC}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build
    COMMAND_ERROR_IS_FATAL ANY
    )
elseif(WIN32)
  
    # if(NOT ENV{NCBITK_ROOT_DIR}) - setting this from Makevars for ncbi_components.cmake
    #     message(FATAL_ERROR "Export NCBI toolkit ROOT path to NCBITK_ROOT_DIR env variable")
    # endif(NOT ENV{NCBITK_ROOT_DIR})

    message(STATUS "Rtools Install Path:" ${RTOOLS_DIR})

    execute_process(
        COMMAND nproc;
        OUTPUT_VARIABLE N_PROC
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/gcc -dumpversion;
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
        COMMAND pacman -Syuu -q --needed --noconfirm
    )

    execute_process(
        COMMAND pacman -Sy -q --needed --noconfirm development mingw-w64-cross-toolchain libraries compression mingw-w64-cross sys-utils mingw-w64-clang-x86_64-toolchain mingw-w64-ucrt-x86_64-toolchain mingw-w64-x86_64-toolchain
    )

    execute_process(
        COMMAND pacman -Sy -q --needed --noconfirm cmake mingw-w64-clang-x86_64-cmake mingw-w64-ucrt-x86_64-cmake git mingw-w64-ucrt-x86_64-libutf8proc mingw-w64-ucrt-x86_64-utf8cpp mingw-w64-clang-x86_64-utf8cpp mingw-w64-clang-x86_64-libutf8proc mingw-w64-x86_64-libutf8proc mingw-w64-x86_64-utf8cpp mingw-w64-ucrt-x86_64-xsimd mingw-w64-x86_64-xsimd mingw-w64-clang-x86_64-xsimd mingw-w64-ucrt-x86_64-re2 mingw-w64-clang-x86_64-re2 mingw-w64-x86_64-re2 mingw-w64-ucrt-x86_64-boost mingw-w64-clang-x86_64-boost mingw-w64-x86_64-boost mingw-w64-ucrt-x86_64-thrift mingw-w64-clang-x86_64-thrift inetutils mingw-w64-ucrt-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-python-win_inet_pton mingw-w64-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-clang mingw-w64-ucrt-x86_64-clang clang llvm mingw-w64-clang-x86_64-llvm mingw-w64-ucrt-x86_64-llvm mingw-w64-x86_64-llvm mingw-w64-ucrt-x86_64-llvm-libs mingw-w64-clang-x86_64-llvm-libs mingw-w64-x86_64-llvm-libs lzop mingw-w64-clang-x86_64-lzo2 mingw-w64-ucrt-x86_64-lzo2 mingw-w64-x86_64-lzo2 liblzo2 liblzo2-devel mingw-w64-x86_64-openmp mingw-w64-ucrt-x86_64-openmp mingw-w64-clang-x86_64-openmp brotli mingw-w64-clang-x86_64-brotli mingw-w64-ucrt-x86_64-brotli mingw-w64-x86_64-brotli brotli-devel mingw-w64-ucrt-x86_64-samtools mingw-w64-clang-x86_64-samtools mingw-w64-x86_64-samtools mingw-w64-ucrt-x86_64-arrow mingw-w64-x86_64-arrow liblz4 lz4 mingw-w64-clang-x86_64-lz4 mingw-w64-ucrt-x86_64-lz4 mingw-w64-x86_64-lz4 zlib zlib-devel mingw-w64-cross-zlib mingw-w64-x86_64-zlib mingw-w64-ucrt-x86_64-zlib mingw-w64-clang-x86_64-zlib mingw-w64-clang-x86_64-python-conan mingw-w64-ucrt-x86_64-python-conan mingw-w64-x86_64-libc++ mingw-w64-ucrt-x86_64-libc++ mingw-w64-clang-x86_64-libc++ mingw-w64-clang-x86_64-headers-git mingw-w64-ucrt-x86_64-headers-git mingw-w64-x86_64-headers-git mingw-w64-ucrt-x86_64-libunwind mingw-w64-x86_64-libunwind mingw-w64-clang-x86_64-libunwind mingw-w64-x86_64-libbacktrace mingw-w64-ucrt-x86_64-libbacktrace mingw-w64-clang-x86_64-libbacktrace mingw-w64-x86_64-python-conan mingw-w64-clang-x86_64-python-conan mingw-w64-ucrt-x86_64-python-conan make mingw-w64-x86_64-clang mingw-w64-x86_64-libxml2 mingw-w64-clang-x86_64-libxml2 mingw-w64-ucrt-x86_64-libxml2 libxml2 mingw-w64-clang-x86_64-lld mingw-w64-ucrt-x86_64-lld mingw-w64-x86_64-lld lld autotools mingw-w64-clang-x86_64-autotools mingw-w64-ucrt-x86_64-autotools mingw-w64-x86_64-autotools libgcrypt-devel libgcrypt mingw-w64-x86_64-libgcrypt mingw-w64-ucrt-x86_64-libgcrypt mingw-w64-clang-x86_64-libgcrypt bison bisonc++ flex reflex flexc++ mingw-w64-ucrt-x86_64-flexdll mingw-w64-x86_64-flexdll msys2-runtime-devel msys2-runtime  msys2-w32api-runtime mingw-w64-cross-binutils mingw-w64-ucrt-x86_64-binutils mingw-w64-x86_64-binutils binutils mingw-w64-clang-x86_64-lmdb mingw-w64-ucrt-x86_64-lmdb mingw-w64-x86_64-lmdb mingw-w64-clang-x86_64-lmdbxx mingw-w64-ucrt-x86_64-lmdbxx mingw-w64-x86_64-lmdbxx ninja mingw-w64-x86_64-arrow mingw-w64-ucrt-x86_64-arrow mingw-w64-x86_64-cmake pcre pcre2 mingw-w64-x86_64-pcre mingw-w64-x86_64-pcre2 mingw-w64-ucrt-x86_64-pcre mingw-w64-ucrt-x86_64-pcre2 mingw-w64-clang-x86_64-pcre mingw-w64-clang-x86_64-pcre2 libpcreposix libpcre2posix pcre2-devel pcre-devel libpcrecpp libpcre libpcre2_32 grep mingw-w64-x86_64-grep mingw-w64-clang-x86_64-grep mingw-w64-ucrt-x86_64-grep
    )

    # execute_process(
    #     COMMAND git clone https://github.com/apache/thrift
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}    
    # )

    # execute_process(
    #     COMMAND bash -c " git tag -l | sort -r -n -t '.' -k2,1 | grep -vi \"thrift\" | head -n 1" ;
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/thrift
    #     OUTPUT_VARIABLE THRIFT_GIT_TAG
    # 	OUTPUT_STRIP_TRAILING_WHITESPACE
    # )

    # message("THRIFT GIT LATEST TAG : " ${THRIFT_GIT_TAG})

    # execute_process(
    #     COMMAND git checkout ${THRIFT_GIT_TAG}
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/thrift    
    #     COMMAND_ERROR_IS_FATAL ANY
    # )

    # execute_process(
    #     COMMAND git clone https://github.com/llvm/llvm-project
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}    
    # )

    # execute_process(
    #     COMMAND bash -c "git tag -l | sort -r -u -h -k3 -t '-' | grep -iv \"dev\" | head -n 1" ;
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/llvm-project
    #     OUTPUT_VARIABLE LLVM_GIT_TAG
    # 	OUTPUT_STRIP_TRAILING_WHITESPACE
    # )

    # message("LLVM GIT LATEST TAG : " ${LLVM_GIT_TAG})

    # execute_process(
    #     COMMAND git clone https://github.com/Tencent/rapidjson
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}    
    # )

    #  git tag -l | sort -r | head -n 1

    #  git submodule update --init

    if (NOT EXISTS "${ROOT_SRC_DIR}/R_src/")
	message(STATUS "Downloading R-src...")	
	execute_process(
			COMMAND bash -c "mkdir -p ${ROOT_SRC_DIR}/R_src/" ;
			WORKING_DIRECTORY ${ROOT_SRC_DIR}
		)
		execute_process(
			COMMAND bash -c "curl https://cran.r-project.org/src/base/R-latest.tar.gz -o R-latest.tar.gz" ;
			WORKING_DIRECTORY ${ROOT_SRC_DIR}/R_src/
			COMMAND_ERROR_IS_FATAL ANY
		)
		execute_process(
			COMMAND bash -c "tar xzf R-latest.tar.gz 2> /dev/null" ;
			WORKING_DIRECTORY ${ROOT_SRC_DIR}/R_src/
		)
	endif()

    if (NOT EXISTS "${ROOT_SRC_DIR}/arrow/")
    execute_process(
        COMMAND git clone https://github.com/apache/arrow
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND_ERROR_IS_FATAL ANY    
    )
    endif()

    execute_process(
        COMMAND bash -c "git tag -l | sort -r -u -h -k3 -t '-' | grep -iv \"dev\" | head -n 1" ; #no dev versions
        # COMMAND bash -c "git tag -l | sort -r -u -h -k3 -t '-' | head -n 1" ; #- include dev versions
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow
        OUTPUT_VARIABLE ARROW_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("ARROW GIT LATEST TAG : " ${ARROW_GIT_TAG})

    execute_process(
        COMMAND git checkout ${ARROW_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow    
        COMMAND_ERROR_IS_FATAL ANY
    )

    # file(COPY ${CMAKE_SOURCE_DIR}/arrow_defs.h
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/arrow/
    # )

    # execute_process(
    #     COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # )

    # message("Copying " ${RTOOLS_DIR}/usr/include/netinet " into " ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/)
    # file(COPY ${RTOOLS_DIR}/usr/include/netinet 
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/
    # )

    # message("Copying " ${RTOOLS_DIR}/usr/include/cygwin " into " ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/)
    # file(COPY ${RTOOLS_DIR}/usr/include/cygwin 
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/
    # )

    # message("Copying " ${RTOOLS_DIR}/usr/include/sys " into " ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/)
    # file(COPY ${RTOOLS_DIR}/usr/include/sys 
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/src/
    # )

    # execute_process(
    #     COMMAND cmake .. -G "MSYS Makefiles" -DARROW_CSV=ON -DARROW_FILESYSTEM=ON -DARROW_HDFS=ON -DARROW_PARQUET=ON -DARROW_DEPENDENCY_USE_SHARED=OFF -DARROW_WITH_BZ2=ON  -DARROW_WITH_LZ4=ON -DARROW_DEPENDENCY_SOURCE=BUNDLED -DARROW_BUILD_STATIC=ON -DARROW_BUILD_SHARED=OFF -DPARQUET_MINIMAL_DEPENDENCY=OFF -DARROW_WITH_BROTLI=ON -DARROW_WITH_ZLIB=OFF -DARROW_WITH_ZSTD=ON -DCMAKE_BUILD_TYPE=Release -DARROW_IPC=ON -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON -DARROW_COMPUTE=ON #-DARROW_THRIFT_USE_SHARED=OFF -DARROW_BZ2_USE_SHARED=OFF -DARROW_LZ4_USE_SHARED=OFF -DARROW_BOOST_USE_SHARED=OFF -DARROW_BROTLI_USE_SHARED=OFF -DARROW_PROTOBUF_USE_SHARED=OFF -DARROW_SNAPPY_USE_SHARED=OFF -DARROW_UTF8PROC_USE_SHARED=OFF -DARROW_ZSTD_USE_SHARED=OFF
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    #     COMMAND_ERROR_IS_FATAL ANY    
    # )

    # execute_process(
    #     COMMAND make -j ${N_PROC}
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    #     COMMAND_ERROR_IS_FATAL ANY
    # )

    # execute_process(
    # COMMAND git clone https://github.com/libunwind/libunwind
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # #COMMAND_ERROR_IS_FATAL ANY
    # )

#    execute_process(
#    COMMAND autoreconf -i && ${CMAKE_CURRENT_BINARY_DIR}/libunwind/configure --prefix=BUILD && make && make install
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libunwind
#    COMMAND_ERROR_IS_FATAL ANY
#    )

    # execute_process(
    # COMMAND git clone https://github.com/bombela/backward-cpp
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # #COMMAND_ERROR_IS_FATAL ANY
    # )

    # set(Backward_DIR ${CMAKE_CURRENT_BINARY_DIR}/backward-cpp)
    # # list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR}/backward-cpp)
    # find_package(Backward)

    if (NOT EXISTS "${ROOT_SRC_DIR}/ncbi-cxx-toolkit-public/")
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND_ERROR_IS_FATAL ANY
    )
    endif()

    execute_process(
        COMMAND bash -c "git tag -l | sort -r | head -n 1" ;
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public
        OUTPUT_VARIABLE NCBI_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("NCBI GIT LATEST TAG : " ${NCBI_GIT_TAG})

    execute_process(
        COMMAND git checkout ${NCBI_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public 
        COMMAND_ERROR_IS_FATAL ANY
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/src " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    # Copy files & folders to the repo clone
    # file(COPY ${CMAKE_SOURCE_DIR}/setup-x86_64.exe 
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cygwin
    # )

    file(COPY ${CMAKE_SOURCE_DIR}/dependencies/QuickBLAST.prj 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}//dependencies/ncbi_crosscomp_inject_pre.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/dependencies/ncbi_crosscomp_inject_post.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/dependencies/ncbi_components.cmake 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    # file(COPY ${CMAKE_SOURCE_DIR}/ncbi_preprocessor_hooks.h 
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    # )

    ## file(COPY ${CMAKE_SOURCE_DIR}/arrow_int64_t_linker.ld
    ## DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    ## )

    file(COPY ${CMAKE_SOURCE_DIR}/src 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/include " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/include
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    # file(COPY ${RTOOLS_DIR}/mingw64/include/arrow
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/include
    # )

    execute_process(
        COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cmake-configure ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/configure
    )

    # execute_process(
    # COMMAND git clone https://github.com/Mbed-TLS/mbedtls
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # #COMMAND_ERROR_IS_FATAL ANY
    # )

    # #Copying mbedtls
    # message("Copying " ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/library/* " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/)
    # message("Copying " ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/include/mbedtls " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/mbedtls/)
    # message("Copying " ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/include/psa " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/mbedtls/)

    # file(GLOB MBEDTLS_HEADERS LIST_DIRECTORIES false
    # ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/library/*.h
    # ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/library/*.c
    # )
    # # message(${MBEDTLS_HEADERS})
    # file(COPY ${MBEDTLS_HEADERS}
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/mbedtls/
    # )
    # file(COPY ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/include/mbedtls
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/mbedtls/
    # )
    # file(COPY ${CMAKE_CURRENT_BINARY_DIR}/mbedtls/include/psa
    # DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/connect/mbedtls/
    # )

    set(CMAKE_CXX_EXTENSIONS ON )
    set(CMAKE_SYSROOT_LINK "/x86_64-w64-mingw32.static.posix/lib")
    # execute_process(
    # COMMAND bash -c "ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/cmkTool.sh GCC ${GCC_VERSION} | grep 'in' | awk -F' ' '{print $NF}' "
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # OUTPUT_VARIABLE TOOLCHAIN_FILE
    # OUTPUT_STRIP_TRAILING_WHITESPACE
    # COMMAND_ERROR_IS_FATAL ANY
    # )

    # if(TOOLCHAIN_FILE STREQUAL "")
    #     set(TOOLCHAIN_FILE toolchain.in)
    #     message("TOOLCHAIN CMAKE : " ${TOOLCHAIN_FILE})
    # else()
    #     execute_process(
    #     COMMAND bash -c "basename ${TOOLCHAIN_FILE} "
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #     OUTPUT_VARIABLE TOOLCHAIN_FILE
    #     OUTPUT_STRIP_TRAILING_WHITESPACE
    #     COMMAND_ERROR_IS_FATAL ANY
    #     )
    # endif()

    # set(ENV{CMAKE_ARGS} "-DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON") #-G 'Unix Makefiles' -G 'MSYS Makefiles' #-DUNIX=1 -DWIN32=0 -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_SYSTEM_PROCESSOR=x86_64 
    
    # # export CMAKE_CMD="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/cmake.exe"
    # # set(CMAKE_C_COMPILER "${RTOOLS_DIR}/x86_64-w64-mingw32.static.posix/bin/gcc.exe")
    # # set(CMAKE_CXX_COMPILER "${RTOOLS_DIR}/x86_64-w64-mingw32.static.posix/bin/g++.exe")

    # message("CMAKE ARGS:" ${CMAKE_ARGS})
    # message("GCC VERSION:" ${GCC_VERSION})

    # #-mingw64
    # execute_process(
    # COMMAND ${RTOOLS_DIR}/msys2_shell.cmd -where ${CMAKE_CURRENT_BINARY_DIR} -no-start -msys -defterm -shell bash ${CMAKE_SOURCE_DIR}/ncbi_crosscomp_inject.sh ${TOOLCHAIN_FILE}
    # #COMMAND ${RTOOLS_DIR}/msys2_shell.cmd -where ${CMAKE_CURRENT_BINARY_DIR} -no-start -mingw64 -defterm -shell bash ${CMAKE_SOURCE_DIR}/ncbi_crosscomp_inject.sh ${TOOLCHAIN_FILE}
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    # COMMAND_ERROR_IS_FATAL ANY
    # )

    # ##set(CMAKE_ARGS "${CMAKE_ARGS} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/toolchain.in")
    # #set(ENV{CMAKE_ARGS} "${CMAKE_ARGS} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/${TOOLCHAIN_FILE}")
    # # set(ENV{MINGW32} 1)


    # execute_process(
    # COMMAND bash cmake-configure GCC ${GCC_VERSION} --with-projects=QuickBLAST.prj --without-debug --with-build-root=BUILD #--with-conan # --with-dll
    # #COMMAND cmake-configure.bat --with-projects=QuickBLAST.prj --without-debug --with-build-root=BUILD --with-conan # --with-dll
    # #COMMAND bash configure --without-debug --with-static --with-static-exe --with-bin-release --with-openmp --with-64 --with-lfs --with-build-root=BUILD --with-projects=QuickBLAST.prj #--with-conan --with-autodep
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    # COMMAND_ERROR_IS_FATAL ANY
    # )

    # execute_process(
    # COMMAND make -j ${N_PROC}
    # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build
    # COMMAND_ERROR_IS_FATAL ANY
    # )

    # execute_process(
    #     COMMAND rm -f ${TOOLCHAIN_FILE}
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/
    # )

endif()

