# $Id$

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)
project(QuickBLAST)

if(NOT R_SHELL)
    message(FATAL_ERROR "Install through R. Use devtools package. Does not support generic CMAKE compilation")
endif(NOT R_SHELL)

execute_process(
    COMMAND ${RTOOLS_DIR}/usr/bin/wget --delete-after www.google.com
    ERROR_QUIET
    OUTPUT_QUIET
    RESULT_VARIABLE NO_CONNECTION
)

if(NO_CONNECTION GREATER 0)
    #offline
    message(FATAL_ERROR "Requires an internet connection")
endif()

#message("DTOOLKIT_CPPFLAGS ${TOOLKIT_CPPFLAGS}")
#message("DRCPP_FOLDER ${RCPP_FOLDER}")
#message("DRCPPPROG_FOLDER ${RCPPPROG_FOLDER}")
#return()

# Standard configuration 
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Use standard C++20 without GNU extensions for better portability

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

        if (NOT EXISTS "${ROOT_SRC_DIR}/arrow/")
    execute_process(
        COMMAND git clone https://github.com/apache/arrow arrow
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
     #   COMMAND_ERROR_IS_FATAL ANY    
    )
    endif()

    #execute_process(
    #    COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r -u -h -k3 -t '-' | grep -iv \"dev\" | head -n 1" ; #no dev versions
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow
    #    OUTPUT_VARIABLE ARROW_GIT_TAG
    #	OUTPUT_STRIP_TRAILING_WHITESPACE
    #)

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l 'apache-arrow-*.*.*' --sort=-v:refname | grep -iv \"dev\" | grep -iv \"rc\" | grep -iv \"js\" | head -n 1"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow
        OUTPUT_VARIABLE ARROW_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    message("ARROW GIT LATEST TAG : " ${ARROW_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${ARROW_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow    
        COMMAND_ERROR_IS_FATAL ANY
    )

    if (NOT EXISTS "${ROOT_SRC_DIR}/ncbi-cxx-toolkit-public/")
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND_ERROR_IS_FATAL ANY
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l --sort=-v:refname | grep -iv \"dev\" | head -n 1" ; #no dev versions
        #COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r | head -n 1" ;
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public
        OUTPUT_VARIABLE NCBI_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("NCBI GIT LATEST TAG : " ${NCBI_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${NCBI_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public 
        COMMAND_ERROR_IS_FATAL ANY
    )

    #execute_process(
    #    COMMAND git clone https://github.com/NixOS/patchelf patchelf
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #    #COMMAND_ERROR_IS_FATAL ANY    
    #)
    #execute_process(
    #    COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l --sort=-v:refname | grep -iv \"dev\" | head -n 1" ; #no dev versions
    #    #COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r -u -h -k3 -t '-' | grep -iv \"dev\" | head -n 1" ; #no dev versions
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/patchelf
    #    OUTPUT_VARIABLE PATCHELF_GIT_TAG
    #	OUTPUT_STRIP_TRAILING_WHITESPACE
    #)
    #message("PATCHELF GIT LATEST TAG : " ${PATCHELF_GIT_TAG})
    #execute_process(
    #    COMMAND git checkout -q ${PATCHELF_GIT_TAG}
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/patchelf
    #    COMMAND_ERROR_IS_FATAL ANY
    #)
    
    message("Copying " ${CMAKE_SOURCE_DIR}/src " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/QuickBLAST.prj 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_crosscomp_inject_pre.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_crosscomp_inject_post.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_components.cmake 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/ 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/src/include " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/src/include
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    execute_process(
        COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cmake-configure ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/configure
    )

    set(CMAKE_CXX_EXTENSIONS OFF)

    execute_process(
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    message("If build fails, change -DARROW_DEPENDENCY_SOURCE=AUTO TO -DARROW_DEPENDENCY_SOURCE=BUNDLED or install arrow with package manager (Default paths are seached for lib)")
  
  	execute_process(
      COMMAND pkg-config --cflags libR;
      OUTPUT_VARIABLE R_FLAGS
      OUTPUT_STRIP_TRAILING_WHITESPACE
  	)
  	message("R flags:" ${R_FLAGS})

#-fdiagnostics-color=always -Wall
#-std=gnu++20
#-std=gnu2x
#-O3 -fPIC  -fdiagnostics-color=always

if(APPLE)
    # Apple Clang (macOS) compatible flags
    set(ARROW_CXX_FLAGS " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -pthread")
    set(ARROW_C_FLAGS  " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -pthread")
else()
    # GCC compatible flags
    set(ARROW_CXX_FLAGS " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -mtune=generic -march=native -maccumulate-outgoing-args -pthread")
    set(ARROW_C_FLAGS  " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc -mtune=generic -maccumulate-outgoing-args -pthread ")
endif()

#set(ARROW_CXX_FLAGS " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -mtune=generic -march=native -maccumulate-outgoing-args -pthread")
#set(ARROW_C_FLAGS  " -fpermissive ${R_FLAGS} -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc -mtune=generic -maccumulate-outgoing-args -pthread ")

set(ARROW_EXE_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")
set(ARROW_LINKER_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")


    execute_process(
    COMMAND cmake .. -G "Unix Makefiles" 
    "-DCMAKE_CXX_FLAGS=${ARROW_CXX_FLAGS}"
    "-DCMAKE_C_FLAGS=${ARROW_C_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${ARROW_EXE_FLAGS}"
    "-DCMAKE_SHARED_LINKER_FLAGS=${ARROW_LINKER_FLAGS}"
    -DARROW_DEPENDENCY_USE_SHARED=OFF -DARROW_BUILD_SHARED=OFF -DARROW_BUILD_STATIC=ON -DPARQUET_BUILD_SHARED=OFF -DPARQUET_BUILD_STATIC=ON -DARROW_BZ2_USE_SHARED=OFF -DARROW_COMPUTE=ON -DARROW_CSV=ON -DARROW_DEPENDENCY_SOURCE=AUTO -DARROW_ENABLE_THREADING=ON -DARROW_FILESYSTEM=ON -DARROW_HDFS=ON -DARROW_IPC=ON -DARROW_JEMALLOC=ON -DARROW_JSON=ON -DARROW_PARQUET=ON -DARROW_POSITION_INDEPENDENT_CODE=ON -DARROW_TENSORFLOW=ON -DARROW_VERBOSE_THIRDPARTY_BUILD=ON -DARROW_WITH_BACKTRACE=ON -DARROW_WITH_BROTLI=ON -DARROW_WITH_BZ2=ON -DARROW_WITH_LZ4=ON -DARROW_WITH_RE2=ON -DARROW_RPATH_ORIGIN=ON -DARROW_WITH_SNAPPY=ON -DARROW_WITH_UTF8PROC=ON -DARROW_WITH_ZLIB=OFF -DARROW_WITH_ZSTD=ON -DCMAKE_BUILD_TYPE=release -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_SKIP_BUILD_RPATH=FALSE -DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE -DCMAKE_INSTALL_RPATH="$ORIGIN"	-DCMAKE_POSITION_INDEPENDENT_CODE=ON #-DARROW_ZSTD_USE_SHARED=OFF -DARROW_JEMALLOC_USE_SHARED=OFF -DARROW_LLVM_USE_SHARED=OFF -DARROW_LZ4_USE_SHARED=OFF -DARROW_OPENSSL_USE_SHARED=OFF -DARROW_THRIFT_USE_SHARED=OFF -DARROW_UTF8PROC_USE_SHARED=OFF -DARROW_BOOST_USE_SHARED=OFF -DARROW_BROTLI_USE_SHARED=OFF -DARROW_PROTOBUF_USE_SHARED=OFF -DARROW_SNAPPY_USE_SHARED=OFF # -DARROW_CUDA=ON -DARROW_DEPENDENCY_SOURCE=BUNDLED  # -DARROW_BUILD_INTEGRATION=ON
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
    COMMAND make -j ${N_PROC} 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    #COMMAND_ERROR_IS_FATAL ANY
    )

    ## Clone the repository using Git
    #execute_process(
    #COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public ncbi-cxx-toolkit-public
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ##COMMAND_ERROR_IS_FATAL ANY
    #)

    #message("Copying " ${CMAKE_SOURCE_DIR}/src/ " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    ## Copy files & folders to the repo clone
    #file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/QuickBLAST.prj 
    #DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    #)

    #file(COPY ${CMAKE_SOURCE_DIR}/src/ 
    #DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    #)

    #message("Copying " ${CMAKE_SOURCE_DIR}/src/include " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    #file(COPY ${CMAKE_SOURCE_DIR}/src/include
    #DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    #)
    
    #execute_process(
    #COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cmake-configure
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #)

#-std=gnu2x
#-std=gnu++20
#-O3 -fPIC  -fdiagnostics-color=always

if(APPLE)
set(ENV{CMAKE_ARGS} "-DRCPP_FOLDER=${RCPP_FOLDER} \
    -DCMAKE_C_FLAGS=' -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native' \
    -DCMAKE_CXX_FLAGS=' -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native' \
    -DRCPPPROG_FOLDER=${RCPPPROG_FOLDER} \
    -DRCPPTHREAD_FOLDER=${RCPPTHREAD_FOLDER} \
    -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
    -DCMAKE_BUILD_RPATH='-Wl,-rpath,@loader_path' \
    -DCMAKE_INSTALL_RPATH='-Wl,-rpath,@loader_path' \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF")
#    set(ENV{CMAKE_ARGS} "-DRCPP_FOLDER=${RCPP_FOLDER} \
#    -DCMAKE_C_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc' \
#    -DCMAKE_CXX_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -march=native' \
#    -DRCPPPROG_FOLDER=${RCPPPROG_FOLDER} \
#    -DRCPPTHREAD_FOLDER=${RCPPTHREAD_FOLDER} \
#    -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
#    -DCMAKE_BUILD_RPATH='-Wl,-rpath,@loader_path' \
#    -DCMAKE_INSTALL_RPATH='-Wl,-rpath,@loader_path' \
#    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
#    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF")
else()
    set(ENV{CMAKE_ARGS} "-DRCPP_FOLDER=${RCPP_FOLDER} \
    -DCMAKE_C_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc' \
    -DCMAKE_CXX_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -march=native' \
    -DRCPPPROG_FOLDER=${RCPPPROG_FOLDER} \
    -DRCPPTHREAD_FOLDER=${RCPPTHREAD_FOLDER} \
    -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_BUILD_RPATH='\$$ORIGIN' \
    -DCMAKE_INSTALL_RPATH='\$$ORIGIN' \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF")
endif()
    
    #-std=gnu2x
    #-std=gnu++20 
    
    #message("CMAKE ARGS: $ENV{CMAKE_ARGS}")
    
    #execute_process(
    #  COMMAND bash src/build-system/cmake/toolchains/cmkTool.sh GCC ${GCC_VERSION}
    #  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    #  OUTPUT_VARIABLE TOOLCHAIN_FILE
    #  OUTPUT_STRIP_TRAILING_WHITESPACE
    #  COMMAND_ERROR_IS_FATAL ANY
    #)

    #message(STATUS "Discovered Toolchain File: ${TOOLCHAIN_FILE}")
  
    execute_process(
    COMMAND ./cmake-configure GCC ${GCC_VERSION} --with-projects=QuickBLAST.prj --without-debug --with-build-root=BUILD --without-dll
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    COMMAND_ERROR_IS_FATAL ANY
    )

#    set(TOOLCHAIN_FILE "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build/ncbi_toolchain.cmake")
#    # Append the raw text to the toolchain file
#    file(APPEND "${TOOLCHAIN_FILE}" [=[
# --- Injected RPATH settings for QuickBLAST ---
# Use CACHE ... FORCE to ensure these settings override the NCBI internal defaults
#set(CMAKE_BUILD_WITH_INSTALL_RPATH ON CACHE BOOL "" FORCE)
#set(CMAKE_INSTALL_RPATH "$ORIGIN" CACHE STRING "" FORCE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF CACHE BOOL "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -Wl,-rpath,\$$ORIGIN -Wl,-rpath,@loader_path/")
#set(CMAKE_SHARED_LINKER_FLAGS_INIT "${CMAKE_SHARED_LINKER_FLAGS_INIT} -Wl,-rpath,\$$ORIGIN -Wl,-rpath,@loader_path/")
#    ]=])

    execute_process(
    COMMAND make -j ${N_PROC} datatool QuickBLASTcpp makeblastdb
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build
    COMMAND_ERROR_IS_FATAL ANY
    )
    
elseif(WIN32)
    message(STATUS "Rtools Install Path:" ${RTOOLS_DIR})

    #execute_process(
    #    COMMAND pacman -Syuu -q --needed --noconfirm
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm development mingw-w64-cross-toolchain libraries compression mingw-w64-cross sys-utils mingw-w64-clang-x86_64-toolchain mingw-w64-ucrt-x86_64-toolchain mingw-w64-x86_64-toolchain
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm cmake gcc mingw-w64-cross-mingw64-gcc mingw-w64-cross-ucrt64-gcc mingw-w64-ucrt-x86_64-gcc mingw-w64-x86_64-gcc mingw-w64-clang-x86_64-cmake mingw-w64-ucrt-x86_64-cmake git mingw-w64-ucrt-x86_64-libutf8proc mingw-w64-ucrt-x86_64-utf8cpp mingw-w64-clang-x86_64-utf8cpp mingw-w64-clang-x86_64-libutf8proc mingw-w64-x86_64-libutf8proc mingw-w64-x86_64-utf8cpp mingw-w64-ucrt-x86_64-xsimd mingw-w64-x86_64-xsimd mingw-w64-clang-x86_64-xsimd mingw-w64-ucrt-x86_64-re2 mingw-w64-clang-x86_64-re2 mingw-w64-x86_64-re2 mingw-w64-ucrt-x86_64-boost mingw-w64-clang-x86_64-boost mingw-w64-x86_64-boost mingw-w64-ucrt-x86_64-thrift mingw-w64-clang-x86_64-thrift inetutils mingw-w64-ucrt-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-python-win_inet_pton mingw-w64-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-clang mingw-w64-ucrt-x86_64-clang clang llvm mingw-w64-clang-x86_64-llvm mingw-w64-ucrt-x86_64-llvm mingw-w64-x86_64-llvm mingw-w64-ucrt-x86_64-llvm-libs mingw-w64-clang-x86_64-llvm-libs mingw-w64-x86_64-llvm-libs lzop mingw-w64-clang-x86_64-lzo2 mingw-w64-ucrt-x86_64-lzo2 mingw-w64-x86_64-lzo2 liblzo2 liblzo2-devel mingw-w64-x86_64-openmp mingw-w64-ucrt-x86_64-openmp mingw-w64-clang-x86_64-openmp brotli mingw-w64-clang-x86_64-brotli mingw-w64-ucrt-x86_64-brotli mingw-w64-x86_64-brotli brotli-devel mingw-w64-ucrt-x86_64-samtools mingw-w64-clang-x86_64-samtools mingw-w64-x86_64-samtools mingw-w64-ucrt-x86_64-arrow mingw-w64-x86_64-arrow liblz4 lz4 mingw-w64-clang-x86_64-lz4 mingw-w64-ucrt-x86_64-lz4 mingw-w64-x86_64-lz4 zlib zlib-devel mingw-w64-cross-zlib mingw-w64-x86_64-zlib mingw-w64-ucrt-x86_64-zlib mingw-w64-clang-x86_64-zlib mingw-w64-clang-x86_64-python-conan mingw-w64-ucrt-x86_64-python-conan mingw-w64-x86_64-libc++ mingw-w64-ucrt-x86_64-libc++ mingw-w64-clang-x86_64-libc++ mingw-w64-clang-x86_64-headers-git mingw-w64-ucrt-x86_64-headers-git mingw-w64-x86_64-headers-git mingw-w64-ucrt-x86_64-libunwind mingw-w64-x86_64-libunwind mingw-w64-clang-x86_64-libunwind mingw-w64-x86_64-libbacktrace mingw-w64-ucrt-x86_64-libbacktrace mingw-w64-clang-x86_64-libbacktrace mingw-w64-x86_64-python-conan mingw-w64-clang-x86_64-python-conan mingw-w64-ucrt-x86_64-python-conan make mingw-w64-x86_64-clang mingw-w64-x86_64-libxml2 mingw-w64-clang-x86_64-libxml2 mingw-w64-ucrt-x86_64-libxml2 libxml2 mingw-w64-clang-x86_64-lld mingw-w64-ucrt-x86_64-lld mingw-w64-x86_64-lld lld autotools mingw-w64-clang-x86_64-autotools mingw-w64-ucrt-x86_64-autotools mingw-w64-x86_64-autotools libgcrypt-devel libgcrypt mingw-w64-x86_64-libgcrypt mingw-w64-ucrt-x86_64-libgcrypt mingw-w64-clang-x86_64-libgcrypt bison bisonc++ flex reflex flexc++ mingw-w64-ucrt-x86_64-flexdll mingw-w64-x86_64-flexdll mingw-w64-cross-binutils mingw-w64-ucrt-x86_64-binutils mingw-w64-x86_64-binutils binutils mingw-w64-clang-x86_64-lmdb mingw-w64-ucrt-x86_64-lmdb mingw-w64-x86_64-lmdb mingw-w64-clang-x86_64-lmdbxx mingw-w64-ucrt-x86_64-lmdbxx mingw-w64-x86_64-lmdbxx ninja mingw-w64-x86_64-arrow mingw-w64-ucrt-x86_64-arrow mingw-w64-x86_64-cmake pcre pcre2 mingw-w64-x86_64-pcre mingw-w64-x86_64-pcre2 mingw-w64-ucrt-x86_64-pcre mingw-w64-ucrt-x86_64-pcre2 mingw-w64-clang-x86_64-pcre mingw-w64-clang-x86_64-pcre2 libpcreposix libpcre2posix pcre2-devel pcre-devel libpcrecpp libpcre libpcre2_32 grep mingw-w64-x86_64-grep mingw-w64-clang-x86_64-grep mingw-w64-ucrt-x86_64-grep
    #)

#msys2-runtime-devel msys2-runtime  msys2-w32api-runtime

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm mingw-w64-x86_64-cmake
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm sqlite libsqlite cmake gcc gcc-libs pkgconf gnutls libgnutls liblzo2 mingw-w64-x86_64-pkgconf mingw-w64-x86_64-lzo2 mingw-w64-x86_64-gnutls mingw-w64-ucrt-x86_64-pkgconf mingw-w64-ucrt-x86_64-lzo2 mingw-w64-ucrt-x86_64-gnutls
    #)
    
    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm wget git
    #)
    
    execute_process(
        COMMAND nproc;
        OUTPUT_VARIABLE N_PROC
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND gcc -dumpversion || ${RTOOLS_DIR}/usr/bin/gcc -dumpversion;
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )

    set(ROOT_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/")
    message(STATUS "ROOT_SRC_DIR: ${ROOT_SRC_DIR}")	
    #if (NOT EXISTS "${ROOT_SRC_DIR}/R_src/")
    #    message(STATUS "Downloading R-src...")	
    #    execute_process(
	#		COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "mkdir -p ${ROOT_SRC_DIR}/R_src/" ;
	#		WORKING_DIRECTORY ${ROOT_SRC_DIR}
	#	)
	#	execute_process(
	#		COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "curl https://cran.r-project.org/src/base/R-latest.tar.gz -o R-latest.tar.gz" ;
	#		WORKING_DIRECTORY ${ROOT_SRC_DIR}/R_src/
	#		COMMAND_ERROR_IS_FATAL ANY
	#	)
	#	execute_process(
	#		COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "tar xzf R-latest.tar.gz 2> /dev/null" ;
	#		WORKING_DIRECTORY ${ROOT_SRC_DIR}/R_src/
	#	)
	#endif()

    #execute_process(
	#	# COMMAND bash -c "cygpath ${FOLDER_MODE} $(echo $(realpath $(find . -maxdepth 1 -iname 'R-*' -type d ))/src/include/)" ;
	#	COMMAND bash -c "cygpath ${FOLDER_MODE} $(echo $(realpath $(find . -maxdepth 1 -iname 'R-*' -type d ))/src/include/)" ;
	#	OUTPUT_VARIABLE R_INC_UNIX
    #	OUTPUT_STRIP_TRAILING_WHITESPACE
	#	WORKING_DIRECTORY ${ROOT_SRC_DIR}/R_src/
	#	COMMAND_ERROR_IS_FATAL ANY
	#)
	#message("R INC (SRC) : " ${R_INC_UNIX})

    ## Safely find the R source directory using native CMake
    #file(GLOB R_SRC_DIRS "${ROOT_SRC_DIR}/R_src/R-*")
    #if(R_SRC_DIRS)
    #    # Grab the first match
    #    list(GET R_SRC_DIRS 0 MATCHED_R_DIR)
    #    # Construct the absolute path
    #    set(R_INC_UNIX "${MATCHED_R_DIR}/src/include")
    #    message(STATUS "R INC (SRC) : ${R_INC_UNIX}")
    #else()
    #    message(FATAL_ERROR "Could not find R source directory in ${ROOT_SRC_DIR}/R_src/")
    #endif()
    #
    #if(R_SRC_DIRS)
    #    # Grab the first match
    #    list(GET R_SRC_DIRS 0 MATCHED_R_DIR)
    #    # Construct the absolute path
    #    set(R_ROOT_SRC_DIR "${MATCHED_R_DIR}/")
    #    message(STATUS "R ROOT (SRC) : ${R_ROOT_SRC_DIR}")
    #else()
    #    message(FATAL_ERROR "Could not find R source directory in ${ROOT_SRC_DIR}/R_src/")
    #endif()

    if (NOT EXISTS "${ROOT_SRC_DIR}/arrow/")
    execute_process(
        COMMAND git clone https://github.com/apache/arrow arrow
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        #COMMAND_ERROR_IS_FATAL ANY    
    )
    endif()

    #execute_process(
    #    COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r -u -h -k3 -t '-' | grep -iv \"dev\" | head -n 1" ; #no dev versions
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow
    #    OUTPUT_VARIABLE ARROW_GIT_TAG
    #	OUTPUT_STRIP_TRAILING_WHITESPACE
    #)

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l 'apache-arrow-*.*.*' --sort=-v:refname | grep -iv \"dev\" | grep -iv \"rc\" | grep -iv \"js\" | head -n 1"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow
        OUTPUT_VARIABLE ARROW_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    message("ARROW GIT LATEST TAG : " ${ARROW_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${ARROW_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow    
        COMMAND_ERROR_IS_FATAL ANY
    )

    if (NOT EXISTS "${ROOT_SRC_DIR}/ncbi-cxx-toolkit-public/")
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND_ERROR_IS_FATAL ANY
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l --sort=-v:refname | grep -iv \"dev\" | head -n 1" ; #no dev versions
        #COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r | head -n 1" ;
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public
        OUTPUT_VARIABLE NCBI_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("NCBI GIT LATEST TAG : " ${NCBI_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${NCBI_GIT_TAG}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public 
        COMMAND_ERROR_IS_FATAL ANY
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/src/ " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/QuickBLAST.prj 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_crosscomp_inject_pre.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_crosscomp_inject_post.sh 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/dependencies/ncbi_components.cmake 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    file(COPY ${CMAKE_SOURCE_DIR}/src/ 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/src/
    )

    message("Copying " ${CMAKE_SOURCE_DIR}/src/include " into " ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/)

    file(COPY ${CMAKE_SOURCE_DIR}/src/include
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    )

    execute_process(
        COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/cmake-configure ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/configure
    )

    set(CMAKE_CXX_EXTENSIONS ON )
    set(CMAKE_SYSROOT_LINK "/x86_64-w64-mingw32.static.posix/lib")
    set(CMAKE_SYSROOT_LINK "/x86_64-w64-mingw32.static.posix/lib")
    set(CMAKE_SYSROOT_COMPILE "/x86_64-w64-mingw32.static.posix/include")
	set(CMAKE_SYSROOT "/x86_64-w64-mingw32.static.posix")

    set(FOLDER_MODE "-m") #FOR cygpath
	function(cygpath_cmd ret_var _folder_path)
		block(SCOPE_FOR VARIABLES)
			# message(STATUS ${_folder_path})
			execute_process(
				COMMAND bash -c "cygpath ${FOLDER_MODE} '${_folder_path}'" ;
				OUTPUT_VARIABLE _return_path
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)
			# message(STATUS ${ret_var} ${_return_path})
			set(${ret_var} ${_return_path})
			return(PROPAGATE ${ret_var} ${_return_path})
		endblock()
	endfunction(cygpath_cmd)

    #TRYING with NATIVE CMake calls instead of calling it through the Makevars.ucrt shell
    include(ExternalProject)

    # Generate a native CMake script to safely prepend the C++20 polyfill
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/patch_arrow.cmake" [=[
    file(READ "${ARROW_FILE}" FILE_CONTENT)
    set(POLYFILL [==[
    #include <format>
    #include <chrono>
    #include <string_view>
    
    namespace arrow { namespace internal { struct OffsetZone; } }
    
    // 1. Catch Arrow's custom OffsetZone
    template <typename D> struct std::formatter<std::chrono::zoned_time<D, const arrow::internal::OffsetZone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    
    // 2. Catch GCC 14's failing standard time_zone with ratio<86400> (Days)
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::ratio<86400>>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    
    // 3. Catch GCC 14's failing standard time_zone with ratio<1> (Seconds)
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::ratio<1>>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    
    // 4. Catch millisecond precision
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::milli>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    
    // 5. Catch microsecond precision
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::micro>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    
    // 6. Catch nanosecond precision
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::nano>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, char> { 
        template <typename C> auto format(const auto&, C& ctx) const { return std::formatter<std::string_view, char>::format("ArrowTime", ctx); } 
    };
    ]==])
    file(WRITE "${ARROW_FILE}" "${POLYFILL}\n${FILE_CONTENT}")
    ]=])

    ExternalProject_Add(build_arrow
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/arrow_ep"
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp"
        BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2"
        
        PATCH_COMMAND ${CMAKE_COMMAND} -DARROW_FILE=${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/src/arrow/compute/kernels/scalar_temporal_unary.cc -P "${CMAKE_CURRENT_BINARY_DIR}/patch_arrow.cmake"
        
        CMAKE_GENERATOR "MSYS Makefiles"
        
        # Pass all your definitions here
        CMAKE_ARGS
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DARROW_BUILD_INTEGRATION=OFF
            -DARROW_BUILD_TESTS=OFF
            -DARROW_DEPENDENCY_USE_SHARED=OFF
            -DARROW_WITH_BZ2=ON
            -DARROW_WITH_LZ4=ON
            -DARROW_BUILD_SHARED=OFF
            -DPARQUET_BUILD_SHARED=OFF
            -DPARQUET_BUILD_STATIC=ON  
            -DARROW_BUILD_STATIC=ON
            -DBUILD_SHARED_LIBS=OFF
            -DBUILD_STATIC_LIBS=ON
            -DPARQUET_MINIMAL_DEPENDENCY=OFF
            -DARROW_WITH_BROTLI=ON
            -DARROW_WITH_RE2=OFF
            -DARROW_WITH_ZLIB=ON
            -DARROW_WITH_ZSTD=ON
            -DARROW_IPC=ON
            -DARROW_FILESYSTEM=ON
            -DARROW_COMPUTE=ON
            -DARROW_ENABLE_THREADING=ON
            -DARROW_MIMALLOC=OFF
            -DARROW_JEMALLOC=OFF
            -DARROW_DATASET=OFF
            -DARROW_ACERO=OFF
            -DARROW_PARQUET=ON
            #-DARROW_WITH_THRIFT=ON
            #-DThrift_SOURCE=BUNDLED
            #-DTHRIFT_BUILD_COMPILER=OFF
            -DARROW_CSV=ON
            -DARROW_HDFS=ON
            -DARROW_VERBOSE_THIRDPARTY_BUILD=ON
            -DARROW_POSITION_INDEPENDENT_CODE=ON
            -DARROW_JSON=ON
            -DARROW_SUBSTRAIT=ON
            -DARROW_WITH_SNAPPY=ON
            -DARROW_WITH_UTF8PROC=ON
            -DARROW_RPATH_ORIGIN=ON
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
            -DCMAKE_C_COMPILER=gcc.exe
            -DCMAKE_CXX_COMPILER=g++.exe
            "-DCMAKE_C_FLAGS= -pedantic -Wformat -Werror=format-security -UHAVE_AF_UNIX_H -D__MINGW32__ -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc -mtune=generic -maccumulate-outgoing-args -pthread -mconsole -mstackrealign "
            "-DCMAKE_CXX_FLAGS= -pedantic -Wformat -Werror=format-security -UHAVE_AF_UNIX_H -D__MINGW32__ -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -mtune=generic -march=native -maccumulate-outgoing-args -pthread -mconsole -mstackrealign"
            "-DCMAKE_SHARED_LINKER_FLAGS=-Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,-shared -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,--dll -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--no-whole-archive -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group "
        
        #-Wall
        #-DARROW_EXPORTING
        #-std=gnu2x
        #-std=gnu++20
        #-O3 -fPIC  -fdiagnostics-color=always

        # Map the Makevars build command
        BUILD_COMMAND $(MAKE) -j ${N_PROC} 
        #arrow parquet
        
        # Map the Makevars copy/install command
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow.a"
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow_msys2.a"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/release/libparquet.a"
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/release/libparquet_msys2.a"
    )

    # =========================================================================
    # AUTOMATE CURL DISCOVERY & WINDOWS DEPENDENCY EXTRACTION
    # =========================================================================
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CURL_PKG REQUIRED libcurl)

    # 1. Find the absolute path to the main curl library
    find_library(AUTO_CURL_LIB 
        NAMES libcurl.a curl
        PATHS ${CURL_PKG_LIBRARY_DIRS}
        NO_DEFAULT_PATH
    )

    # 2. Extract the include directory
    list(GET CURL_PKG_INCLUDE_DIRS 0 AUTO_CURL_INC)

    # 3. Extract dependencies (ws2_32, crypt32, etc.) and format as "-lws2_32 -lcrypt32"
    set(AUTO_CURL_DEPS "")
    foreach(_lib IN LISTS CURL_PKG_LIBRARIES)
        # Skip the main library, we only want the dependencies
        if(NOT _lib STREQUAL "curl")
            string(APPEND AUTO_CURL_DEPS " -l${_lib}")
        endif()
    endforeach()

    message(STATUS "Automated CURL Library: ${AUTO_CURL_LIB}")
    message(STATUS "Automated CURL Include: ${AUTO_CURL_INC}")
    message(STATUS "Automated CURL Deps: ${AUTO_CURL_DEPS}")

    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/sqlite/libsqlite3.a")    
        file(DOWNLOAD
            "https://sqlite.org/download.html"
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/download.html"
            SHOW_PROGRESS
            STATUS download_status
            TIMEOUT 10
        )
        list(GET download_status 0 download_code)
        if(NOT "${download_code}" STREQUAL "0")
            message(FATAL_ERROR "SQLite version check failed: ${download_status}")
        endif()
        execute_process(
            COMMAND bash -c "grep -Eo '[0-9]{4}/sqlite-amalgamation-[0-9]+\\\\.zip' ${CMAKE_CURRENT_SOURCE_DIR}/BUILD/download.html | head -n 1"
            OUTPUT_VARIABLE SQLITE_LINK
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if (SQLITE_LINK)
            message(STATUS "SQLite version link found: https://sqlite.org/${SQLITE_LINK}")
        else()
            message(FATAL_ERROR "SQLite version link not found")
        endif()
        execute_process(
            COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/BUILD/sqlite/"
            )
        file(DOWNLOAD
            "https://sqlite.org/${SQLITE_LINK}"
            "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/sqlite/sqlite.zip"
            SHOW_PROGRESS
            STATUS download_status
            TIMEOUT 10
        )
        list(GET download_status 0 download_code)
        if(NOT "${download_code}" STREQUAL "0")
            message(FATAL_ERROR "SQLite download failed: ${download_status}")
        endif()
    endif()
    # =========================================================================

    # Setup the massive CMAKE_ARGS string for the NCBI script
    #--log-level=VERBOSE 
    set(NCBI_CMAKE_ARGS " \
        --fresh \
        -G 'MSYS Makefiles' \
        -DCMAKE_MAKE_PROGRAM=make \
        -DCMAKE_SYSTEM_NAME=Windows \
        -DMINGW=1 \
        -DWIN32=1 \
        -DUNIX=0 \
        -DCYGWIN=0 \
        -DMSYS=0 \
        -DNCBI_COMPILER_MSVC=0 \
        -DNCBI_COMPILER_GCC=1 \
        -DCMAKE_CXX_COMPILER_VERSION=${GCC_VER} \
        -DSHLIB_OPENMP_CXXFLAGS='${SHLIB_OPENMP_CXXFLAGS}' \
        -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
        -DROOT_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR} \
        -DNCBITK_ROOT_DIR=${CMAKE_CURRENT_SOURCE_DIR}/BUILD/ncbi-cxx-toolkit-public \
        -DR_HOME='$ENV{R_HOME}' \
        -DR_ARCH='$ENV{R_ARCH}' \
        -DRCPP_FOLDER='${RCPP_FOLDER}' \
        -DRCPPPROG_FOLDER='${RCPPPROG_FOLDER}' \
        -DRCPPTHREAD_FOLDER='${RCPPTHREAD_FOLDER}' \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_STATIC_LIBS=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY \
        -DNCBI_COMPONENT_CURL_FOUND=YES \
        -DNCBI_COMPONENT_CURL_LIBS='${AUTO_CURL_LIB}' \
        -DNCBI_COMPONENT_CURL_INCLUDE='${AUTO_CURL_INC}' \
        -DCMAKE_C_STANDARD_LIBRARIES='${AUTO_CURL_DEPS}' \
        -DCMAKE_CXX_STANDARD_LIBRARIES='${AUTO_CURL_DEPS}' \
        -DCMAKE_C_FLAGS='${CMAKE_C_FLAGS} -pedantic -Wno-attributes -Wno-date-time -Wno-unused-parameter -Wno-ignored-attributes -Wformat -Werror=format-security -Wa,-mbig-obj -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -static-libgcc -mtune=generic -maccumulate-outgoing-args -pthread -mconsole -mstackrealign ' \
        -DCMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS} -DLIBICONV_STATIC -pedantic -Wno-attributes -Wno-date-time -Wno-unused-parameter -Wno-ignored-attributes -Wformat -Werror=format-security -Wa,-mbig-obj -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -static-libstdc++ -mtune=generic -march=native -maccumulate-outgoing-args -pthread -mconsole -mstackrealign ' \
        -DCMAKE_SHARED_LINKER_FLAGS=' -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,-shared -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,--dll -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--no-whole-archive -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group ' \
        -DCMAKE_EXE_LINKER_FLAGS=' -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group ' \
        -DCMAKE_USE_WIN32_THREADS_INIT=1 \
        -DCMAKE_USE_PTHREADS_INIT=0 \
        -DTHREADS_PREFER_PTHREAD_FLAG=OFF \
    ")
    #-ftrack-macro-expansion=0 -pipe 
    #-std=gnu2x
    #-std=gnu++20
    #-O3 -fPIC

    ExternalProject_Add(ncbi_toolkit
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/ncbi_ep"
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/ncbi-cxx-toolkit-public"
        BUILD_IN_SOURCE TRUE
        
        # STEP 1: Run pre-inject and cmake-configure
        CONFIGURE_COMMAND bash -c "SHLIB_OPENMP_CXXFLAGS='${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}' RTOOLS_DIR='${RTOOLS_DIR}' ./ncbi_crosscomp_inject_pre.sh && export CMAKE_ARGS=\"${NCBI_CMAKE_ARGS}\" && export CMAKE='${RTOOLS_DIR}/ucrt64/bin/cmake.exe' && ./cmake-configure GCC ${GCC_VER} --with-projects=QuickBLAST.prj --without-dll --with-composite --without-debug --with-build-root=BUILD"
        #--with-dll
        # STEP 2: Use double quotes, escape inner quotes, and use $$ to pass command substitution safely to Bash
        BUILD_COMMAND $(MAKE) -C BUILD/build -j ${N_PROC} datatool all makeblastdb QuickBLASTcpp
        
        # STEP 3: Run the post-injection script
        INSTALL_COMMAND bash -c "RTOOLS_DIR='${RTOOLS_DIR}' ./ncbi_crosscomp_inject_post.sh"
    )

    execute_process(COMMAND Rscript -e "cat(Sys.getenv('R_HOME'))" OUTPUT_VARIABLE R_USR_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(R_INC "${R_USR_HOME}/include")
    cygpath_cmd(R_INC "${R_USR_HOME}/include")
    message("R inc:" ${R_INC})
    
    execute_process(COMMAND Rscript -e "cat(file.path(Sys.getenv('R_HOME'),'bin',Sys.getenv('R_ARCH')))" OUTPUT_VARIABLE R_BIN OUTPUT_STRIP_TRAILING_WHITESPACE)
    cygpath_cmd(R_BIN "${R_BIN}")
    message("R BIN:" ${R_BIN})

    #ExternalProject_Add(configure_R
    #    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/configure_r_ep"
    #    SOURCE_DIR "${R_ROOT_SRC_DIR}"
    #    BUILD_IN_SOURCE TRUE
    #    
    #    # STEP 1: Run pre-inject and cmake-configure
    #    CONFIGURE_COMMAND bash -c "SHLIB_OPENMP_CXXFLAGS='${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}' RTOOLS_DIR='${RTOOLS_DIR}' ./configure --enable-R-static-lib --enable-lto --enable-static --with-gnu-ld --with-internal-tzcode --without-x --disable-rpath "
    #
    #    # STEP 2: Use double quotes, escape inner quotes, and use $$ to pass command substitution safely to Bash
    #    BUILD_COMMAND bash -c "echo ${R_INC}..."
    #    
    #    # STEP 3: Run the post-injection script
    #    INSTALL_COMMAND bash -c "echo Skipping R INSTALL..."
    #)

    # Ensure Arrow builds before NCBI
    add_dependencies(ncbi_toolkit build_arrow)
    #add_dependencies(ncbi_toolkit_static ncbi_toolkit)

    # =========================================================================
    # 3. DEFINE THE MAIN QUICKBLAST TARGET
    # =========================================================================

    # Gather your main Rcpp / C++ source files
    file(GLOB QUICKBLAST_SOURCES "src/*.cpp" "src/*.c")
    
    execute_process(COMMAND Rscript -e "cat(system.file(package = 'Rcpp'))" OUTPUT_VARIABLE RCPP_FOLDER OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND Rscript -e "cat(system.file(package = 'RcppProgress'))" OUTPUT_VARIABLE RCPPPROG_FOLDER OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND Rscript -e "cat(system.file(package = 'RcppThread'))" OUTPUT_VARIABLE RCPPTHREAD_FOLDER OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Create the shared library target for R
    #add_library(QuickBLAST STATIC ${QUICKBLAST_SOURCES})
    add_library(QuickBLAST SHARED ${QUICKBLAST_SOURCES})

    add_dependencies(QuickBLAST ncbi_toolkit)
    #add_dependencies(QuickBLAST ncbi_toolkit_static)
    #add_dependencies(QuickBLAST configure_R)
    
    # =========================================================================
    # 4. LINK ARROW (IMPORTED TARGET)
    # =========================================================================

    # Create a "fake" target for Arrow so CMake knows how to link it
    add_library(Arrow::msys2 STATIC IMPORTED GLOBAL)
    set_target_properties(Arrow::msys2 PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow_msys2.a"
    )
    # Tell CMake: "Do not try to link Arrow until build_arrow is 100% finished"
    add_dependencies(Arrow::msys2 build_arrow)

    # =========================================================================
    # 5. LINK NCBI (IMPORTED TARGETS)
    # =========================================================================

    # Example: Linking the core NCBI static libraries (update names based on what NCBI outputs)
    add_library(NCBI::core STATIC IMPORTED GLOBAL)
    set_target_properties(NCBI::core PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/ncbi-cxx-toolkit-public/BUILD/lib/libncbi_core.a"
    )
    # Tell CMake: "Do not try to link NCBI until setup_ncbi is 100% finished"
    add_dependencies(NCBI::core setup_ncbi)
    
    # =========================================================================
    # FORCE QUICKBLAST TO MATCH NCBI TOOLKIT'S WINDOWS ABI
    # =========================================================================
    target_compile_definitions(QuickBLAST PUBLIC
        NCBI_OS_MSWIN=1
        WIN32=1
        _WIN32=1
        UNICODE
        _UNICODE
        NCBI_WIN32_THREADS=1
    )
    target_compile_features(QuickBLAST PRIVATE cxx_std_20)
    # Tell the compiler that iconv is a static library
    target_compile_definitions(QuickBLAST PUBLIC LIBICONV_STATIC)
    target_compile_definitions(QuickBLAST PRIVATE ARROW_STATIC)
    target_compile_definitions(QuickBLAST PRIVATE PARQUET_STATIC)

    # Forcefully remove POSIX thread macros if CMake's FindThreads injected them into QuickBLAST
    get_target_property(QB_DEFS QuickBLAST COMPILE_DEFINITIONS)
    if(QB_DEFS)
        list(REMOVE_ITEM QB_DEFS "NCBI_POSIX_THREADS" "NCBI_USE_MUTEX_PTHREADS")
        set_target_properties(QuickBLAST PROPERTIES COMPILE_DEFINITIONS "${QB_DEFS}")
    endif()

    set_target_properties(QuickBLAST PROPERTIES 
        PREFIX "" 
        OUTPUT_NAME "QuickBLAST"
        SUFFIX ".dll"
    )
    # =========================================================================
    # 6. BRING IT ALL TOGETHER
    # =========================================================================

    # Apply all include directories needed for your QuickBLAST code
    target_include_directories(QuickBLAST PRIVATE
        #"${R_INC_UNIX}"
        "${R_INC}"
        "${R_INC}/R_ext/"
        #"${R_ROOT_SRC_DIR}/src/gnuwin32/fixed/h/"
        "${RCPP_FOLDER}/include"
        "${RCPPPROG_FOLDER}/include"
        "${RCPPTHREAD_FOLDER}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/src"
        "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/arrow/cpp/BUILD_MSYS2/src"
        "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/ncbi-cxx-toolkit-public/include" # To find QuickBLAST/commons.hpp
    )

    # 1. Tell CMake where the directory of the built libraries is
    target_link_directories(QuickBLAST PRIVATE
        "${CMAKE_CURRENT_BINARY_DIR}/sqlite"
        "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/bin"
        "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/bin"
        "${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD_MSYS2/release"
        "${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD_MSYS2/lib"
        "${R_BIN}"
    )

    set(ICONV_STATIC_LIB "${RTOOLS_DIR}/x86_64-w64-mingw32.static.posix/lib/libiconv.a")

    # 2. Link the libraries by name (do not include "lib" or ".a")
    target_link_libraries(QuickBLAST PUBLIC
        -Wl,--enable-auto-import
        -Wl,--start-group
        bz2
        arrow
        arrow_bundled_dependencies
        parquet
        thrift
        "${CMAKE_CURRENT_BINARY_DIR}/sqlite/libsqlite3.a"
        NCBI::core
        ncbi_algo
        ncbi_align_format
        ncbi_blastinput
        ncbi_dbapi_driver
        ncbi_eutils
        ncbi_general
        ncbi_misc
        ncbi_pub
        ncbi_seq
        ncbi_seqext
        ncbi_trackmgr
        ncbi_web
        ncbi_xdbapi_ftds
        ncbi_xloader_blastdb
        ncbi_xloader_genbank
        ncbi_xobjsimple
        ncbi_xreader_cache
        ncbi_xreader_id1
        ncbi_xreader_id2
        ncbi_xreader
        xxconnect2
        ct_ftds14
        efetch
        lmdb
        ncbi_xloader_blastdb_rmt
        psg_client
        QuickBLASTcpp
        sqlitewrapp
        tds_ftds14
        utrtprof
        xflatfile
        zcf
        zstd
        R
        #Required Native Windows System Libraries
        dbghelp
        advapi32
        psapi
        iphlpapi
        oleaut32
        ole32
        user32
        netapi32
        ws2_32      # Networking
        wsock32     # Networking legacy
        crypt32     # Cryptography
        secur32     # Security
        normaliz    # Unicode normalization
        bcrypt
        sqlite3
        z
        lz4
        brotlidec
        brotlienc
        brotlicommon
        iconv
        #${ICONV_STATIC_LIB}
        -Wl,--end-group
    )

    # 1. Identify the R architecture directory
    if(DEFINED ENV{R_ARCH})
        set(R_ARCH_DIR "/$ENV{R_ARCH}")
    else()
        set(R_ARCH_DIR "")
    endif()

    # Destination inside your R package source
    set(DEV_LIB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/libs${R_ARCH}")
    set(INST_LIB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/inst/libs${R_ARCH}")

    add_custom_command(TARGET QuickBLAST POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${DEV_LIB_DEST}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INST_LIB_DEST}"
    )

    # Source paths (adjust these if your BUILD directories differ)
    set(NCBI_BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/bin")
    set(ARROW_BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD_MSYS2/release")

    ## 2. List of specific NCBI files to strip and copy
    set(NCBI_FILES
    #    "libncbi_algo.dll"
    #    "libncbi_align_format.dll"
    #    "libncbi_blastinput.dll"
    #    "libncbi_core.dll"
    #    "libncbi_dbapi_driver.dll"
    #    "libncbi_eutils.dll"
    #    "libncbi_general.dll"
    #    "libncbi_misc.dll"
    #    "libncbi_pub.dll"
    #    "libncbi_seq.dll"
    #    "libncbi_seqext.dll"
    #    "libncbi_trackmgr.dll"
    #    "libncbi_web.dll"
    #    "libncbi_xdbapi_ftds.dll"
    #    "libncbi_xloader_blastdb.dll"
    #    "libncbi_xloader_genbank.dll"
    #    "libncbi_xobjsimple.dll"
    #    "libncbi_xreader.dll"
    #    "libncbi_xreader_cache.dll"
    #    "libncbi_xreader_id1.dll"
    #    "libncbi_xreader_id2.dll"
    #    "libxxconnect2.dll"
        "makeblastdb.exe"
    )

    # --- 1. Define Destinations (Make sure these variable names match exactly) ---
    set(DEV_LIB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/libs${R_ARCH}")
    set(INST_LIB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/inst/libs${R_ARCH}")

    # --- 2. Create BOTH directories before copying ---
    add_custom_command(TARGET QuickBLAST POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${DEV_LIB_DEST}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INST_LIB_DEST}"
        COMMENT "Syncing QuickBLAST and NCBI libraries to R package directories..."
    )

    # --- 3. Process NCBI files ---
    foreach(FILENAME ${NCBI_FILES})
        add_custom_command(TARGET QuickBLAST POST_BUILD
            # Copy to the root libs/ folder for devtools::load_all()
            COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                    "${NCBI_BIN_DIR}/${FILENAME}" "${DEV_LIB_DEST}/${FILENAME}"
                    
            # Copy to inst/libs/ for the final package installation
            COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                    "${NCBI_BIN_DIR}/${FILENAME}" "${INST_LIB_DEST}/${FILENAME}"
        )
    endforeach()

    ## Process Arrow DLLs
    ## Note: Since CMake's copy command doesn't support wildcards directly,
    ## we use the system shell provided by Rtools/MinGW
    #add_custom_command(TARGET QuickBLAST POST_BUILD
    #    COMMAND bash -c "cp -uf ${ARROW_BIN_DIR}/*.dll ${INST_LIBS_DEST}/"
    #    COMMENT "Deploying Arrow dependencies"
    #)

endif()

