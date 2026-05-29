# $Id$

#Apache Arrow requires a minimum CMake version of 3.25
cmake_minimum_required(VERSION 3.25) 
cmake_policy(VERSION 3.25)

set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)
set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF)
set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF)

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
    
    if (NOT EXISTS "${ROOT_SRC_DIR}/libuv/")
    execute_process(
    COMMAND git clone https://github.com/libuv/libuv libuv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND_ERROR_IS_FATAL ANY
    )
    endif()
    execute_process(
  	  COMMAND bash -c "realpath ../../../libuv/include/" ;
       OUTPUT_VARIABLE UV_INC
       OUTPUT_STRIP_TRAILING_WHITESPACE
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
  
  	#execute_process(
    #  COMMAND pkg-config --cflags libR;
    #  OUTPUT_VARIABLE R_FLAGS
    #  OUTPUT_STRIP_TRAILING_WHITESPACE
  	#)
  	#message("R flags:" ${R_FLAGS})

#-fdiagnostics-color=always -Wall
#-std=gnu++20
#-std=gnu2x
#-O3 -fPIC  -fdiagnostics-color=always

if(APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin" OR CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    # Apple Clang (macOS) compatible flags
    #${R_FLAGS}
    set(ARROW_CXX_FLAGS " -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -pthread")
    set(ARROW_C_FLAGS  " -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -pthread")
else()
    # GCC compatible flags
    #${R_FLAGS} -static-libgcc  -static-libstdc++ 
    set(ARROW_CXX_FLAGS " -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -mtune=generic -march=native -maccumulate-outgoi")
    set(ARROW_C_FLAGS  " -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native -mtune=generic -maccumulate-outgoin")
endif()

set(ARROW_EXE_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")
set(ARROW_LINKER_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")


    execute_process(
    COMMAND cmake .. -G "Unix Makefiles" 
    "-DCMAKE_CXX_FLAGS=${ARROW_CXX_FLAGS}"
    "-DCMAKE_C_FLAGS=${ARROW_C_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${ARROW_EXE_FLAGS}"
    "-DCMAKE_SHARED_LINKER_FLAGS=${ARROW_LINKER_FLAGS}"
    -DARROW_DEPENDENCY_USE_SHARED=OFF
    -DARROW_BUILD_SHARED=OFF
    -DARROW_BUILD_STATIC=ON
    -DPARQUET_BUILD_SHARED=OFF
    -DPARQUET_BUILD_STATIC=ON
    -DARROW_BZ2_USE_SHARED=OFF
    -DARROW_COMPUTE=ON
    -DARROW_PYTHON=OFF
    -DARROW_CSV=ON
    -DARROW_DEPENDENCY_SOURCE=BUNDLED
    -DARROW_ENABLE_THREADING=ON
    -DARROW_FILESYSTEM=ON
    -DARROW_HDFS=ON
    -DARROW_IPC=ON
    -DARROW_JEMALLOC=ON
    -DARROW_JSON=ON
    -DARROW_PARQUET=ON
    -DARROW_POSITION_INDEPENDENT_CODE=ON
    -DARROW_TENSORFLOW=ON
    -DARROW_VERBOSE_THIRDPARTY_BUILD=ON
    -DARROW_WITH_BACKTRACE=ON
    -DARROW_WITH_BROTLI=ON
    -DARROW_WITH_BZ2=ON
    -DARROW_WITH_LZ4=ON
    -DARROW_WITH_RE2=ON
    -DARROW_RPATH_ORIGIN=ON
    -DARROW_WITH_SNAPPY=ON
    -DARROW_WITH_UTF8PROC=ON
    -DARROW_WITH_ZLIB=OFF
    -DARROW_WITH_ZSTD=ON
    -DCMAKE_BUILD_TYPE=release
    -DBUILD_SHARED_LIBS=OFF
    -DBUILD_STATIC_LIBS=ON
    -DCMAKE_SKIP_BUILD_RPATH=FALSE
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE
    -DCMAKE_INSTALL_RPATH="$ORIGIN"
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DARROW_BROTLI_USE_SHARED=OFF
    -DARROW_LZ4_USE_SHARED=OFF
    -DARROW_ZSTD_USE_SHARED=OFF
    -DARROW_LLVM_USE_SHARED=OFF
    -DARROW_LZ4_USE_SHARED=OFF
    -DARROW_OPENSSL_USE_SHARED=OFF
    -DARROW_THRIFT_USE_SHARED=OFF
    -DARROW_UTF8PROC_USE_SHARED=OFF
    -DARROW_BOOST_USE_SHARED=OFF
    -DARROW_BROTLI_USE_SHARED=OFF
    -DARROW_PROTOBUF_USE_SHARED=OFF
    -DARROW_SNAPPY_USE_SHARED=OFF 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    )

#-DARROW_ZSTD_USE_SHARED=OFF -DARROW_JEMALLOC_USE_SHARED=OFF # -DARROW_CUDA=ON -DARROW_DEPENDENCY_SOURCE=BUNDLED  # -DARROW_BUILD_INTEGRATION=ON

    execute_process(
    COMMAND make -j ${N_PROC} 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    )

    # ==============================================================================
    # SYNCHRONOUS SQLITE PRE-BUILD (Forces SQLite to exist before NCBI configures)
    # ==============================================================================
    set(SQLITE_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/BUILD/sqlite")
    set(SQLITE_TARGET_LIB "${SQLITE_ROOT_DIR}/libsqlite3.a")
    
    if(NOT EXISTS "${SQLITE_TARGET_LIB}")
        message(STATUS "Pre-building SQLite3 statically during CMake configuration...")
        
        file(MAKE_DIRECTORY "${SQLITE_ROOT_DIR}")
        
        # 1. Download and Extract
        file(DOWNLOAD "https://sqlite.org/download.html" "${SQLITE_ROOT_DIR}/download.html" STATUS dl_status TIMEOUT 10)
        execute_process(COMMAND bash -c "grep -Eo '[0-9]{4}/sqlite-amalgamation-[0-9]+\\.zip' ${SQLITE_ROOT_DIR}/download.html | head -n 1" OUTPUT_VARIABLE SQLITE_LINK OUTPUT_STRIP_TRAILING_WHITE)
        file(DOWNLOAD "https://sqlite.org/${SQLITE_LINK}" "${SQLITE_ROOT_DIR}/sqlite.zip" STATUS dl_status TIMEOUT 10)
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "sqlite.zip" WORKING_DIRECTORY "${SQLITE_ROOT_DIR}")
        file(GLOB SQLITE_SRC_DIR "${SQLITE_ROOT_DIR}/sqlite-amalgamation-*")
        
        # 2. Write a minimal, dedicated CMakeLists.txt just for SQLite
        file(WRITE "${SQLITE_SRC_DIR}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.10)
    project(sqlite3_prebuild C)
    add_library(sqlite3 STATIC sqlite3.c)
    # Ensure the output is strictly named libsqlite3.a even on Windows/MinGW
    set_target_properties(sqlite3 PROPERTIES OUTPUT_NAME \"sqlite3\" PREFIX \"lib\")
    target_compile_definitions(sqlite3 PRIVATE _FORTIFY_SOURCE=2 SQLITE_ENABLE_UNLOCK_NOTIFY=1 SQLITE_ENABLE_MEMSYS5 SQLITE_ENABLE_MEMSYS3 SQLITE_ENABLE_FTS3 SQLITE_ENABLE_FTS4 SQLITE_ENABLE_FTS5)
    # Safe, cross-platform flags
    target_compile_options(sqlite3 PRIVATE -Wformat -Werror=format-security -pedantic -pthread)
    if(WIN32)
        target_compile_options(sqlite3 PRIVATE -Wno-format-y2k -static-libgcc -mconsole -mstackrealign)
    endif()
    ")
    
        # 3. SYNCHRONOUSLY configure and build SQLite using CMake
        execute_process(
            COMMAND ${CMAKE_COMMAND} -S . -B build_dir 
            WORKING_DIRECTORY "${SQLITE_SRC_DIR}"
            COMMAND_ERROR_IS_FATAL ANY
        )
        execute_process(
            COMMAND ${CMAKE_COMMAND} --build build_dir
            WORKING_DIRECTORY "${SQLITE_SRC_DIR}"
            COMMAND_ERROR_IS_FATAL ANY
        )
    
        # 4. Copy the freshly compiled static library to the exact location ncbi_components.cmake expects!
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E copy "${SQLITE_SRC_DIR}/build_dir/libsqlite3.a" "${SQLITE_TARGET_LIB}"
        )
        
        message(STATUS "SQLite3 pre-build complete! Library secured at: ${SQLITE_TARGET_LIB}")
    endif()
    # ==============================================================================
    file(GLOB SQLITE_INC_DIR "${SQLITE_ROOT_DIR}/sqlite-amalgamation-*")
    
#-std=gnu2x
#-std=gnu++20
#-O3 -fPIC  -fdiagnostics-color=always

if(APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin" OR CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
set(ENV{CMAKE_ARGS} "-DRCPP_FOLDER=${RCPP_FOLDER} \
    -DCMAKE_C_FLAGS=' -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native' \
    -DCMAKE_CXX_FLAGS=' -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native' \
    -DRCPPPROG_FOLDER=${RCPPPROG_FOLDER} \
    -DRCPPTHREAD_FOLDER=${RCPPTHREAD_FOLDER} \
    -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
    -DR_INCLUDE_DIR='${R_INCLUDE_DIR}' \
    -DCMAKE_BUILD_RPATH='-Wl,-rpath,@loader_path' \
    -DCMAKE_INSTALL_RPATH='-Wl,-rpath,@loader_path' \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF")
else()
    set(ENV{CMAKE_ARGS} "-DRCPP_FOLDER=${RCPP_FOLDER} \
    -DCMAKE_C_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native ' \
    -DCMAKE_CXX_FLAGS=' -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -march=native' \
    -DRCPPPROG_FOLDER=${RCPPPROG_FOLDER} \
    -DRCPPTHREAD_FOLDER=${RCPPTHREAD_FOLDER} \
    -DTOOLKIT_CPPFLAGS='${TOOLKIT_CPPFLAGS}' \
    -DR_INCLUDE_DIR='${R_INCLUDE_DIR}' \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_BUILD_RPATH='\$$ORIGIN' \
    -DCMAKE_INSTALL_RPATH='\$$ORIGIN' \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF")
endif()
    #-static-libstdc++   -static-libgcc
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
    
    #if(APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin" OR CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    #  execute_process(
    #    COMMAND ./configure --with-projects=QuickBLAST.prj --without-debug --with-build-root=BUILD --without-dll --with-sqlite3="${SQLITE_ROOT_DIR}"
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
    #    COMMAND_ERROR_IS_FATAL ANY
    #  )
    #  execute_process(
    #    COMMAND make -j ${N_PROC} all_p
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build
    #    COMMAND_ERROR_IS_FATAL ANY
    #  )
    #else()
      execute_process(
        COMMAND ${CMAKE_COMMAND} -E env 
          "CC=${CMAKE_C_COMPILER}" 
          "CXX=${CMAKE_CXX_COMPILER}" 
          "CMAKE_PREFIX_PATH=${SQLITE_ROOT_DIR}"
          "CMAKE_INCLUDE_PATH=${SQLITE_INC_DIR};${UV_INC}"
          "CMAKE_LIBRARY_PATH=${SQLITE_ROOT_DIR}"
      ./cmake-configure 
          --with-projects=QuickBLAST.prj 
          --without-debug 
          --with-build-root=BUILD 
          --without-dll
          "-DSQLITE3_LIBRARY=${SQLITE_TARGET_LIB}"
          "-DSQLITE3_INCLUDE_DIR=${SQLITE_INC_DIR}"
          "-DNCBI_COMPONENT_SQLITE3_FOUND=YES"
          "-DNCBI_COMPONENT_SQLITE3_INCLUDE=${SQLITE_INC_DIR}"
          "-DNCBI_COMPONENT_SQLITE3_LIBS=${SQLITE_TARGET_LIB}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/
            COMMAND_ERROR_IS_FATAL ANY
        )
    #endif()
    execute_process(
          COMMAND make -j ${N_PROC} datatool QuickBLASTcpp makeblastdb
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build
          COMMAND_ERROR_IS_FATAL ANY
        )

#    set(TOOLCHAIN_FILE "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build/ncbi_toolchain.cmake")

    
elseif(WIN32)
    message(STATUS "Rtools Install Path:" ${RTOOLS_DIR})

    #execute_process(
    #    COMMAND pacman -Syuu -q --needed --noconfirm
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm development mingw-w64-cross-toolchain libraries compression mingw-w64-cross sys-utils mingw-w64-clang-x86_64-toolchain mingw-w64-ucrt-x86_64-to
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm cmake gcc mingw-w64-cross-mingw64-gcc mingw-w64-cross-ucrt64-gcc mingw-w64-ucrt-x86_64-gcc mingw-w64-x86_64-gcc mingw-w64-clang-x86_64-cmake mi
    #)

#msys2-runtime-devel msys2-runtime  msys2-w32api-runtime

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm mingw-w64-x86_64-cmake
    #)

    #execute_process(
    #    COMMAND pacman -Sy -q --needed --noconfirm sqlite libsqlite cmake gcc gcc-libs pkgconf gnutls libgnutls liblzo2 mingw-w64-x86_64-pkgconf mingw-w64-x86_64-lzo2 mingw-w64-x86_64-gnutls min
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
		# Execute cygpath command
		execute_process(
			COMMAND bash -c "cygpath ${FOLDER_MODE} '${_folder_path}'" 
			OUTPUT_VARIABLE _return_path
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		# Return the result to the caller's scope
		set(${ret_var} ${_return_path} PARENT_SCOPE)
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
    template <> struct std::formatter<std::chrono::zoned_time<std::chrono::duration<long long int, std::ratio<86400>>, const std::chrono::time_zone*>, char> : std::formatter<std::string_view, cha>
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

endif()
