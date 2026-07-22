# $Id$

#Apache Arrow requires a minimum CMake version of 3.25
cmake_minimum_required(VERSION 3.25) 
cmake_policy(VERSION 3.25)

set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)
set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF)
set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF)

#set(ENV{CMAKE_EXPORT_NO_PACKAGE_REGISTRY} "1")
#set(ENV{XDG_DATA_HOME} "${CMAKE_CURRENT_BINARY_DIR}/tmp_home")
#set(ENV{HOME} "${CMAKE_CURRENT_BINARY_DIR}/tmp_home")

project(QuickBLAST)

#Assuming CMakeLists.txt is in tools/ and project root is in ../
set(PRJ_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")


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

# =========================================================================
# CRAN COMPLIANCE: SILENCE 3RD PARTY SOURCE WARNINGS
# =========================================================================
# CRAN scrapes the install log for warnings. Since these flags are inside 
# CMake and not Makevars, they bypass CRAN's static scanners.

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Silence One Definition Rule (LTO) warnings in NCBI C++ files
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-odr")
    
    # Silence stringop-overflow warnings in NCBI C files (blast_aalookup.c)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-stringop-overflow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
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
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD
     #   COMMAND_ERROR_IS_FATAL ANY    
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l 'apache-arrow-*.*.*' --sort=-v:refname | grep -iv \"dev\" | grep -iv \"rc\" | grep -iv \"js\" | head -n 1"
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow
        OUTPUT_VARIABLE ARROW_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    message("ARROW GIT LATEST TAG : " ${ARROW_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${ARROW_GIT_TAG}
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow    
        COMMAND_ERROR_IS_FATAL ANY
    )
    
    set(UV_BASE_DIR "${PRJ_ROOT_DIR}/src/BUILD/libuv/")
    
  	if (NOT EXISTS "${UV_BASE_DIR}")
        message(STATUS "Cloning libuv into ${UV_BASE_DIR}...")
        execute_process(
            COMMAND git clone https://github.com/libuv/libuv libuv
            WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/
            COMMAND_ERROR_IS_FATAL ANY
        )
    endif()
    
    # 3. Use CMake to get the absolute path securely (No bash required!)
    get_filename_component(UV_INC "${UV_BASE_DIR}/include" ABSOLUTE)
  	
  	#if(APPLE)
    #    # --- macOS (Apple Silicon & Intel Homebrew) ---
    #    find_path(UV_INC NAMES uv.h PATHS /opt/homebrew/include /usr/local/include)
    #    if(NOT UV_INC)
    #        message(FATAL_ERROR "Could not find libuv.h! Ensure brew install libuv ran.")
    #    endif()
    #    
    #    find_path(NGHTTP2_INC NAMES nghttp2/nghttp2.h PATHS /opt/homebrew/include /usr/local/include)
    #    if(NOT NGHTTP2_INC)
    #        message(FATAL_ERROR "Could not find nghttp2.h! Ensure brew install nghttp2 ran.")
    #    endif()
    #    
    #elseif(UNIX AND NOT APPLE)
    #    # --- Linux (Ubuntu/Debian via apt-get or CentOS via yum) ---
    #    # CMake will natively search /usr/include and /usr/local/include
    #    find_path(UV_INC NAMES uv.h)
    #    if(NOT UV_INC)
    #        message(FATAL_ERROR "Could not find libuv.h! Ensure 'libuv1-dev' is installed.")
    #    endif()
    #    
    #    find_path(NGHTTP2_INC NAMES nghttp2/nghttp2.h)
    #    if(NOT NGHTTP2_INC)
    #        message(FATAL_ERROR "Could not find nghttp2.h! Ensure 'libnghttp2-dev' is installed.")
    #    endif()
    #endif()
    #
    #message(STATUS "Resolved UV_INC: ${UV_INC}")
    #message(STATUS "Resolved NGHTTP2_INC: ${NGHTTP2_INC}")
  	
  	# Add Homebrew to search paths on macOS (Apple Silicon & Intel)
    if(APPLE)
        find_program(BREW_EXECUTABLE brew)
        if(BREW_EXECUTABLE)
            # Get the homebrew prefix (/opt/homebrew on M-series, /usr/local on Intel)
            execute_process(
                COMMAND ${BREW_EXECUTABLE} --prefix 
                OUTPUT_VARIABLE BREW_PREFIX 
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            
            # Point pkg-config to Homebrew's pkgconfig directory
            set(ENV{PKG_CONFIG_PATH} "${BREW_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
            
            # Point standard CMake find_package calls to Homebrew
            list(APPEND CMAKE_PREFIX_PATH "${BREW_PREFIX}")
        endif()
    endif()
  	
  	# --- Dependency Resolution via pkg-config ---
    find_package(PkgConfig QUIET)
    if(NOT PKG_CONFIG_FOUND)
        message(FATAL_ERROR 
            "----------------------------------------------------------\n"
            "Configuration failed: pkg-config is not installed.\n"
            "Please install pkg-config to build QuickBLAST.\n"
            "----------------------------------------------------------"
        )
    endif()
    
    # 1. Check for libuv
    pkg_check_modules(UV QUIET libuv)

    if(UV_FOUND)
        message(STATUS "Found system libuv: ${UV_INCLUDE_DIRS}")
    else()
        message(STATUS "System libuv not found via pkg-config. Falling back to downloaded libuv at ${UV_INC}...")
        
        # Verify the clone actually succeeded and the include directory exists
        if(EXISTS "${UV_INC}")
            set(UV_FOUND TRUE)
            set(UV_INCLUDE_DIRS "${UV_INC}")
        else()
            message(FATAL_ERROR 
                "----------------------------------------------------------\n"
                "Configuration failed: System libuv not found and downloaded fallback is missing!\n"
                "----------------------------------------------------------"
            )
        endif()
    endif()
    
    # 2. Check for nghttp2
    pkg_check_modules(NGHTTP2 QUIET libnghttp2)
    
    if(NGHTTP2_FOUND)
        message(STATUS "Found system nghttp2 via pkg-config: ${NGHTTP2_INCLUDE_DIRS}")
    else()
        message(STATUS "nghttp2 not found via pkg-config. Attempting manual search for macOS...")
        
        find_path(NGHTTP2_INCLUDE_DIRS 
            NAMES nghttp2/nghttp2.h
            HINTS 
            /opt/homebrew/include
            /usr/local/include
            /opt/R/arm64/include
            /opt/R/x86_64/include
            /opt/homebrew/opt/nghttp2/include   # Added Keg-only paths
            /usr/local/opt/nghttp2/include
        )
        
        find_library(NGHTTP2_LIBRARIES 
            NAMES nghttp2
            HINTS
            /opt/homebrew/lib
            /usr/local/lib
            /opt/R/arm64/lib
            /opt/R/x86_64/lib
            /opt/homebrew/opt/nghttp2/lib       # Added Keg-only paths
            /usr/local/opt/nghttp2/lib
        )
        
        if(NGHTTP2_INCLUDE_DIRS AND NGHTTP2_LIBRARIES)
            set(NGHTTP2_FOUND TRUE)
            message(STATUS "Found nghttp2 manually. Headers: ${NGHTTP2_INCLUDE_DIRS}, Libs: ${NGHTTP2_LIBRARIES}")
        else()
            message(FATAL_ERROR "Configuration failed: nghttp2 headers or libs were not found!")
        endif()
    endif()
    
    # Map the pkg-config outputs to your existing CMake variables so 
    # your execute_process() block downstream doesn't need to change.
    # pkg-config returns a list (even if it's just one path), so we extract the first item.
    list(GET UV_INCLUDE_DIRS 0 UV_INC)
    list(GET NGHTTP2_INCLUDE_DIRS 0 NGHTTP2_INC)
    
    message(STATUS "Resolved UV_INC via pkg-config: ${UV_INC}")
    message(STATUS "Resolved NGHTTP2_INC via pkg-config: ${NGHTTP2_INC}")
  	
    if (NOT EXISTS "${ROOT_SRC_DIR}/ncbi-cxx-toolkit-public/")
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/
    COMMAND_ERROR_IS_FATAL ANY
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l --sort=-v:refname | grep -iv \"dev\" | head -n 1" ; #no dev versions
        #COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r | head -n 1" ;
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public
        OUTPUT_VARIABLE NCBI_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("NCBI GIT LATEST TAG : " ${NCBI_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${NCBI_GIT_TAG}
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public 
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
    
    message("Copying " ${PRJ_ROOT_DIR}/src " into " ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/)

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/QuickBLAST.prj 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_crosscomp_inject_pre.sh 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_crosscomp_inject_post.sh 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_components.cmake 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(GLOB SRC_CONTENTS "${PRJ_ROOT_DIR}/src/*")
    foreach(ITEM IN LISTS SRC_CONTENTS)
        get_filename_component(FILENAME "${ITEM}" NAME)
        
        # Skip ONLY the root CMakeLists.txt so we don't overwrite NCBI's
        if(NOT FILENAME STREQUAL "CMakeLists.txt" AND NOT FILENAME STREQUAL "BUILD")
            file(COPY "${ITEM}" DESTINATION "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/src/")
        endif()
    endforeach()

    message("Copying " ${PRJ_ROOT_DIR}/src/include " into " ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/)

    file(COPY ${PRJ_ROOT_DIR}/src/include
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    execute_process(
        COMMAND chmod +x ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/cmake-configure ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/configure
    )

    set(CMAKE_CXX_EXTENSIONS OFF)

    execute_process(
    COMMAND mkdir -p ${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD
    WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow/
    )

    message("If build fails, change -DARROW_DEPENDENCY_SOURCE=AUTO TO -DARROW_DEPENDENCY_SOURCE=BUNDLED or install arrow with package manager (Default paths are seached for lib)")

if(APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin" OR CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    # Apple Clang (macOS) compatible flags
    set(ARROW_CXX_FLAGS " -w -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -pthread")
    set(ARROW_C_FLAGS  " -w -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -pthread")
else()
    # GCC compatible flags
    set(ARROW_CXX_FLAGS " -w -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -mtune=generic   -pthread")
    set(ARROW_C_FLAGS  " -w -fpermissive -pedantic -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -mtune=generic  -pthread")
endif()

set(ARROW_EXE_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")
set(ARROW_LINKER_FLAGS  "-Wl,-rpath,\$ORIGIN -Wl,-rpath,@loader_path/")


  execute_process(
    COMMAND cmake .. -G "Unix Makefiles" 
    # --- CRAN Policy Compliance Flags ---
    "-DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON"
    "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF"
    "-DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF"
    # --- MAC OS & DEPENDENCY FIXES ---
    "-DCMAKE_PREFIX_PATH=/opt/homebrew;/usr/local;/opt/R/arm64;/opt/R/x86_64"
    "-DARROW_DEPENDENCY_SOURCE=BUNDLED" 
    "-DBoost_SOURCE=SYSTEM"           # Skip building heavy Boost
    # ------------------------------------    
    "-DCMAKE_CXX_FLAGS=${ARROW_CXX_FLAGS}"
    "-DCMAKE_C_FLAGS=${ARROW_C_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${ARROW_EXE_FLAGS}"
    "-DCMAKE_SHARED_LINKER_FLAGS=${ARROW_LINKER_FLAGS}"
    "-DARROW_SIMD_LEVEL=NONE" 
    "-DARROW_RUNTIME_SIMD_LEVEL=NONE" 
    "-DARROW_JEMALLOC=OFF" 
    "-DARROW_MIMALLOC=OFF" 
    "-DARROW_COMPUTE=OFF" 
    "-DARROW_DATASET=OFF" 
    "-DARROW_JSON=OFF" 
    "-DARROW_WITH_RE2=OFF" 
    "-DARROW_WITH_UTF8PROC=OFF" 
    "-DARROW_BUILD_TESTS=OFF" 
    "-DARROW_DEPENDENCY_USE_SHARED=OFF"
    "-DARROW_BUILD_SHARED=OFF"
    "-DARROW_BUILD_STATIC=ON"
    "-DPARQUET_BUILD_SHARED=OFF"
    "-DPARQUET_BUILD_STATIC=ON"
    "-DARROW_PYTHON=OFF"
    "-DARROW_CSV=ON"
    "-DARROW_ENABLE_THREADING=ON"
    "-DARROW_FILESYSTEM=ON"
    "-DARROW_IPC=ON"
    "-DARROW_PARQUET=ON"
    "-DARROW_POSITION_INDEPENDENT_CODE=ON"
    "-DARROW_VERBOSE_THIRDPARTY_BUILD=ON"
    "-DARROW_WITH_BROTLI=ON"
    "-DARROW_WITH_BZ2=ON"
    "-DARROW_WITH_LZ4=ON"
    "-DARROW_RPATH_ORIGIN=ON"
    "-DARROW_WITH_SNAPPY=ON"
    "-DARROW_WITH_ZLIB=OFF"
    "-DARROW_WITH_ZSTD=ON"
    "-DCMAKE_BUILD_TYPE=Release"
    "-DBUILD_SHARED_LIBS=OFF"
    "-DBUILD_STATIC_LIBS=ON"
    "-DCMAKE_SKIP_BUILD_RPATH=FALSE"
    "-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"
    "-DCMAKE_INSTALL_RPATH=$ORIGIN"
    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
    
    # --- UNNECESSARY HEAVYWEIGHTS DISABLED FOR SPEED ---
    "-DARROW_TENSORFLOW=OFF"
    "-DARROW_HDFS=OFF"
    "-DARROW_WITH_BACKTRACE=OFF"

    # --- ENFORCE STATIC LINKING FOR CRAN ---
    "-DARROW_BROTLI_USE_SHARED=OFF"
    "-DARROW_LZ4_USE_SHARED=OFF"
    "-DARROW_ZSTD_USE_SHARED=OFF"
    "-DARROW_LLVM_USE_SHARED=OFF"
    "-DARROW_OPENSSL_USE_SHARED=OFF"
    "-DARROW_THRIFT_USE_SHARED=OFF"
    "-DARROW_UTF8PROC_USE_SHARED=OFF"
    "-DARROW_BOOST_USE_SHARED=OFF"
    "-DARROW_PROTOBUF_USE_SHARED=OFF"
    "-DARROW_SNAPPY_USE_SHARED=OFF"
    "-DARROW_BZ2_USE_SHARED=OFF"
    
    WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    )

#-DARROW_ZSTD_USE_SHARED=OFF -DARROW_JEMALLOC_USE_SHARED=OFF # -DARROW_CUDA=ON -DARROW_DEPENDENCY_SOURCE=BUNDLED  # -DARROW_BUILD_INTEGRATION=ON

    execute_process(
    COMMAND ${CMAKE_MAKE_PROGRAM} -s -j ${N_PROC} arrow_bundled_dependencies_merge parquet_static arrow
    WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD
    COMMAND_ERROR_IS_FATAL ANY
    RESULT_VARIABLE ARROW_BUILD_STATUS
    )

    #if(NOT ARROW_BUILD_STATUS EQUAL 0)
    #    message(STATUS "Parallel Arrow build failed or interrupted. Retrying with standard single-thread build...")
    #    execute_process(
    #        COMMAND ${CMAKE_MAKE_PROGRAM}
    #        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/arrow/cpp/BUILD"
    #        COMMAND_ERROR_IS_FATAL ANY
    #    )
    #endif()
    
    # ==============================================================================
    # SYNCHRONOUS SQLITE PRE-BUILD (Forces SQLite to exist before NCBI configures)
    # ==============================================================================
    set(SQLITE_ROOT_DIR "${PRJ_ROOT_DIR}/src/BUILD/sqlite")
    set(SQLITE_TARGET_LIB "${SQLITE_ROOT_DIR}/libsqlite3.a")
    
    if(NOT EXISTS "${SQLITE_TARGET_LIB}")
        message(STATUS "Pre-building SQLite3 statically during CMake configuration...")
        
        file(MAKE_DIRECTORY "${SQLITE_ROOT_DIR}")
        
        # 1. Download and Extract
        file(DOWNLOAD "https://sqlite.org/download.html" "${SQLITE_ROOT_DIR}/download.html" STATUS dl_status TIMEOUT 10)
        execute_process(COMMAND bash -c "grep -Eo '[0-9]{4}/sqlite-amalgamation-[0-9]+\\.zip' ${SQLITE_ROOT_DIR}/download.html | head -n 1" OUTPUT_VARIABLE SQLITE_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
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
        target_compile_options(sqlite3 PRIVATE -static-libgcc  )
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
    -DCMAKE_C_FLAGS=' -DPCRE2POSIX_EXP_DEFN= -w -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} ' \
    -DCMAKE_CXX_FLAGS=' -w -pedantic -Wformat -Werror=format-security -flto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} ' \
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
    -DCMAKE_C_FLAGS=' -w -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  ' \
    -DCMAKE_CXX_FLAGS=' -w -pedantic -Wformat -Werror=format-security -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} ' \
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
    
    #"CMAKE_PREFIX_PATH=${SQLITE_ROOT_DIR}"
    #"CMAKE_INCLUDE_PATH=${SQLITE_INC_DIR};${UV_INC}"
    #"CMAKE_LIBRARY_PATH=${SQLITE_ROOT_DIR}"
    
      execute_process(
        COMMAND ${CMAKE_COMMAND} -E env 
          "CC=${CMAKE_C_COMPILER}" 
          "CXX=${CMAKE_CXX_COMPILER}" 
          "CFLAGS=-DPCRE2POSIX_EXP_DEFN= -DPCRE2_STATIC=1"
          "CXXFLAGS=-DPCRE2POSIX_EXP_DEFN= -DPCRE2_STATIC=1"
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
          "-DNCBI_COMPONENT_UV_FOUND=YES"
          "-DNCBI_COMPONENT_UV_INCLUDE=${UV_INC}"
          "-DNCBI_COMPONENT_NGHTTP2_FOUND=YES"
          "-DNCBI_COMPONENT_NGHTTP2_INCLUDE=${NGHTTP2_INC}"
          "-DNCBI_COMPONENT_PSGLoader_DISABLED=YES"
            WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
            COMMAND_ERROR_IS_FATAL ANY
            OUTPUT_FILE "${PRJ_ROOT_DIR}/src/BUILD/ncbi_cmake_configure.log"
            ERROR_FILE "${PRJ_ROOT_DIR}/src/BUILD/ncbi_cmake_configure.log"
            ERROR_QUIET
        )
    #endif()
    execute_process(
          COMMAND ${CMAKE_MAKE_PROGRAM} -s -j ${N_PROC} datatool QuickBLASTcpp makeblastdb 
          WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/BUILD/build
          COMMAND_ERROR_IS_FATAL ANY
          RESULT_VARIABLE NCBI_BUILD_STATUS
        )
      
    #if(NOT NCBI_BUILD_STATUS EQUAL 0)
    #    message(STATUS "NCBI build failed or interrupted. Retrying with standard single-thread build...")
    #    execute_process(
    #        COMMAND ${CMAKE_MAKE_PROGRAM} datatool QuickBLASTcpp makeblastdb 
    #        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ncbi-cxx-toolkit-public/BUILD/build"
    #        COMMAND_ERROR_IS_FATAL ANY
    #    )
    #endif()

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

    set(ROOT_SRC_DIR "${PRJ_ROOT_DIR}/BUILD/")
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
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD
        #COMMAND_ERROR_IS_FATAL ANY    
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l 'apache-arrow-*.*.*' --sort=-v:refname | grep -iv \"dev\" | grep -iv \"rc\" | grep -iv \"js\" | head -n 1"
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow
        OUTPUT_VARIABLE ARROW_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    message("ARROW GIT LATEST TAG : " ${ARROW_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${ARROW_GIT_TAG}
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/arrow    
        COMMAND_ERROR_IS_FATAL ANY
    )

    if (NOT EXISTS "${ROOT_SRC_DIR}/ncbi-cxx-toolkit-public/")
    execute_process(
    COMMAND git clone https://github.com/ncbi/ncbi-cxx-toolkit-public ncbi-cxx-toolkit-public
    WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD
    #COMMAND_ERROR_IS_FATAL ANY
    )
    endif()

    execute_process(
        COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l --sort=-v:refname | grep -iv \"dev\" | head -n 1" ; #no dev versions
        #COMMAND ${RTOOLS_DIR}/usr/bin/bash -c "git tag -l | sort -r | head -n 1" ;
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public
        OUTPUT_VARIABLE NCBI_GIT_TAG
    	OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message("NCBI GIT LATEST TAG : " ${NCBI_GIT_TAG})

    execute_process(
        COMMAND git checkout -q ${NCBI_GIT_TAG}
        WORKING_DIRECTORY ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public 
        COMMAND_ERROR_IS_FATAL ANY
    )

    message("Copying " ${PRJ_ROOT_DIR}/src/ " into " ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/)

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/QuickBLAST.prj 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_crosscomp_inject_pre.sh 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_crosscomp_inject_post.sh 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(COPY ${PRJ_ROOT_DIR}/src/dependencies/ncbi_components.cmake 
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    file(GLOB SRC_CONTENTS "${PRJ_ROOT_DIR}/src/*")
    foreach(ITEM IN LISTS SRC_CONTENTS)
        get_filename_component(FILENAME "${ITEM}" NAME)
        
        # Skip ONLY the root CMakeLists.txt so we don't overwrite NCBI's
        if(NOT FILENAME STREQUAL "CMakeLists.txt" AND NOT FILENAME STREQUAL "BUILD")
            file(COPY "${ITEM}" DESTINATION "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/src/")
        endif()
    endforeach()

    message("Copying " ${PRJ_ROOT_DIR}/src/include " into " ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/)

    file(COPY ${PRJ_ROOT_DIR}/src/include
    DESTINATION ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/
    )

    execute_process(
        COMMAND chmod +x ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/cmake-configure ${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/configure
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
    file(WRITE "${PRJ_ROOT_DIR}/src/BUILD/patch_arrow.cmake" [=[
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
        PREFIX "${PRJ_ROOT_DIR}/src/BUILD/arrow_ep"
        SOURCE_DIR "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp"
        BINARY_DIR "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2"
        
        PATCH_COMMAND ${CMAKE_COMMAND} -DARROW_FILE=${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/src/arrow/compute/kernels/scalar_temporal_unary.cc -P "${PRJ_ROOT_DIR}/src/BUILD/patch_arrow.cmake"
        
        CMAKE_GENERATOR "MSYS Makefiles"
        
        # Pass all your definitions here
        CMAKE_ARGS
            -DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON
            -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF
            -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DARROW_BUILD_INTEGRATION=OFF
            -DARROW_SIMD_LEVEL=NONE 
          	-DARROW_RUNTIME_SIMD_LEVEL=NONE
          	-DARROW_DEPENDENCY_SOURCE=BUNDLED
          	-DBoost_SOURCE=SYSTEM          
          	-DARROW_JEMALLOC=OFF 
          	-DARROW_MIMALLOC=OFF 
          	-DARROW_COMPUTE=OFF 
          	-DARROW_DATASET=OFF 
          	-DARROW_JSON=OFF 
          	-DARROW_WITH_RE2=OFF 
          	-DARROW_WITH_UTF8PROC=OFF 
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
#            -DPARQUET_MINIMAL_DEPENDENCY=ON
            -DARROW_WITH_BROTLI=ON
#            -DARROW_WITH_RE2=OFF
            -DARROW_WITH_ZLIB=ON
            -DARROW_WITH_ZSTD=ON
            -DARROW_IPC=ON
            -DARROW_FILESYSTEM=ON
            -DARROW_PYTHON=OFF
            -DARROW_ENABLE_THREADING=ON
            -DARROW_ACERO=OFF
            -DARROW_PARQUET=ON
            #-DARROW_WITH_THRIFT=ON
            #-DThrift_SOURCE=BUNDLED
            #-DTHRIFT_BUILD_COMPILER=OFF
            -DARROW_CSV=ON
            -DARROW_HDFS=ON
            -DARROW_VERBOSE_THIRDPARTY_BUILD=ON
            -DARROW_POSITION_INDEPENDENT_CODE=ON
#            -DARROW_JSON=ON
            -DARROW_SUBSTRAIT=ON
            -DARROW_WITH_SNAPPY=ON
#            -DARROW_WITH_UTF8PROC=ON
            -DARROW_RPATH_ORIGIN=ON
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
            -DCMAKE_C_COMPILER=gcc.exe
            -DCMAKE_CXX_COMPILER=g++.exe
            "-DCMAKE_C_FLAGS= -pedantic -Wformat -Werror=format-security -UHAVE_AF_UNIX_H -D__MINGW32__ -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -mtune=generic  -pthread   "
            "-DCMAKE_CXX_FLAGS= -pedantic -Wformat -Werror=format-security -UHAVE_AF_UNIX_H -D__MINGW32__ -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS} -mtune=generic   -pthread  "
            "-DCMAKE_SHARED_LINKER_FLAGS=-Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,-shared -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,--dll -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--no-whole-archive -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group "

        # Map the Makevars build command
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -s -j ${N_PROC} arrow_bundled_dependencies_merge parquet_static arrow
        
        # Map the Makevars copy/install command
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow.a"
            "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow_msys2.a"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release/libparquet.a"
            "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release/libparquet_msys2.a"
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

    if(NOT EXISTS "${PRJ_ROOT_DIR}/src/BUILD/sqlite/libsqlite3.a")    
        file(DOWNLOAD
            "https://sqlite.org/download.html"
            "${PRJ_ROOT_DIR}/src/BUILD/download.html"
            SHOW_PROGRESS
            STATUS download_status
            TIMEOUT 10
        )
        list(GET download_status 0 download_code)
        if(NOT "${download_code}" STREQUAL "0")
            message(FATAL_ERROR "SQLite version check failed: ${download_status}")
        endif()
        execute_process(
            COMMAND bash -c "grep -Eo '[0-9]{4}/sqlite-amalgamation-[0-9]+\\\\.zip' ${PRJ_ROOT_DIR}/src/BUILD/download.html | head -n 1"
            OUTPUT_VARIABLE SQLITE_LINK
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if (SQLITE_LINK)
            message(STATUS "SQLite version link found: https://sqlite.org/${SQLITE_LINK}")
        else()
            message(FATAL_ERROR "SQLite version link not found")
        endif()
        execute_process(
            COMMAND bash -c "mkdir -p ${PRJ_ROOT_DIR}/src/BUILD/sqlite/"
            )
        file(DOWNLOAD
            "https://sqlite.org/${SQLITE_LINK}"
            "${PRJ_ROOT_DIR}/src/BUILD/sqlite/sqlite.zip"
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
        -DROOT_SRC_DIR=${PRJ_ROOT_DIR} \
        -DR_INCLUDE_DIR='${R_INCLUDE_DIR}' \
        -DNCBITK_ROOT_DIR=${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public \
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
        -DCMAKE_C_FLAGS='${CMAKE_C_FLAGS} -pedantic -Wformat -Werror=format-security -Wa,-mbig-obj -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -mtune=generic  -pthread  -w ' \
        -DCMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS} -DLIBICONV_STATIC -pedantic -Wformat -Werror=format-security -Wa,-mbig-obj -D_FORTIFY_SOURCE=2 -ffat-lto-objects -flto=auto -fstack-protector-strong ${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}  -mtune=generic   -pthread  -w ' \
        -DCMAKE_SHARED_LINKER_FLAGS=' -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,-shared -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,--dll -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--no-whole-archive -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group ' \
        -DCMAKE_EXE_LINKER_FLAGS=' -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--pic-executable -Wl,--subsystem,console -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 -lsqlite3 -liconv -lintl -lbcrypt -Wl,--end-group ' \
        -DCMAKE_USE_WIN32_THREADS_INIT=1 \
        -DCMAKE_USE_PTHREADS_INIT=0 \
        -DTHREADS_PREFER_PTHREAD_FLAG=OFF \
    ")
    #-static-libstdc++ -static-libgcc
    #-ftrack-macro-expansion=0 -pipe 
    #-std=gnu2x
    #-std=gnu++20
    #-O3 -fPIC

    ExternalProject_Add(ncbi_toolkit
        PREFIX "${PRJ_ROOT_DIR}/src/BUILD/ncbi_ep"
        SOURCE_DIR "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public"
        BUILD_IN_SOURCE TRUE
        
        # STEP 1: Run pre-inject and cmake-configure
        CONFIGURE_COMMAND bash -c "SHLIB_OPENMP_CXXFLAGS='${SHLIB_OPENMP_CXXFLAGS} ${TOOLKIT_CPPFLAGS}' RTOOLS_DIR='${RTOOLS_DIR}' R_INCLUDE_DIR='${R_INCLUDE_DIR}' ./ncbi_crosscomp_inject_pre.sh && export CMAKE_ARGS=\"${NCBI_CMAKE_ARGS}\" && export CMAKE='${RTOOLS_DIR}/ucrt64/bin/cmake.exe' && ./cmake-configure GCC ${GCC_VER} --with-projects=QuickBLAST.prj --without-dll --with-composite --without-debug --with-build-root=BUILD"
        #--with-dll
        # STEP 2: Use double quotes, escape inner quotes, and use $$ to pass command substitution safely to Bash
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -C BUILD/build -j ${N_PROC} datatool all makeblastdb QuickBLASTcpp
        
        # STEP 3: Run the post-injection script
        INSTALL_COMMAND bash -c "echo ncbi-cxx-toolkit build complete!"
        #INSTALL_COMMAND bash -c "RTOOLS_DIR='${RTOOLS_DIR}' ./ncbi_crosscomp_inject_post.sh"
        
        #OUTPUT_FILE "${PRJ_ROOT_DIR}/src/BUILD/ncbi_cmake_configure.log"
        #ERROR_FILE "${PRJ_ROOT_DIR}/src/BUILD/ncbi_cmake_configure.log"
        #LOG_CONFIGURE ON
        #LOG_MERGED_STDOUTERR ON
    )

    execute_process(COMMAND Rscript -e "cat(Sys.getenv('R_HOME'))" OUTPUT_VARIABLE R_USR_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(R_INC "${R_USR_HOME}/include")
    cygpath_cmd(R_INC "${R_USR_HOME}/include")
    message("R inc:" ${R_INC})
    
    execute_process(COMMAND Rscript -e "cat(file.path(Sys.getenv('R_HOME'),'bin',Sys.getenv('R_ARCH')))" OUTPUT_VARIABLE R_BIN OUTPUT_STRIP_TRAILING_WHITESPACE)
    cygpath_cmd(R_BIN "${R_BIN}")
    message("R BIN:" ${R_BIN})

    # Ensure Arrow builds before NCBI
    add_dependencies(ncbi_toolkit build_arrow)
    #add_dependencies(ncbi_toolkit_static ncbi_toolkit)

    # =========================================================================
    # 3. DEFINE THE MAIN QUICKBLAST TARGET
    # =========================================================================

    # Gather your main Rcpp / C++ source files
    file(GLOB QUICKBLAST_SOURCES "${PRJ_ROOT_DIR}/src/*.cpp" "${PRJ_ROOT_DIR}/src/*.c")
    
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
        IMPORTED_LOCATION "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release/libarrow_msys2.a"
    )
    # Tell CMake: "Do not try to link Arrow until build_arrow is 100% finished"
    add_dependencies(Arrow::msys2 build_arrow)

    # =========================================================================
    # 5. LINK NCBI (IMPORTED TARGETS)
    # =========================================================================

    # Example: Linking the core NCBI static libraries (update names based on what NCBI outputs)
    add_library(NCBI::core STATIC IMPORTED GLOBAL)
    set_target_properties(NCBI::core PROPERTIES
        IMPORTED_LOCATION "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/BUILD/lib/libncbi_core.a"
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

  # if(WIN32 OR MINGW)
  #    # Prevent MinGW from exceeding the 64k DLL export limit by hiding static lib symbols
  #    # Tell CMake explicitly NOT to auto-export everything
  #    set_target_properties(QuickBLAST PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF)
  #    # Tell the MinGW linker not to re-export symbols from all the NCBI static archives
  #    target_link_options(QuickBLAST PRIVATE "-Wl,--exclude-libs=ALL")
  #    # Force RcppExport to act as a Windows DLLExport macro for this target
  #    target_compile_definitions(QuickBLAST PRIVATE 
  #       # RcppExport=extern\ \"C\"\ __declspec(dllexport)
  #        "RcppExport=extern \"C\" __declspec(dllexport)"
  #    )
  # endif()
  
  if(WIN32 OR MINGW)
    # 1. Keep this to ensure static archive symbols stay private
    target_link_options(QuickBLAST PRIVATE "-Wl,--exclude-libs=ALL")
    # 2. Automatically generate a .def file exporting !~~ONLY~~ the R (side) functions
    set(DEF_FILE "${PRJ_ROOT_DIR}/src/BUILD/QuickBLAST.def")
    file(WRITE "${DEF_FILE}" 
    "EXPORTS
        R_init_QuickBLAST
        _QuickBLAST_isQuickBLASTLoaded
        _QuickBLAST_RecordBatchVectorToFlattenedDFList
        _QuickBLAST_CreateQuickBLASTInstance
        _QuickBLAST_GetInstanceCount
        _QuickBLAST_GetInstanceID
        _QuickBLAST_GetQuickBLASTInstance
        _QuickBLAST_DeleteQuickBLASTInstance
        _QuickBLAST_GetQuickBLASTOptions
        _QuickBLAST_SetQuickBLASTOptions
        _QuickBLAST_BLAST2Seqs
        _QuickBLAST_BLAST2Folders
        _QuickBLAST_BLAST1Folder
        _QuickBLAST_BLAST2Files
        _QuickBLAST_RemoteBLAST
        _QuickBLAST_isBLASTDB
        _QuickBLAST_MakeBLASTDB
        _QuickBLAST_BLAST2DBs
        _QuickBLAST_GetFASTAHeaders
        _QuickBLAST_BLASTFile2DB
    ")
    # 3. Add the .def file directly to the target sources
    # The Windows linker will read this and open the door exclusively for R!
    target_sources(QuickBLAST PRIVATE "${DEF_FILE}")
    # NOTE: Completely REMOVE the target_compile_definitions(QuickBLAST PRIVATE "RcppExport=...") line!
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
        "${R_INCLUDE_DIR}"
        "${PRJ_ROOT_DIR}/src/include/"
        "${RCPP_FOLDER}/include"
        "${RCPPPROG_FOLDER}/include"
        "${RCPPTHREAD_FOLDER}/include"
        "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/src"
        "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/src"
        "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/include" # To find QuickBLAST/commons.hpp
    )
    
    target_compile_options(QuickBLAST PRIVATE
        "-include" "rcpp_workaround.hpp"
        ${SHLIB_OPENMP_CXXFLAGS}
    )
    
    # 1. Tell CMake where the directory of the built libraries is
    target_link_directories(QuickBLAST PRIVATE
        "${PRJ_ROOT_DIR}/src/BUILD/sqlite"
        "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/BUILD/bin"
        "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/BUILD/lib"
        "${PRJ_ROOT_DIR}/src/BUILD/lib"
        "${PRJ_ROOT_DIR}/src/BUILD/bin"
        "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release"
        "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/lib"
        "${R_BIN}"
    )

    set(ICONV_STATIC_LIB "${RTOOLS_DIR}/x86_64-w64-mingw32.static.posix/lib/libiconv.a")

    # 2. Link the libraries by name (do not include "lib" or ".a")
    target_link_libraries(QuickBLAST PUBLIC
        ${PKG_LIBS}
        -Wl,--enable-auto-import
        -Wl,--start-group
        bz2
        arrow
        arrow_bundled_dependencies
        parquet
        thrift
        "${PRJ_ROOT_DIR}/src/BUILD/sqlite/libsqlite3.a"
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
        #R
        "${R_BIN}/R.dll"
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
    if(NOT DEFINED R_ARCH)
        if(DEFINED ENV{R_ARCH})
            set(R_ARCH "$ENV{R_ARCH}")
        else()
            # Dynamically discover and ask the R executable for its active sub-architecture layout
            find_program(R_EXECUTABLE NAMES R R.exe)
            if(R_EXECUTABLE)
                execute_process(
                    COMMAND ${R_EXECUTABLE} CMD config R_ARCH
                    OUTPUT_VARIABLE R_ARCH_OUT
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                )
                set(R_ARCH "${R_ARCH_OUT}")
            else()
                set(R_ARCH "")
            endif()
        endif()
    endif()
    
    # ---------------------------------------------------------------------
    # Sanitize MSYS2 path conversion on Windows
    # If MSYS2 translated "/x64" into "C:/rtools45/x64", strip it back down.
    # ---------------------------------------------------------------------
    if(R_ARCH MATCHES "(x64|i386)$")
        set(R_ARCH "/${CMAKE_MATCH_1}")
    endif()
    
    # Destination inside your R package source
    set(DEV_LIB_DEST "${PRJ_ROOT_DIR}/libs${R_ARCH}")
    
    # ---------------------------------------------------------------------
    # Dynamically target R's temporary staging folder during installation
    # ---------------------------------------------------------------------
    if(DEFINED ENV{R_PACKAGE_DIR})
        # R CMD INSTALL is running - Staged installation
        file(TO_CMAKE_PATH "$ENV{R_PACKAGE_DIR}" R_PKG_DIR)
        set(INST_LIB_DEST "${R_PKG_DIR}/libs${R_ARCH}")
        set(INST_BIN_DEST "${R_PKG_DIR}/bin${R_ARCH}")
    else()
        # Fallback for local development (devtools::load_all)
        set(INST_LIB_DEST "${PRJ_ROOT_DIR}/inst/libs${R_ARCH}")
        set(INST_BIN_DEST "${PRJ_ROOT_DIR}/inst/bin${R_ARCH}")
    endif()

    add_custom_command(TARGET QuickBLAST POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INST_LIB_DEST}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${INST_BIN_DEST}"
        COMMENT "Syncing QuickBLAST and NCBI libraries to R package directories..."
    )

    # Source paths (adjust these if your BUILD directories differ)
    set(NCBI_BIN_DIR "${PRJ_ROOT_DIR}/src/BUILD/ncbi-cxx-toolkit-public/BUILD/bin")
    set(ARROW_BIN_DIR "${PRJ_ROOT_DIR}/src/BUILD/arrow/cpp/BUILD_MSYS2/release")

    ## 2. List of specific NCBI files to strip and copy
    set(NCBI_FILES
        "makeblastdb.exe"
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
    )

    foreach(FILENAME ${NCBI_FILES})
        
        # Dynamically route Executables to bin/ and DLLs to libs/
        if(FILENAME MATCHES "\\.exe$" OR FILENAME MATCHES "makeblastdb")
            set(FINAL_DEST "${INST_BIN_DEST}")
            set(DEV_DEST "${PRJ_ROOT_DIR}/inst/bin${R_ARCH}")
        else()
            set(FINAL_DEST "${INST_LIB_DEST}")
            set(DEV_DEST "${DEV_LIB_DEST}") # Assuming this maps to src/ or inst/libs for devtools
        endif()

        if(DEFINED ENV{R_PACKAGE_DIR})
            # STRICT CRAN MODE: Only copy to the staging directory to prevent source pollution
            add_custom_command(TARGET QuickBLAST POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                        "${NCBI_BIN_DIR}/${FILENAME}" "${FINAL_DEST}/${FILENAME}"
            )
        else()
            # LOCAL DEV MODE: Copy to devtools directories as well
            add_custom_command(TARGET QuickBLAST POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory "${DEV_DEST}"
                
                COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                        "${NCBI_BIN_DIR}/${FILENAME}" "${DEV_DEST}/${FILENAME}"
                        
                COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                        "${NCBI_BIN_DIR}/${FILENAME}" "${FINAL_DEST}/${FILENAME}"
            )
        endif()

    endforeach()

    ## Process Arrow DLLs
    ## Note: Since CMake's copy command doesn't support wildcards directly,
    ## we use the system shell provided by Rtools/MinGW
    #add_custom_command(TARGET QuickBLAST POST_BUILD
    #    COMMAND bash -c "cp -uf ${ARROW_BIN_DIR}/*.dll ${INST_LIBS_DEST}/"
    #    COMMENT "Deploying Arrow dependencies"
    #)

endif()
