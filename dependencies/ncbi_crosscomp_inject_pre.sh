#!/bin/bash

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$(
    cd -- "$(dirname "$0")" >/dev/null 2>&1
    pwd -P
)"

FOLDER_MODE="-u" #for cygpath
NCBI_DIR=$(cygpath $FOLDER_MODE $(pwd))
SCRIPTPATH=$(cygpath $FOLDER_MODE $SCRIPTPATH)
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath $FOLDER_MODE $RTOOLS_DIR)

OPENMP_LIB_FLAG=""
if [[ ! -z $SHLIB_OPENMP_CXXFLAGS ]]; then
    OPENMP_LIB_FLAG="-lgomp.dll -lgomp"
    # SHLIB_OPENMP_CXXFLAGS="-D_OPENMP $SHLIB_OPENMP_CXXFLAGS"
fi

echo "PRE CONFIGURE INJECT..."
echo "BUILD DIR: $NCBI_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"
echo "OpenMP Flags : $SHLIB_OPENMP_CXXFLAGS"

cd $NCBI_DIR

export CC="/usr/bin/gcc"
export CXX="/usr/bin/g++"

COMPILER_VER=$($CC -dumpversion | tr -d '[:space:]')
COMPILER_VER_NOPUNCT=$($CC -dumpversion | tr -d '[:space:]' | tr -d '[:punct:]')

toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh GCC $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')

if [[ -z $toolchain_file ]]; then
    echo "EXITING. COULD NOT GUESS TOOLCHAIN"
    exit 0
fi

SHELL_TYPE=$(echo $MSYSTEM_CARCH | tr "[:upper:]" "[:lower:]")
echo "TOOLCHAIN FILE : $toolchain_file"
echo "CC : $CC"
echo "CXX : $CXX"
echo "COMPILER VER : $COMPILER_VER"

ln -sf $CC $(realpath $CC)-$COMPILER_VER
ln -sf $CXX $(realpath $CXX)-$COMPILER_VER

rm -f src/build-system/cmake/toolchains/$toolchain_file

touch src/build-system/cmake/toolchains/$toolchain_file

#############GCC
echo "set(NCBI_PTBCFG_FLAGS_DEFINED YES)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "include_guard(GLOBAL)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_C_FLAGS_INIT         \"-gdwarf-4 -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_FLAGS_INIT       \"-gdwarf-4 -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_SSE       \"-msse4.2\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_COVERAGE  \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_COVERAGE     \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_MAXDEBUG  \"-fsanitize=address -fstack-check\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_MAXDEBUG   \"-fsanitize=address\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_STATICCOMPONENTS \"-static-libgcc -static-libstdc++\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_HOST_SYSTEM_NAME cygwin)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SYSTEM_NAME msys)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_BUILD_TYPE Release)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(BUILD_SHARED_LIBS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_POSITION_INDEPENDENT_CODE ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER GNU)" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(NCBI_COMPILER_VERSION $COMPILER_VER_NOPUNCT)" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_C \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_CXX \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_C \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_C_COMPILER $CC)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_COMPILER $CXX)" >>src/build-system/cmake/toolchains/$toolchain_file

#############COMPILER - MSYS2 - GCC#################################
echo "set(CMAKE_C_FLAGS_RELEASE \" -std=gnu11 -O2 -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -fstack-protector-strong -fdiagnostics-color=always -pthread -mthreads -mtune=generic -maccumulate-outgoing-args -Wall -Wno-format-y2k \")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_CXX_FLAGS_RELEASE \" -std=gnu++17 -O2 -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -fstack-protector-strong -fdiagnostics-color=always $SHLIB_OPENMP_CXXFLAGS -pthread -mthreads -mtune=generic -maccumulate-outgoing-args -Wall -Wno-format-y2k \")" >>src/build-system/cmake/toolchains/$toolchain_file
###################################################

##############LINKER - GCC##############
echo "set(CMAKE_EXE_LINKER_FLAGS_INIT  \" -Wl,--demangle -Wl,--pic-executable -Wl,--support-old-code -Wl,--as-needed \")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_SHARED_LINKER_FLAGS_INIT  \" -Wl,--demangle -Wl,--pic-executable -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--dll -Wl,-shared \")" >>src/build-system/cmake/toolchains/$toolchain_file
################################

echo "include($NCBI_DIR/ncbi_components.cmake)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_EXTENSIONS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(MINGW32 1)" >>src/build-system/cmake/toolchains/$toolchain_file

##################################################################################################
